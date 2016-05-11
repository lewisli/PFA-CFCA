%% LibyanCaseDemo.m
%
% Demonstration of PFA-CFCA on responses obtained from analog of the 
% WinstersHall Concession C97-I.
%
% Author: Lewis Li
% Date: May 10th 2016

close all; clear all;
addpath('../cfca');
addpath('../util');

%% Load data structures
load('../../data/PriorRuns/PriorData.mat');

% Set aside one realization that we will deem the "reference";
TruthRealization = 12;
FontSize = 32;

% Plot to verify data structures/choice of input/output
h1  = PlotInputResponse( HistoricalStruct,TruthRealization,FontSize);
h2  = PlotInputResponse( ForecastStruct,TruthRealization,FontSize);

%% Picking the basis functions
% This is a trial-error process, where you will pick the order and number
% of knots to use as the splines. Generally the longer the time forecast,
% the more knots will be required. The best way to do this is to pick some
% splines, and graphically view the resulting fit for some randomly selected
% realizations. This is implemented in ComputeHarmonicScores
addpath('cfca');
HistoricalStruct.spline=[6 20]; % 6th order B-Spline with 20 knots
histPCA = ComputeHarmonicScores(HistoricalStruct,4);

ForecastStruct.spline = [6 20]; % 6th order B-Spline with 20 knots
predPCA = ComputeHarmonicScores(ForecastStruct,3);


%% Perform CFCA
% The eigenvalue tolerance sets the number of eigenvalues we will keep
% after FPCA. Keeping too many eigenvalues may need to highly oscillating
% posterior times series; while keeping too little results in oversmoothed
% and unrealistic models
EigenvalueTolerance = 0.99;
OutlierPercentile = 95;

% Run CFCA: The results are the mean/covariance of h*(in Gaussian space)
% and Dc,Hc are the prior model coefficients in the canonical space
PlotLevel = 1;
FontSize = 24;
[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,PlotLevel,FontSize);

%% Sample from CFCA posterior and transform forecasts back into time domain
NumPosteriorSamples = 100;

[SampledPosteriorRealizations,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);

% Compute quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    ForecastStruct.data, SampledPosteriorRealizations);

% Plot sampled responses and quantiles
PlotPosteriorSamplesAndQuantiles(ForecastStruct,TruthRealization, ...
    SampledPosteriorRealizations,PriorQuantiles,PosteriorQuantiles);

display(['Average Posterior Distance: ' num2str(mean(PosteriorQuantiles(3,:) - ...
    PosteriorQuantiles(1,:)))]);


%% Boot strap to estimate confidence
addpath('../bootstrap');
addpath('../cfca');
NumBootstrap = 1000;
BootstrapSize = 500;

[omega, d_hat,d_empirical] = FunctionalBootstrap(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    NumBootstrap, BootstrapSize);

%% Compare with Rejection Sampling
load('../../data/RejectionSampling/RJ.mat');
[PriorQuantiles, RJQuantiles] = ComputeQuantiles(ForecastStruct.data,...
    RJStruct.data);

FontSize=32;
hold on;
h1 = plot(ForecastStruct.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(ForecastStruct.time,RJQuantiles,'k:','LineWidth',3);
h3 = plot(ForecastStruct.time,PosteriorQuantiles,'b--','LineWidth',3);

legend([h1(1), h2(1),h3(1)],'Prior','Direct Forecasting','Rejection Sampling');
xlabel('t(days)');ylabel(['Forecasted: ' ForecastStruct.name]);axis square;
title('Quantiles');
set(gca,'FontSize',FontSize);
xlim([round(ForecastStruct.time(1)) round(ForecastStruct.time(end)+1)])
set(gcf,'color','w');

%% Forecast just using P5
load('../../data/PriorRuns/PriorP5Data.mat');
h1  = PlotInputResponse( HistoricalStruct,TruthRealization,FontSize);
h2  = PlotInputResponse( ForecastStruct,TruthRealization,FontSize);
EigenvalueTolerance = 0.95;
OutlierPercentile = 100;
predPCA = ComputeHarmonicScores(ForecastStruct,3);

[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ...
    ComputeCFCAPosterior(HistoricalStruct, ForecastStruct, ...
    TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,1,FontSize);

NumPosteriorSamples = 100;
[SampledPosteriorRealizations,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);

% Compute quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    ForecastStruct.data, SampledPosteriorRealizations);

% Plot sampled responses and quantiles
PlotPosteriorSamplesAndQuantiles(ForecastStruct,TruthRealization, ...
    SampledPosteriorRealizations,PriorQuantiles,PosteriorQuantiles);

%% Boot strap for just P5
addpath('../bootstrap');
addpath('../cfca');
NumBootstrap = 1000;
BootstrapSize = 500;

[omega, d_hat,d_empirical] = FunctionalBootstrap(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    NumBootstrap, BootstrapSize);

%% Plot bootstrap results
addpath('../../data/BootstrapResults/');
load('BootstrapTest.mat');

OmegaValues = cell2mat(BootstrapResults(:,2));
PosteriorUncertainty = cell2mat(BootstrapResults(:,3));

scatter(OmegaValues,PosteriorUncertainty,60,'filled');
set(gcf,'color','w');
xlabel('P90-P10 stb oil/day','FontSize',FontSize);
ylabel('\omega','FontSize',FontSize);
set(gca,'FontSize',FontSize);

text(OmegaValues(5)*0.95,PosteriorUncertainty(5)*1.2,...
    BootstrapResults{5,1},'FontSize',FontSize);
text(OmegaValues(4)*1.01,PosteriorUncertainty(4)*1.01,...
    BootstrapResults{4,1},'FontSize',FontSize);
text(OmegaValues(end)*1.01,PosteriorUncertainty(end)*1.01,...
    'P1,P2,P3,P4,P5','FontSize',FontSize);
text(OmegaValues(17)*1.01,PosteriorUncertainty(17)*1.01,...
    'P1,P2,P4','FontSize',FontSize);
axis tight;


