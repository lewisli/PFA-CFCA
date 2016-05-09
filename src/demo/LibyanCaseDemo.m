close all; clear all;
addpath('../cfca');
addpath('../util');
load('../../data/PriorRuns/Prior.mat');

%% Step 1b) Generate data structures
% We call everything before this step the 'observed' data
HistoricalEnd = 65;

% We call everything after this step the 'forecast'
ForecastStart = 125;

% Total number of days simulated in 3DSL
TotalDaysSimulated=11500;

% Total number of time steps
TotalNumTimeSteps = 200;

% The column in the Data struct that refers to the attribute we want to use
% as the forecast/historical
ForecastColumn = 4;   % Oil Rate
HistoricalColumn = 4; % Oil Rate
TimeColumn = 2;       % Simulation Time

% Object on which forecasting is required (New well to be drilled)
ForecastObjectName = {'PNEW2'};

% Existing wells whose production rates are used as historical data
HistoricalObjectName = {'P1','P2','P3','P4','P5'};

% Generates data structure which will be used later for CFCA
[HistoricalStruct,ForecastStruct] = GenerateDataStructsWithInterpolation(Data,...
    PropertyNames,ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,...
    ForecastStart,TotalNumTimeSteps,[6 20],[6 20],ForecastObjectName,...
    HistoricalObjectName,TotalDaysSimulated);

%%
% Set aside one realization that we will deem the "reference";
TruthRealization = 12;

% Plot to verify data structures/choice of input/output
%h1  = PlotInputResponse( HistoricalStruct,TruthRealization);
h2  = PlotInputResponse( ForecastStruct,TruthRealization);

%% Step 2: Picking the basis splines for functional data analysis
% This is a trial-error process, where you will pick the order and number
% of knots to use as the splines. Generally the longer the time forecast,
% the more knots will be required. The best way to do this is to pick some
% splines, and graphically view the resulting fit for some randomly selected
% realizations. This is implemented in ComputeHarmonicScores

close all;
addpath('cfca');
HistoricalStruct.spline=[6 20]; % Use for rates
histPCA = ComputeHarmonicScores(HistoricalStruct,4);

ForecastStruct.time = linspace(HistoricalStruct.time(end),...
    HistoricalStruct.time(end)+4000,length(ForecastStruct.time));
ForecastStruct.spline = [6 20];
predPCA = ComputeHarmonicScores(ForecastStruct,3);

%%
% The eigenvalue tolerance sets the number of eigenvalues we will keep
% after FPCA. Keeping too many eigenvalues may need to highly oscillating
% posterior times series; while keeping too little results in oversmoothed
% and unrealistic models
%EigenvalueTolerance = 0.995;
EigenvalueTolerance = 0.99;
OutlierPercentile = 95;

% Run CFCA: The results are the mean/covariance of h*(in Gaussian space)
% and Dc,Hc are the prior model coefficients in the canonical space
[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,1);

%% Step 4: Once we have the posterior means and covariances, we
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

%% Step 5: Compare with Rejection Sampling
load('../../data/RejectionSampling/RJ.mat');
[PriorQuantiles, RJQuantiles] = ComputeQuantiles(ForecastStruct.data,RJStruct.data);

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
HistoricalObjectName = {'P5'};

% Generates data structure which will be used later for CFCA
[HistoricalStruct,ForecastStruct] = GenerateDataStructsWithInterpolation(Data,...
    PropertyNames,ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,...
    ForecastStart,TotalNumTimeSteps,[6 20],[6 20],ForecastObjectName,...
    HistoricalObjectName,TotalDaysSimulated);
ForecastStruct.time = linspace(HistoricalStruct.time(end),...
    HistoricalStruct.time(end)+4000,length(ForecastStruct.time));
h1  = PlotInputResponse( HistoricalStruct,TruthRealization);

[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,1);

[SampledPosteriorRealizations,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);

% Compute quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    ForecastStruct.data, SampledPosteriorRealizations);

% Plot sampled responses and quantiles
PlotPosteriorSamplesAndQuantiles(ForecastStruct,TruthRealization, ...
    SampledPosteriorRealizations,PriorQuantiles,PosteriorQuantiles);

%% Boot strap
addpath('../bootstrap');
addpath('../cfca');
NumBootstrap = 1000;
BootstrapSize = 500;

[omega, d_hat,d_empirical] = FunctionalBootstrap(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    NumBootstrap, BootstrapSize);
