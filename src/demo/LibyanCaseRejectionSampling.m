%% LibyanCaseRejectionSampling.m
%
% Demonstration of PFA-CFCA on responses obtained from analog of the
% WinstersHall Concession C97-I. Comparison with 5000 rejection sampling
% runs
%
% Author: Lewis Li
% Date: September 29th 2016

close all; clear all;
addpath('../cfca');
addpath('../util');

CaseName = 'RejectionSampling'
SaveFolder = ['../../figures/' CaseName '/'];

%% Load data structures
load('../../data/PriorRuns/PriorData.mat');

RS_Path = '/home/lewisli/code-dev/directforecasting/data/Compressible/RejectionSample/RejectionSample.mat';
load(RS_Path);

%%

% Set aside one realization that we will deem the "reference";
TruthRealization = 12;
FontSize = 22;

% Plot to verify data structures/choice of input/output
h1  = PlotInputResponse( HistoricalStruct,TruthRealization,FontSize);
h2  = PlotInputResponse( ForecastStruct,TruthRealization,FontSize);



%%
addpath('../util');
ForecastObjectName = {'PNEW2'};

% Best Case. Average Difference = 238.234
HistoricalObjectName = {'P1','P2','P3','P4','P5'};
% This is the time step that we divide the forecast and history
HistoricalEnd = 65;
ForecastStart = 125;
TotalNumTimeSteps = 200;

% The column in the Data struct that refers to the attribute we want to use
% as the forecast/historical
ForecastColumn = 4;
HistoricalColumn = 4;
TimeColumn = 2;

% Generates data structure which will be used later for CFCA
[HistoricalRJ,ForecastRJ] = GenerateDataStructsWithInterpolation(Data,PropertyNames,...
    ForecastColumn,HistoricalColumn,TimeColumn,HistoricalEnd,ForecastStart,...
    TotalNumTimeSteps,[4 8],[4 8],ForecastObjectName,HistoricalObjectName,11500);


%%



%%
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


% Perform CFCA
% The eigenvalue tolerance sets the number of eigenvalues we will keep
% after FPCA. Keeping too many eigenvalues may need to highly oscillating
% posterior times series; while keeping too little results in oversmoothed
% and unrealistic models
close all;
CaseName= 'Case_1/RJComparison';
SaveFolder = ['../../figures/' CaseName '/'];
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end


EigenvalueTolerance = 0.99;
OutlierPercentile = 95;
Epilson = 50;
MeasurementNoise= Epilson;
% Run CFCA: The results are the mean/covariance of h*(in Gaussian space)
% and Dc,Hc are the prior model coefficients in the canonical space
PlotLevel = 1;
FontSize = 24;
[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,PlotLevel,FontSize,SaveFolder,MeasurementNoise);

%%
% Sample from CFCA posterior and transform forecasts back into time domain
NumPosteriorSamples = 100;
close all;
[SampledPosteriorRealizations,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA,0,0,0,SaveFolder);

% Compute quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    ForecastStruct.data, SampledPosteriorRealizations);

% Plot sampled responses and quantiles
PlotPosteriorSamplesAndQuantiles(ForecastStruct,TruthRealization, ...
    SampledPosteriorRealizations,PriorQuantiles,PosteriorQuantiles,SaveFolder);

display(['Average Posterior Distance: ' num2str(mean(PosteriorQuantiles(3,:) - ...
    PosteriorQuantiles(1,:)))]);

display(['Average Prior Distance: ' num2str(mean(PriorQuantiles(3,:) - ...
    PriorQuantiles(1,:)))]);



%
quantiles = [.1,.5,.9];

WellID = 1;
ObservedData = HistoricalStruct.data(TruthRealization,:,WellID);
TrueForecast = ForecastStruct.data(TruthRealization,:);

Epilson = 400;
NumTimeSteps = size(HistoricalStruct.data,2);
errorMat = (diag(ones(NumTimeSteps,1))*Epilson);
errorMatInv = inv(errorMat);

NumRJModels = size(HistoricalRJ.data,1);
likelihood = zeros(NumRJModels,1);
for i = 1:NumRJModels


    diff = HistoricalRJ.data(i,:,WellID) - ObservedData;    
    likelihood(i) = exp(-0.5*diff*errorMatInv*diff'/(NumTimeSteps^2));
end

% Scaled likelihood
[MaxLikelihood,MaxID] = max(likelihood);
KeptModels = find((rand(NumRJModels,1) <= likelihood./MaxLikelihood)==1) ;


RJQuantiles = quantile(ForecastRJ.data(KeptModels,:,:),quantiles);
FigHandle = figure('Position', [100, 100, 600, 600]);
hold on;
h0 = plot(ForecastStruct.time,PriorQuantiles','color',[0.5 0.5 0.5],...
    'LineWidth',3);
h2 = plot(ForecastStruct.time,ForecastStruct.data(TruthRealization,:),...
    'r','LineWidth',3);
h3 = plot(ForecastStruct.time,PosteriorQuantiles,'b--','LineWidth',3);
h4 = plot(ForecastStruct.time,RJQuantiles,'k--','LineWidth',3);

legend([h0(1), h2(1),h3(1),h4(1)],'Prior','Reference','Posterior','Rejection Sampling');
xlabel('t(days)');ylabel(['Forecasted: ' ForecastStruct.name]);axis square;
title('Quantiles');
set(gca,'FontSize',18);
set(gcf,'color','w');
axis tight;

display(['Average Posterior Distance: ' num2str(mean(PosteriorQuantiles(3,:) - ...
    PosteriorQuantiles(1,:)))]);

display(['Average Rejection Sampling Distance: ' num2str(mean(RJQuantiles(3,:) - ...
    RJQuantiles(1,:)))]);

display(['Average Prior Distance: ' num2str(mean(PriorQuantiles(3,:) - ...
    PriorQuantiles(1,:)))]);

export_fig -m2 RJ_Sample1.png




