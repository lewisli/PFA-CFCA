

%results_path = '/home/lewisli/Documents/Research/Direct_Forecasting/Prior_Case_2/';
NumRealizations = 500;
addpath('../util');

load('../../data/Case2/PriorModels/Prior.mat');
CaseName = 'NoiseLevel500';
%% Generate Data Structs
ForecastColumn = 4;
HistoricalColumn = 4;

% This is the time step that we divide the forecast and history
TimeColumn = 2;
NumTimeSteps = 240;
EndTime = 12000;
ForecastObjName={'PNEW'};
HistoricalObjName = {'P1','P2','P3','P4','P5'};
ForecastSpline = [4 20];
HistoricalSpline = [4 20];
ForecastStart = 160;
HistoricalEnd = 80;

[HistoricalStruct, ForecastStruct] = GenerateDataStructsWithInterpolation(Data, ...
    PropertyNames, ForecastColumn, HistoricalColumn, TimeColumn, ...
    HistoricalEnd,ForecastStart, NumTimeSteps, ForecastSpline, ...
    HistoricalSpline, ForecastObjName, HistoricalObjName, EndTime);

ForecastStruct.time = linspace(3500,7500,length(ForecastStruct.time));

%%
% Set aside one realization that we will deem the "reference";
TruthRealization = 120;

close all;
% Plot to verify data structures/choice of input/output
SaveFolder = ['../../figures/' CaseName '/'];
FontSize = 22;

h1  = PlotInputResponse( HistoricalStruct,TruthRealization,FontSize,SaveFolder);
h2  = PlotInputResponse( ForecastStruct,TruthRealization,FontSize,SaveFolder);


%% Picking the basis functions
% This is a trial-error process, where you will pick the order and number
% of knots to use as the splines. Generally the longer the time forecast,
% the more knots will be required. The best way to do this is to pick some
% splines, and graphically view the resulting fit for some randomly selected
% realizations. This is implemented in ComputeHarmonicScores
addpath('../cfca');
HistoricalStruct.spline=[3 40]; % 6th order B-Spline with 20 knots
histPCA = ComputeHarmonicScores(HistoricalStruct,4,SaveFolder);

%ForecastStruct.spline = [6 20]; % 6th order B-Spline with 20 knots
predPCA = ComputeHarmonicScores(ForecastStruct,3,SaveFolder);

%%
% Perform CFCAa
% The eigenvalue tolerance sets the number of eigenvalues we will keep
% after FPCA. Keeping too many eigenvalues may need to highly oscillating
% posterior times series; while keeping too little results in oversmoothed
% and unrealistic models
EigenvalueTolerance = 0.9;
OutlierPercentile = 95;
epsilon = 500;

% Run CFCA: The results are the mean/covariance of h*(in Gaussian space)
% and Dc,Hc are the prior model coefficients in the canonical space
PlotLevel = 1;
FontSize = 24;
[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,PlotLevel,FontSize,SaveFolder,epsilon);

%% Sample from CFCA posterior and transform forecasts back into time domain
close all;
NumPosteriorSamples = 1000;

[SampledPosteriorRealizations,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA,0,0,0,SaveFolder);

% Compute quantiles
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    ForecastStruct.data, SampledPosteriorRealizations);
close all;
% Plot sampled responses and quantiles
PlotPosteriorSamplesAndQuantiles(ForecastStruct,TruthRealization, ...
    SampledPosteriorRealizations,PriorQuantiles,PosteriorQuantiles,SaveFolder);

display(['Average Posterior Distance: ' num2str(mean(PosteriorQuantiles(3,:) - ...
    PosteriorQuantiles(1,:)))]);
