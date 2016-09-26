results_path = '/home/lewisli/Documents/Research/Direct_Forecasting/Prior_Case_2/';
NumRealizations = 500;
WellNames = {'P1','P2','P3','P4','P5','PNEW'};

[Data, PropertyNames] = Process3DSLResults(results_path,NumRealizations,WellNames);

%% Generate Data Structs
ForecastColumn = 4;
HistoricalColumn = 4;

% This is the time step that we divide the forecast and history
TimeColumn = 2;
NumTimeSteps = 130;
EndTime = 12000;
ForecastObjName={'PNEW'};
HistoricalObjName = {'P1','P2','P3','P4','P5'};
ForecastSpline = [4 20];
HistoricalSpline = [4 20];
ForecastStart = 90;
HistoricalEnd = 90;

[HistoricalStruct, ForecastStruct] = GenerateDataStructsWithInterpolation(Data, ...
    PropertyNames, ForecastColumn, HistoricalColumn, TimeColumn, ...
    HistoricalEnd,ForecastStart, NumTimeSteps, ForecastSpline, ...
    HistoricalSpline, ForecastObjName, HistoricalObjName, EndTime);

%%
% Set aside one realization that we will deem the "reference";
TruthRealization = 12;
FontSize = 20;
close all;
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
