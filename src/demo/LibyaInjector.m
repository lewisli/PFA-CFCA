% LibyaInjector.m
% 
% Forecast injector efficiency using production data
project_dir = 'C:\Users\Lewis Li\ResearchData\LibyanCase\Situation1\results\'
NumRealizations = 500;
WellNames = {'P1','P2','P3','P4','P5','I1','I2','I3'};

load([project_dir 'SimulationResults.mat']);

%% Generate Data Structures for Production Rates
addpath('../util');
ForecastColumn = 4;
HistoricalColumn = 4;

% This is the time step that we divide the forecast and history
TimeColumn = 2;
NumTimeSteps = 200;
EndTime = 4000;
ForecastObjName={'I2'};
HistoricalObjName = {'P1','P2','P3','P4','P5'};
ForecastSpline = [4 20];
HistoricalSpline = [4 20];
ForecastStart = 160;
HistoricalEnd = 80;

[HistoricalStruct, ~] = GenerateDataStructsWithInterpolation(Data, ...
    PropertyNames, ForecastColumn, HistoricalColumn, TimeColumn, ...
    HistoricalEnd,ForecastStart, NumTimeSteps, ForecastSpline, ...
    HistoricalSpline, ForecastObjName, HistoricalObjName, EndTime);
%% Plot production rates
h1  = PlotInputResponse( HistoricalStruct,22,22);

%% Read in forecast


ForecastStruct = struct();
ForecastStruct.time = HistoricalStruct.time;
ForecastStruct.ObjNames =  {'I2'};
ForecastStruct.name= 'Efficency';
ForecastStruct.type = 'Forecast';
ForecastStruct.RunNames = HistoricalStruct.RunNames;
ForecastStruct.spline = [4,20];

NumRuns = length(ForecastStruct.RunNames);
NumTimeSteps = length(ForecastStruct.time);
NumInjectors = length(ForecastStruct.ObjNames);

ForecastStruct.data = zeros(NumRuns,NumTimeSteps,NumInjectors);

for i = 1:NumRuns
     RunName = ['Run' num2str(HistoricalStruct.RunNames(i))];
     RunWAF = load([project_dir 'waf/' RunName '.mat']);
     
     for j = 1:NumInjectors
         RealWAFTime = RunWAF.(RunName).(ForecastStruct.ObjNames{j}).('time');
         RealWAF = RunWAF.(RunName).(ForecastStruct.ObjNames{j}).('inj_eff');    
         WAF = interp1(RealWAFTime,RealWAF,HistoricalStruct.time,'pchip',0);
         ForecastStruct.data(i,:,j) = WAF;
     end
end
%% Plot effiency rates
h1  = PlotInputResponse( ForecastStruct,22,22);

%% Picking the basis functions
% This is a trial-error process, where you will pick the order and number
% of knots to use as the splines. Generally the longer the time forecast,
% the more knots will be required. The best way to do this is to pick some
% splines, and graphically view the resulting fit for some randomly selected
% realizations. This is implemented in ComputeHarmonicScores
addpath('../cfca');
close all;
HistoricalStruct.spline=[3 24]; % 6th order B-Spline with 20 knots
histPCA = ComputeHarmonicScores(HistoricalStruct,0);

ForecastStruct.spline = [2 30]; % 6th order B-Spline with 20 knots
predPCA = ComputeHarmonicScores(ForecastStruct,4);


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
TruthRealization = 1;

%%


[ mu_posterior, C_posterior, Dc, Df, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,PlotLevel,FontSize,0,0);

%% Sample from CFCA posterior and transform forecasts back into time domain
NumPosteriorSamples = 1000;

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
