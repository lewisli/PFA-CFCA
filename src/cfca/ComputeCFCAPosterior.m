function [ mu_posterior, C_posterior,Dc,Df,Hc,Hf, B, dobs_c] = ...
    ComputeCFCAPosterior(HistoricalStruct, ForecastStruct, TruthRealization, ...
    EigenTolerance,OutlierPercentile,PlotLevel,FontSize,SavePath,epsilon)
%ComputeCFCAPosterior Computes posterior distribution of forecasts
%conditioned to d_obs
%   Performs CFCA to compute posterior mean and covariance in canonical
%   space.
%
% Inputs:
%   HistoricalStruct: Struct containing historical responses
%   ForecastStrct: Struct containing forecast responses
%   TruthRealization: Realization number corresponding to the truth
%   EigenTolerance: % of variance we will use to pick how many eigenvalues
%   are kept
%   OutlierPercentile: % of models we will keep if we perform outlier
%   detection (default is to keep everything 100)
%   PlotLevel: Plot canonical and functional results
%   FontSize: Font size in plots
%
% Outputs:
%   mu_posterior: posterior mean
%   C_posterior: posterior covariance
%   Dc: historical in canonical coordinates
%   Df: historical in functional coordinates
%   Hc: forecast in canonical coordinates
%   Hf: forecast in functional coordinates
%   B: tranformation matrix from CCA
%   dobs_c: observed data in canonical coordinates
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 5th 2016

if (nargin < 6)
    PlotLevel = 1;
    FontSize=12;
end

if (nargin < 8)
    SaveOn = false;
else
    SaveOn = true;
end


Responses_Historical = HistoricalStruct.data;
NumHistoricalResponses = size(Responses_Historical,3);

AvailableRealizations = setdiff(1:size(Responses_Historical,1),...
    TruthRealization);

PlotLevelPCA = 0;
histPCA = ComputeHarmonicScores(HistoricalStruct,PlotLevelPCA);
predPCA = ComputeHarmonicScores(ForecastStruct,PlotLevelPCA);


%% Plot Eigenvalue Functions
MinEigenValues = 2;

% Get number of required harmonics required for forecasts
nHarmPred = GetNumHarmonics(predPCA{1}, MinEigenValues,EigenTolerance);

% Get number of required harmonics required for forecasts
nHarmHist = 0;
for i = 1:NumHistoricalResponses
    responseNumHarms = GetNumHarmonics(histPCA{i},...
        MinEigenValues,EigenTolerance);
    nHarmHist = max(nHarmHist,responseNumHarms);
end

% Generate matrix of forecast harmonic scores
harmscrpred=predPCA{1}.harmscr;
Hf = harmscrpred(AvailableRealizations,1:nHarmPred);

% Generate matrix of historical harmonic scores
Df = zeros(length(AvailableRealizations),nHarmHist*NumHistoricalResponses);
dobs_f = zeros(1,nHarmHist*NumHistoricalResponses);

% Iterate over each historical response (ex: P1, P2)
for i = 1:NumHistoricalResponses
    harmscrhist=histPCA{i}.harmscr;
    
    % Need to re-arrange harmonic scores into Df such that the first
    % eigenvalues are placed in first
    for j = 1:nHarmHist
        Df(:,(i-1)*nHarmHist + j) = harmscrhist(AvailableRealizations,j);
        dobs_f(:,(i-1)*nHarmHist + j) = harmscrhist(TruthRealization,j);
    end
end

% Run PCA On Df
% Poor naming choice by FDA package... we need to remove it from the path
% in order to run PCA
rmpath('../../thirdparty/fda_matlab');

DfStar = [Df; dobs_f];
[~,score,~] = pca(DfStar);

% Project dobs_f onto PCA basis
dobs_fpca = score(end,:);
score(end,:) = [];

% Perform CCA
Df = score;
[A, B, ~, Dc,Hc] = canoncorr(score,Hf);

% Project dobs_f into canonical space
dobs_c=(dobs_fpca-mean(score))*A;

% Apply an outlier detection
[Dc,Hc] = LSOutlierDetection(Dc,Hc,OutlierPercentile);

if (PlotLevel == 1)
    PlotLowDimModels(Dc,Hc,dobs_c,'c',FontSize);
    if SaveOn == true
        export_fig([SavePath 'Dc_Hc'], '-png','-m3');
    end
    PlotLowDimModels(Df,Hf,dobs_fpca,'f',FontSize);
    if SaveOn == true
        export_fig([SavePath 'Df_Hf'], '-png','-m3');
    end
end

% Perform a normal score transform
Hc_gauss = NormalScoreTransform(Hc,0);
C_H = cov(Hc_gauss);
H_CG_Mean = mean(Hc_gauss)';

% Find best linear bit between Dc and Hc_gauss
G = Dc'/Hc_gauss';
DDiff= Dc'-G*Hc_gauss';
C_T = DDiff*DDiff'/length(Dc);


if epsilon ==0
    C_Dc = zeros(size(C_T));
else
    C_Df = EstimateFunctionalErrorCovariance( HistoricalStruct,EigenTolerance,epsilon);
    C_Dc = A'*C_Df*A;
end


% Perform Gaussian Regression
mu_posterior = H_CG_Mean + C_H*G'*pinv(G*C_H*G' + C_T+C_Dc)*(dobs_c'-G*H_CG_Mean);
C_posterior = C_H - C_H*G'*inv(G*C_H*G' + C_T+C_Dc)*G*C_H;
%C_posterior = inv(G'*pinv(C_T)*G + inv(C_H));

end

