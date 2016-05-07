function [ omega,d_hat ] = BootstrapSignifcance(HistoricalStruct,...
    ForecastStruct,TruthRealization,EigenvalueTolerance,OutlierPercentile,...
    NumBootstrap, BootstrapSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath('../util');

% Compute empirical difference measure
[ mu_posterior, C_posterior, Dc, Hc, Hf, B, dobs_c] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile);

ReferenceForecastFirstStep = ForecastStruct.data(TruthRealization,1);
NumPosteriorSamples = 100;
RJTolerance = 0.001;

CurveValidity = 0;  % Use when either is fine
predPCA = ComputeHarmonicScores(ForecastStruct,0);

% Generate posterior samples
[SampledPosteriorRealizations,~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA,ReferenceForecastFirstStep,RJTolerance,1);

close all;
d_empirical = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));

d_hat = zeros(NumBootstrap,1);

h = waitbar(0,'Running bootstrap...');

for b = 1:NumBootstrap
    
    
    close all;
    
    BootstrapIndex = randsample(length(HistoricalStruct.data),...
        BootstrapSize);
    
    
    HistoricalStructBootstrap = HistoricalStruct;
    HistoricalStructBootstrap.data = HistoricalStructBootstrap.data(...
        BootstrapIndex,:,:);
    HistoricalStructBootstrap.data = [HistoricalStructBootstrap.data;...
        HistoricalStruct.data(TruthRealization,:,:)];
    
    ForecastStructBootstrap = ForecastStruct;
    ForecastStructBootstrap.data = ForecastStructBootstrap.data(...
        BootstrapIndex,:,:);
    ForecastStructBootstrap.data = [ForecastStructBootstrap.data;...
        ForecastStruct.data(TruthRealization,:,:)];
    
    BootstrapTruthIndex = BootstrapSize + 1;
    
    [ mu_posterior, C_posterior, Dc, Hc, Hf, B, dobs_c] = ...
        ComputeCFCAPosterior(HistoricalStructBootstrap, ...
        ForecastStructBootstrap, BootstrapTruthIndex, EigenvalueTolerance,...
        OutlierPercentile,0);
    
    [~,~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA,ReferenceForecastFirstStep,RJTolerance,1);
    
    d_hat(b) = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));
    
    waitbar(b / NumBootstrap);

end
close(h);


omega = sum(d_hat<d_empirical)/NumBootstrap;

end

