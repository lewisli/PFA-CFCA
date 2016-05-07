function [ omega,d_hat,d_empirical ] = FunctionalBootstrap(HistoricalStruct,...
    ForecastStruct,TruthRealization,EigenvalueTolerance,OutlierPercentile,...
    NumBootstrap,RJTolerance, BootstrapSize)
%FunctionalBootstrap Bootstrap on the functional components
%   Detailed explanation goes here
addpath('../thirdparty/fda_matlab');
% Compute empirical difference measure
[ mu_posterior, C_posterior, ~, ~,Hc, Hf, B, ~] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    OutlierPercentile,0);

ReferenceForecastFirstStep = ForecastStruct.data(TruthRealization,1);
NumPosteriorSamples = BootstrapSize;

CurveValidity = 0;  % Use when either is fine
predPCA = ComputeHarmonicScores(ForecastStruct,0);

% Generate posterior samples
[~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);

close all;
d_empirical = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));
d_hat = zeros(NumBootstrap,1);

h = waitbar(0,'Running bootstrap...');

for b = 1:NumBootstrap
    %close gcf;
    BootstrapIndex = randsample(length(HistoricalStruct.data),...
        BootstrapSize,true);
    
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
    
    [ mu_posterior, C_posterior, ~, ~,Hc, Hf, B, ~] = ...
        ComputeCFCAPosterior(HistoricalStructBootstrap, ...
        ForecastStructBootstrap, BootstrapTruthIndex, EigenvalueTolerance,...
        OutlierPercentile,0);
    
    [~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA,ReferenceForecastFirstStep,RJTolerance,...
    CurveValidity);
    
    d_hat(b) = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));
    
    waitbar(b / NumBootstrap);
    
%     figure(b)
%     hold on;
% [fpost,xipost] = ksdensity(Hf_post(:,1));
% [fprior,xiprior] = ksdensity(Hf(:,1));
% plot(xipost,fpost,'b','LineWidth',3)
% plot(xiprior,fprior,'r','LineWidth',3)
% hlegend = legend('Posterior','Prior');
% set(hlegend,'FontSize',15);
% set(hlegend,'Location','best');
% savefig([num2str(b) '.fig']);

end
close(h);


omega = sum(d_hat<d_empirical)/NumBootstrap;

end

