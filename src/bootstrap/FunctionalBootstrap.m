function [ omega,d_hat,d_empirical ] = FunctionalBootstrap(HistoricalStruct,...
    ForecastStruct,TruthRealization,EigenvalueTolerance,NumBootstrap,...
    BootstrapSize)
%FunctionalBootstrap Bootstrap on the functional components
%   Perform a bootstrap test of confidence on the posterior quantiles
%
% Inputs:
%   HistoricalStruct: Struct containing historical data
%   ForecastStruct: Struct containing forecast data
%   TruthRealization: Realization corresponding to truth or reference
%   EigenvalueTolerance: Determines the number of dimensions after PFCA we
%   will use
%   NumBoostrap: Number of times we will perform the bootstrap
%   BootstrapSize: Number of elements we will bootstrap each time
%
% Author: Lewis Li (lewisli@stanford.edu, contact@lewisli.me)
% Date: May 9th 2016

addpath('../thirdparty/fda_matlab');

% Compute empirical difference measure
[ mu_posterior, C_posterior, ~, ~,Hc, Hf, B, ~] = ComputeCFCAPosterior(...
    HistoricalStruct, ForecastStruct, TruthRealization, EigenvalueTolerance,...
    100,0);

% Generate posterior samples for empirical case
NumPosteriorSamples = BootstrapSize;
predPCA = ComputeHarmonicScores(ForecastStruct,0);
[~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);
d_empirical = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));
d_hat = zeros(NumBootstrap,1);
close all;


%h = waitbar(0,'Running bootstrap...');


for b = 1:NumBootstrap
    
    % Sample with replacement from prior set of models
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
    
    % Put truth back into bootstrap samples
    BootstrapTruthIndex = BootstrapSize + 1;
    
    % Compute posterior
    [ mu_posterior, C_posterior, ~, ~,Hc, Hf, B, ~] = ...
        ComputeCFCAPosterior(HistoricalStructBootstrap, ...
        ForecastStructBootstrap, BootstrapTruthIndex, EigenvalueTolerance,...
        100,0);
    
    % Generate posterior samples
    [~,Hf_post]= SampleCanonicalPosterior(...
    mu_posterior,C_posterior,NumPosteriorSamples,Hc,B,Hf,...
    ForecastStruct.time,predPCA);
    
    % Compute bootstrap distances
    d_hat(b) = ComputeCDFDistance(Hf_post(:,1),Hf(:,1));
    
    %waitbar(b / NumBootstrap);
end

%close(h);


omega = sum(d_hat<d_empirical)/NumBootstrap;

end

