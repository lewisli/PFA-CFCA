function [ C_f ] = EstimateFunctionalErrorCovariance( HistoricalStruct, EigenTolerance,epsilon)
%EstimateFunctionalErrorCovariance Uses a Monte Carlo approach to compute
%   Detailed explanation goes here
addpath('../cfca');
% Step 1: Extract

NummRealizations = size(HistoricalStruct.data,1);
NumTimeSteps = size(HistoricalStruct.data,2);
NumHistoricalResponses = size(HistoricalStruct.data,3);


noiseless_scores = ComputeMixedPCAScores(HistoricalStruct,EigenTolerance);

% 
df_error = zeros(size(noiseless_scores));

for i = 1:NummRealizations
     fprintf('Working on realization %d\n',i)
     NoisyStruct = HistoricalStruct;
     
     NoisyStruct.data(i,:,:) = HistoricalStruct.data(i,:,:) + ...
         randn(1,NumTimeSteps,NumHistoricalResponses)*epsilon;
     
     noisy_score = ComputeMixedPCAScores(NoisyStruct,EigenTolerance);
     df_error(i,:) = noisy_score(i,:) - noiseless_scores(i,:);

end
 
C_f = cov(df_error);

end

