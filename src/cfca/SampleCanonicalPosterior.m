function [ h_output,hf_out] = SampleCanonicalPosterior( ...
    mu_posterior, C_posterior, NumEstimates,Hc,B,Hf,Time_Forecast,predPCA,...
    ReferenceForecastFirstStep, RJTolerance,DirectionalValidity)
%SampleCanonicalPosterior Generate samples of the canonical posterior (that
%is sampled from f(h|d_obs), and reconstructs the time series
%   Detailed explanation goes here
%
% Outputs:
%   h_output: Posterior samples in time domain
%   hc_output:  Posterior samples in canonical coefficients
%   HpostCoef: Posterior samples in functional coefficients

if (nargin < 9)
    ReferenceForecastFirstStep=0;
    DirectionalValidity = 0;
    RJTolerance = 100;
end
addpath('../../thirdparty/fda_matlab');
h_output = zeros(NumEstimates,length(Time_Forecast));
hc_out = zeros(NumEstimates,length(mu_posterior));
hf_out = zeros(NumEstimates,length(mu_posterior));
NumValid = 1;

rng('shuffle');

CycleSize = 10000;

NumCycles = 0;

if (all(eig(C_posterior) >= 0))
    C_posterior = (C_posterior + C_posterior')/2;
    %display(C_posterior)
    %eig(C_posterior)
    
    while (NumValid < NumEstimates)
        % The posterior distribution is a normal distribution in canonical space
        PosteriorSamples = mvnrnd(mu_posterior',C_posterior,CycleSize)';
        
        % The h_c are not Gaussian, hence we need to backtransform
        PosteriorSamplesTransformed = BackTransform(PosteriorSamples,Hc);
        
        % We need to undo the canonical transform...
        HpostCoef = PosteriorSamplesTransformed'*pinv(B)+repmat(mean(Hf,1)',...
            1,CycleSize)'; % H_f_posterior
        
        % Finally, we reconstruct the time series (mean_FDA + sum(HpostCoef*PhiH))
        numPredCoeffs = size(Hf,2);
        
        % Principal components for H
        PhiH = eval_fd(Time_Forecast,predPCA{1}.harmfd);
        
        % Re-construct time series
        h_reconstructed = repmat(eval_fd(Time_Forecast,predPCA{1}.meanfd),...
            1,CycleSize) + PhiH(:,1:numPredCoeffs)*HpostCoef(:,1:numPredCoeffs)';
        
        close all;
        
        % Compute difference between first forecasted time and observed value
        % at that time
        Difference = abs(repmat(ReferenceForecastFirstStep,1,CycleSize)-...
            h_reconstructed(1,:))./ReferenceForecastFirstStep;
        
        % A model is valid if the starting point is within a tolerance of the
        % observed and if there
        AllowUpOnly = 1;
        AllowDownOnly = -1;
        
        DirectionalValidityFlag = ones(CycleSize,1);
        
        if (DirectionalValidity == AllowUpOnly)
            for i = 1:CycleSize
                DirectionalValidityFlag(i) = (max(diff(h_reconstructed(:,i))) < 0);
            end
        elseif (DirectionalValidity == AllowDownOnly)
            for i = 1:CycleSize
                DirectionalValidityFlag(i) = (min(diff(h_reconstructed(:,i))) > 0);
            end
        elseif (DirectionalValidity == 0)
            for i = 1:CycleSize
                DirectionalValidityFlag(i) = 1;
            end
        end
        
        %     CycleValid = logical(Difference<RJTolerance) & ...
        %         logical(DirectionalValidityFlag');
        
        CycleValid = logical(DirectionalValidityFlag');
        
        NumCycleValid = sum(CycleValid);
        
        % This means we need another cycle of model sampling
        if (NumValid+NumCycleValid < NumEstimates+1)
            h_output(NumValid:NumValid+NumCycleValid-1,:) = ...
                h_reconstructed(:,CycleValid)';
            
            hc_out(NumValid:NumValid+NumCycleValid-1,:) = ...
                PosteriorSamplesTransformed(:,CycleValid)';
            
            hf_out(NumValid:NumValid+NumCycleValid-1,:) = ...
                HpostCoef(CycleValid,:);
            
            NumValid = NumValid+NumCycleValid;
            NumCycles = NumCycles + 1;
        else
            NumCycles = NumCycles + 1;
            RemainingModels = NumEstimates - NumValid + 1;
            %display(['Finished after ' num2str(NumCycles) ' cycles.']);
            
            ValidModels = h_reconstructed(:,CycleValid)';
            ValidHc = PosteriorSamplesTransformed(:,CycleValid)';
            ValidHf = HpostCoef(CycleValid,:);
            
            h_output(NumValid:end,:) = ValidModels(1:RemainingModels,:);
            hc_out(NumValid:end,:) = ValidHc(1:RemainingModels,:);
            
            NumOut = size(hf_out(NumValid:end,:),2);
            hf_out(NumValid:end,:) = ValidHf(1:RemainingModels,1:NumOut);
            
            return
        end
        
        display(['After ' num2str(NumCycles) ' cycles. Rejection sampler has found '...
            num2str(NumValid) ' models']);
    end
end
end

