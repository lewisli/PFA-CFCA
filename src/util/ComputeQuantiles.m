function [ PriorQuantiles,PosteriorQuantiles ] = ComputeQuantiles(Prior,...
    Posterior )
%ComputeQuantiles Summary of this function goes here
%   Detailed explanation goes here

quantiles = [.1,.5,.9];

PriorQuantiles = quantile(Prior,quantiles);
PosteriorQuantiles = quantile(Posterior,quantiles);

end

