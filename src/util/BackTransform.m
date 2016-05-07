% function Output = BackTransform(Input, Match)
% Author: Lewis Li (lewisli@stanford.edu)
% Date:    Feburary 7th 2016
%
% Function   : BackTransform
%
% Description:  Computes a histogram transform on Gaussian variable Input such 
%               that the resulting variable given in Output has a histogram that
%               matches that of the variable given in Match.
%
% Parameters :
% 	Input 	   		Gaussian variable to be transformed
%	Match:	   		Variable's who's histogram we are trying to match
%
% Return     :
%	Output:	The transformed variable

function [ Output ] = BackTransform(Input, Match )

Output = zeros(size(Input));

for i = 1:size(Input,1)
    figure
    OriginalScores = Match(:,i);
    TransformedSamples = Input(i,:);
    BackTransformedValue = zeros(size(TransformedSamples));
    
    [F,x] = ecdf(OriginalScores);
    
    % Get the cdf values at each value assuming a gaussian with mean/std
    % given by OriginalScores
    FStar = normcdf(TransformedSamples,mean(OriginalScores),std(OriginalScores));
    
    % for each FStar, find closest F
    for j = 1:length(FStar)
        
        [c index] = min(abs(F-FStar(j)));
        
        if (index == 1)
            BackTransformedValue(j) = x(index);
        elseif (index == length(x))
            BackTransformedValue(j) = x(end);
        elseif (F(index) < FStar(j))
            BackTransformedValue(j) = 0.5*(x(index)+x(index-1));
        else
            BackTransformedValue(j) = 0.5*(x(index)+x(index+1));
        end
    end
    
    subplot(1,3,1)
    hist(BackTransformedValue);
    title(['BackTransformed for h^c_' num2str(i)])
    axis square;
    subplot(1,3,2)
    hist(OriginalScores);
    title(['Original for h^c_' num2str(i)])
    axis square;
    subplot(1,3,3)
    qqplot(OriginalScores,BackTransformedValue)
    axis square;
    
    Output(i,:) = BackTransformedValue;
end

end

