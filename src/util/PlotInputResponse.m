function [ h ] = PlotInputResponse( DataStruct,TruthRealization,FontSize)
%PlotInputResponse: Plots input response for historical and forecasts
%
% Inputs:
%   HistoricalStruct: Struct containing historical data
%   ForecastStruct: Struct containing forecast data
%   TruthRealization: Index of realization taken to be the truth
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: March 4th 2016

if (nargin < 4)
    FontSize=12;
end

NumRealizations = size(DataStruct.data,1);
NumResponses = size(DataStruct.data,3);

h = figure('Units', 'normalized', 'Position', [0,0,1,1]);

MaxCols = min(3,NumResponses);
NumRows = ceil(NumResponses/MaxCols);

for i = 1:NumResponses
    subplot(NumRows,MaxCols,i); 
     hold on; title([ DataStruct.ObjNames{i}],...
         'FontSize',FontSize); axis square;
    
    for j = 1:NumRealizations
        h1 = plot(DataStruct.time, DataStruct.data(j,:,i)',...
            'color',[0.5 0.5 0.5]);
    end
    xlabel('Time (days)','FontSize',FontSize);
    ylabel(DataStruct.name,'FontSize',FontSize);
    
    if (TruthRealization>0)
        h2 = plot(DataStruct.time, ...
            DataStruct.data(TruthRealization,:,i)','r','LineWidth',3);
        hlegend = legend([h1,h2],'Prior','Reference');
        axis tight;
        set(gcf,'color','w');
        set(gca,'FontSize',FontSize);
        set(hlegend,'Location','southwest');
    end
end

end

