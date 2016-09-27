function [ h ] = PlotLowDimModels( D,H,DObs,Type,FontSize)
%PlotLowDimModels Plot a low dimensional projection of the reservoir models
%   Plot the models in some low dimension space according to data and
%   forecast
% Inputs:
%   D: Data
%   H: Forecast
%   DObs: Observed data
%   Type: c or f for canonical space or functional space
% Return:
%   h: handle to figure
%
% Author: Lewis Li (lewisli@stanford.edu)
% Date: May 2nd 2016

if (nargin < 4)
    FontSize=12;
end

ScatterSize=100;
ObservedLineThickness=5;

h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
subplot(131);
hold on;
scatter(D(:,1),H(:,1),ScatterSize,'filled')
plot([DObs(1),DObs(1)],[min(H(:,1)),max(H(:,1))],'r-',...
    'LineWidth',ObservedLineThickness);
text(DObs(1) + abs(DObs(1))*0.25,min(H(:,1)) + ...
    abs(min(H(:,1)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
xlabel(['d_',num2str(1),'^' Type],'FontSize',FontSize);
ylabel(['h_',num2str(1),'^' Type],'FontSize',FontSize);
coeff = corrcoef(D(:,1),H(:,1));
% title(['d^' Type '_1 vs h^' Type '_1. \rho = ' ...
%     num2str(coeff(2)) ],'FontSize',FontSize);
set(gca,'FontSize',FontSize);
axis square; axis tight;

subplot(132);
hold on;
scatter(D(:,2),H(:,2),ScatterSize,'filled')
plot([DObs(2),DObs(2)],[min(H(:,2)),max(H(:,2))],'r-',...
    'LineWidth',ObservedLineThickness);
text(DObs(2) + abs(DObs(2))*0.25,min(H(:,2)) + ...
    abs(min(H(:,2)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
xlabel(['d_',num2str(2),'^' Type],'FontSize',FontSize);
ylabel(['h_',num2str(2),'^' Type],'FontSize',FontSize);
coeff = corrcoef(D(:,2),H(:,2));
% title(['d^' Type '_2 vs h^' Type '_2. \rho = ' ...
%     num2str(coeff(2)) ],'FontSize',FontSize);
set(gca,'FontSize',FontSize);
axis square; axis tight;

subplot(133);
hold on;
scatter(D(:,1),H(:,2),ScatterSize,'filled')
plot([DObs(1),DObs(1)],[min(H(:,2)),max(H(:,2))],'r-',...
    'LineWidth',ObservedLineThickness);
text(DObs(1) + abs(DObs(1))*0.25,min(H(:,2)) + ...
    abs(min(H(:,2)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
xlabel(['d_',num2str(1),'^' Type],'FontSize',FontSize);
ylabel(['h_',num2str(2),'^' Type],'FontSize',FontSize);
coeff = corrcoef(D(:,1),H(:,2));
% title(['d^' Type '_1 vs h^' Type '_2. \rho = ' ...
%     num2str(coeff(2)) ],'FontSize',FontSize);
set(gca,'FontSize',FontSize);
axis square; axis tight;
set(gcf,'color','w');

end

