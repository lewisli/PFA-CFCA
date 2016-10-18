clear all;
close all;
addpath('../util');
addpath('../cfca');

addpath('/media/Scratch2/Data/SpatialForecasting/ensembled_realizations/')

% Load time data
load date.mat

% Load production data
load production_summary.mat


map_size = [200,200];
num_realizations = 950;
num_time_steps = length(date_common);
data_col = 1;
load facies_ensemble.mat

%% Load facies maps and plot some realizations
h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
for i = 1:4
    subplot(2,2,i)
    real = randi(num_realizations);
    imagesc(reshape(facies_ensemble(:,real),map_size));
    title(['Real ' num2str(real)]);
    set(gcf,'color','w');
    set(gca,'FontSize',22);
    %axis square
end
export_fig -m2 Oct13/GuangTI.png

%%

Data_Guang=cell(num_realizations,1);

for i = 1:num_realizations
    Data_Guang{i}.OIL1 = [production_summary_samet{i}(:,1) date_common];
    Data_Guang{i}.OIL2 = [production_summary_samet{i}(:,2) date_common];
    Data_Guang{i}.Field = [production_summary_samet{i}(:,1)+production_summary_samet{i}(:,2) date_common];
    Data_Guang{i}.Name = ['Run' num2str(i)];
end

%%
ForecastColumn = 1;
HistoricalColumn = 1;
TimeColumn = 2;
NumTimeSteps = 51;
EndTime = [3650];
ForecastObjName={'OIL1'};
HistoricalObjName = {'OIL1','OIL2'};
ForecastSpline = [4 20];
HistoricalSpline = [4 20];
ForecastStart = 4;
HistoricalEnd = 51;
PropertyNames = {'Oil Rate (stb/day)','Time (days)'};
[HistoricalStruct, ForecastStruct] = GenerateDataStructsWithInterpolation(Data_Guang, ...
    PropertyNames, ForecastColumn, HistoricalColumn, TimeColumn, ...
    HistoricalEnd,ForecastStart, NumTimeSteps, ForecastSpline, ...
    HistoricalSpline, ForecastObjName, HistoricalObjName, EndTime);


%%
TruthRealization = 0;
FontSize=22;
SavePath ='Oct13/PriorForecasts.png'
h1  = PlotInputResponse( HistoricalStruct,TruthRealization,FontSize,SavePath);

%%

HistoricalStruct.spline=[3 20]; % 6th order B-Spline with 20 knots
histPCA = ComputeHarmonicScores(HistoricalStruct,4);

%% Compute compression of images
load timelapse_ensemble.mat;
%%
Year = 8;
StartStep = Year*prod(map_size)+1;
EndStep = (Year+1)*prod(map_size);

saturation_maps = zeros(prod(map_size),num_realizations);

for Real = 1:num_realizations
    saturation_maps(:,Real) = timelapse_ensemble_tosave(StartStep:EndStep,Real);
end

%%
h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
for i= 1:2
    real = randi(num_realizations);
    subplot(1,2,i);
    imagesc(reshape(saturation_maps(:,real),map_size),[0,1]);
    title(['Real ' num2str(real)],'FontSize',FontSize);
    set(gcf,'color','w');
    set(gca,'FontSize',FontSize);
    colorbar;
end
export_fig -m2 Oct13/Satmaps.png


%% Compute eigen-images 
norm_saturation_maps = saturation_maps-repmat(mean(saturation_maps,2),1,num_realizations);
S = norm_saturation_maps'*norm_saturation_maps;
[V,D] = eig(S);
Eigenimages = normc(norm_saturation_maps*V);

%%
h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
for i= 1:2
    subplot(1,2,i);
    imagesc(reshape(Eigenimages(:,end-i+1),map_size),[0,0.05]);
    title(['Eigenimage ' num2str(i)],'FontSize',FontSize);
    set(gcf,'color','w');
    set(gca,'FontSize',FontSize);
    colorbar;
end
export_fig -m2 Oct13/eigenimages.png


%%
mean_image = mean(saturation_maps,2);

%%


%% Verify dimension reduction on saturation maps (reconstruct with only 90

num_weights = 100;
weights = zeros(num_realizations,num_weights);

for i = 1:num_realizations
   weight = Eigenimages'*(saturation_maps(:,i)-mean_image);
   weights(i,:)  = weight(end-num_weights+1:end);
end

realization = randi(num_realizations);
realization = 12;
EvStartIndex = num_realizations-num_weights + 1;
EvEndIndex = num_realizations;
h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
reconstructed = Eigenimages(:,EvStartIndex:EvEndIndex)*weights(realization,:)'+mean_image;
subplot(121);
imagesc(reshape(reconstructed,map_size),[0,1]);
title(['Reconstructed With ' num2str(num_weights) ' Eigenimages'],'fontsize',FontSize);
subplot(122);
imagesc(reshape(saturation_maps(:,realization),200,200),[0,1]);
title(['Original'],'fontsize',FontSize);
set(gcf,'color','w');
export_fig -m2 Oct13/100Reconstruction.png
%%
plot(cumsum(flipud(diag(D)))/sum(diag(D)),'linewidth',3);
set(gcf,'color','w');
xlabel('Eigenimage Components','FontSize',FontSize);
ylabel('Proportion of variance','FontSize',FontSize);
set(gca,'fontsize',FontSize);
export_fig -m2 Oct13/SatMapEigenvariance.png
%%



%% Eigenvalues are ordered from smallest to largest, we only care about largest
weightsFlipped = fliplr(weights);


%% Fundamentally we are going to be trying to find a correlation between the 
% harmonic scores in histPCA and 
close all;
TruthRealization = randi(num_realizations);
%TruthRealization = 12;

AvailableRealizations = setdiff(1:num_realizations,TruthRealization);

EigenvaluesToKeep = 100;
NumHistoricalResponses = 2
nHarmHist = 5;
for i = 1:NumHistoricalResponses
    harmscrhist=histPCA{i}.harmscr;
    
    % Need to re-arrange harmonic scores into Df such that the first
    % eigenvalues are placed in first
    for j = 1:nHarmHist
        Df(:,(i-1)*nHarmHist + j) = harmscrhist(AvailableRealizations,j);
        dobs_f(:,(i-1)*nHarmHist + j) = harmscrhist(TruthRealization,j);
    end
end


Hf = weightsFlipped(AvailableRealizations,1:EigenvaluesToKeep);


% Maximize correlation between 
[A,B,r,Dc,Hc] = canoncorr(Df,Hf);

% Project dobs_f into canonical space
dobs_c=(dobs_f-mean(Df))*A;

% Plot what things look like in canonical space
%PlotLowDimModels(Dc,Hc,dobs_c,'c',FontSize);
%export_fig -m2 Oct13/5EigenvalueCorr.png

%
% Compute Posterior
% Perform a normal score transform
epsilon = 0;
Hc_gauss = NormalScoreTransform(Hc,0);
C_H = cov(Hc_gauss);
H_CG_Mean = mean(Hc_gauss)';

% Find best linear bit between Dc and Hc_gauss
G = Dc'/Hc_gauss';
DDiff= Dc'-G*Hc_gauss';
C_T = DDiff*DDiff'/length(Dc);

C_Dc = zeros(size(C_T));

% Perform Gaussian Regression
mu_posterior = H_CG_Mean + C_H*G'*pinv(G*C_H*G' + C_T+C_Dc)*(dobs_c'-G*H_CG_Mean);
C_posterior = C_H - C_H*G'*inv(G*C_H*G' + C_T+C_Dc)*G*C_H;


% Reconstruct images
CycleSize = 200;
PosteriorSamples = mvnrnd(mu_posterior',C_posterior,CycleSize)';

% The h_c are not Gaussian, hence we need to backtransform
PosteriorSamplesTransformed = BackTransform(PosteriorSamples,Hc);

% We need to undo the canonical transform...
HpostCoef = PosteriorSamplesTransformed'*pinv(B)+repmat(mean(Hf,1)',...
            1,CycleSize)'; % H_f_posterior
        
        
% Need to flip to match the eigenimage ordering
sampledScores = fliplr(HpostCoef);

% Verify dimension reduction on saturation maps (reconstruct with only 90%)
realization = randi(CycleSize);
EvStartIndex = num_realizations-EigenvaluesToKeep + 1;
EvEndIndex = num_realizations;

% take mean of reconstructed images
reconstructed_posterior = zeros(prod(map_size),CycleSize);

for i = 1:CycleSize
    reconstructed_posterior(:,i) = ...
        Eigenimages(:,EvStartIndex:EvEndIndex)*sampledScores(i,:)'+mean_image;
end
%plot(cumsum(flipud(diag(D)))/sum(diag(D)),'linewidth',3);
h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
subplot(131);
imagesc(reshape(mean(reconstructed_posterior,2),map_size),[0,1]);
axis tight;
title('Posterior Mean','FontSize',FontSize);
subplot(132);
imagesc(reshape(reconstructed_posterior(:,100),map_size),[0,1]);
title('Posterior Realization','FontSize',FontSize);
axis tight;
subplot(133);
imagesc(reshape(saturation_maps(:,TruthRealization),map_size),[0,1]);
title('Truth','FontSize',FontSize);
axis tight;

export_fig -m2 Oct13/DF1_3.png
