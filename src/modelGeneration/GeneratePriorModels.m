% Generate3DSLModels.m
% Script to generate multiple StudioSL runs starting from some base case
clear all; close all;
clc
%% Set Paths
% Path to baseline case
BaseCaseDatPath = '../../data/3DSLFiles/CompressibleBaseCase.dat';

% Allocate a cell array that we will use to store the baseline
s=cell(GetNumberOfLines(BaseCaseDatPath),1);

fid = fopen(BaseCaseDatPath);
lineCt = 1;
tline = fgetl(fid);

while ischar(tline)
    s{lineCt} = (tline);
    lineCt = lineCt + 1;
    tline = fgetl(fid);
end

% Case Name
CaseName = 'Prior';
NbSimu = 500;     % Number of simulations
NbParams = 12;    % Number of parameters we will vary

% Fixed seed for reproducibility
rng('shuffle');

%% Set Parameter Ranges
% ParameterRanges is a struct that contains
Normal = 0;
Uniform = 1;

ParameterRanges = struct();
ParameterRanges.('Swir') = [0.2 0.05 Normal]; % Irreducible water saturation
ParameterRanges.('Swor') = [0.2 0.05 Normal]; % Irreducible oil saturation
ParameterRanges.('krw_end') = [0.3 0.1 Normal]; % End point water rel perm
ParameterRanges.('kro_end') = [0.7 0.1 Normal]; % End point oil rel perm
ParameterRanges.('no') = [2.5 0.2 Normal]; % Oil exponent
ParameterRanges.('nw') = [2 0.2 Normal]; % Oil exponent
ParameterRanges.('FaultMulti1') = [0.2 0.8 Uniform]; % Fault 1 trans multiplier
ParameterRanges.('FaultMulti2') = [0.2 0.8 Uniform]; % Fault 2 trans multiplier
ParameterRanges.('FaultMulti3') = [0.2 0.8 Uniform]; % Fault 3 trans multiplier
ParameterRanges.('FaultMulti4') = [0.2 0.8 Uniform]; % Fault 4 trans multiplier
ParameterRanges.('Viscosity') = [4 0.2 Normal];      % Oil viscosity
ParameterRanges.('OWC') = [1061 1076 Uniform];       % Oil water contact

%% Generate actual parameters
ParameterNames = fieldnames(ParameterRanges);
ParameterMatrix = struct();

rng('shuffle');

for i = 1:numel(ParameterNames)
    ParameterName = ParameterNames{i};
    ParameterRange = ParameterRanges.(ParameterName);
    
    if (i < 0)
        if (ParameterRange(3) == Normal)
            Value = ParameterRange(1) + ParameterRange(2)* ParametersValues(:,i);
        elseif(ParameterRange(3) == Uniform)
            Value = ParameterRange(1) + ParametersValues(:,i)*...
                (ParameterRange(2) - ParameterRange(1));
        end
    else
        if (ParameterRange(3) == Normal)
            Value = ParameterRange(1) + ParameterRange(2)*randn(NbSimu,1);
        elseif(ParameterRange(3) == Uniform)
            Value = ParameterRange(1) + rand(NbSimu,1)*...
                (ParameterRange(2) - ParameterRange(1));
        end
    end
    
    ParameterMatrix.(ParameterName) = Value;
    
    % Figure ID
   % i
    figureid = i-floor((i-1)/4)*4
    if (mod(i-1,4) == 0)
         figure(floor(i/4)+1);
         
         figurepic=floor(i/4)+1;
    end
    
    
    display([num2str(i) ' to ' num2str(figureid) 'on' num2str(figurepic)]);
    subplot(2,2,figureid);
    
    if (ParameterRange(3) == Normal)
        xx = linspace(ParameterRange(1)-ParameterRange(2),...
            ParameterRange(1)+ParameterRange(2));
        yy = normpdf(xx,ParameterRange(1),ParameterRange(2));
    elseif (ParameterRange(3) == Uniform)
        xx = linspace(ParameterRange(1),ParameterRange(2));
        yy = ones(size(xx))/(ParameterRange(2)-ParameterRange(1));
    end
    
    FontSize = 14;
    plot(xx,yy,'LineWidth',3);
    set(gcf,'color','w');
    set(gca,'FontSize',FontSize);
    xlabel(ParameterName,'FontSize',FontSize);
    ylabel('PDF','FontSize',FontSize);
    axis tight;
end

%save('/home/lewisli/code-dev/directforecasting/data/Compressible/RejectionSample/ParameterMatrix.mat','ParameterMatrix');

%% Construct Rel Perm curves to verify consistency
clear SW_corey_m;
RelPermEntries = 25;
SW_corey_m = zeros(NbSimu,RelPermEntries);
krw_model = zeros(NbSimu,RelPermEntries);
kro_model = zeros(NbSimu,RelPermEntries);

for i=1:NbSimu
    SW_corey_m(i,:) = linspace(ParameterMatrix.('Swir')(i), ...
        1 - ParameterMatrix.('Swor')(i),RelPermEntries);
end

for i=1:NbSimu
    krw_model(i,:) = ParameterMatrix.('krw_end')(i) .* ...
        ((SW_corey_m(i,:)-ParameterMatrix.('Swir')(i))./...
        (1-ParameterMatrix.('Swir')(i)-...
        ParameterMatrix.('Swor')(i))).^ParameterMatrix.('nw')(i);
    kro_model(i,:) = ParameterMatrix.('kro_end')(i) .* ...
        ((1 - SW_corey_m(i,:) - ParameterMatrix.('Swor')(i))./...
        (1 - ParameterMatrix.('Swir')(i) - ParameterMatrix.('Swor')(i))).^...
        ParameterMatrix.('no')(i);
    
    plot(SW_corey_m(i,:), krw_model(i,:),'b-', SW_corey_m(i,:), ...
        kro_model(i,:), 'r-');
    xlabel('S_w');
    ylabel('Relative Permeability');
    xlim([0 1]); ylim([0 1]); grid on; hold on;
end

%% Writing out the 'dat' files
TrialName = 'RejectionSample';
BaselineRunDirectory = ['/media/Scratch2/Data/3DSLRuns/Compressible/' ...
    TrialName '/'];

for k=1:NbSimu
    
    FolderNameIteration = [BaselineRunDirectory 'Run', num2str(k)];
    
    %creating an new folder for this iteration
    %checking if there is already a folder with that name
    if exist(FolderNameIteration,'dir') ~= 7
        mkdir(FolderNameIteration);
    end
    
    % Create new file
    file_name = [FolderNameIteration '/Run', num2str(k), '.dat'];
    fileID = fopen(file_name,'w+');
    
    % Loading everything before the Faultmultiplier
    MULTFLTIndex = SearchCellArray('MULTFLT',s);
    for j=1:MULTFLTIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    
    % Fault Multiplier
    F1Mult=['fault_1' blanks(1) ...
        num2str(ParameterMatrix.('FaultMulti1')(k)) blanks(1) '/'];
    F2Mult=['fault_2' blanks(1) ...
        num2str(ParameterMatrix.('FaultMulti2')(k)) blanks(1) '/'];
    F3Mult=['fault_3' blanks(1) ...
        num2str(ParameterMatrix.('FaultMulti3')(k)) blanks(1) '/'];
    F4Mult=['fault_4' blanks(1) ...
        num2str(ParameterMatrix.('FaultMulti4')(k)) blanks(1) '/'];
    
    fprintf(fileID,'%s',F1Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F2Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F3Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s',F4Mult);
    fprintf(fileID,'\n');
    fprintf(fileID,'/ \n');
    
    
    
    PVMULTIndex = SearchCellArray('PVMULT',s);
    CVISCOSITIESIndex = SearchCellArray('CVISCOSITIES',s);
    
    % Writing everything before Viscosity
    for j=PVMULTIndex:CVISCOSITIESIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Printing the PVTs out - in this case only viscosity is changed
    Visc=[num2str(ParameterMatrix.('Viscosity')(k)) blanks(1) '0.1 0.4 /'];
    fprintf(fileID,'%s',Visc);
    
    KRWOIndex = SearchCellArray('KRWO',s);
    
    
    % Writing everything before Viscosity
    for j=CVISCOSITIESIndex+1:KRWOIndex
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Writing out the Relperms
    
    formatSpecRelPerm = '%4.4f %4.8f %4.8f %s\n';
    fprintf(fileID,'%s\n', '--    Sw        krw       kro      Pc');
    for j=1:RelPermEntries
        fprintf(fileID,formatSpecRelPerm,SW_corey_m(k,j),...
            krw_model(k,j),kro_model(k,j));
        fprintf(fileID,'\n');
    end
    fprintf(fileID, '/\n');
    fprintf(fileID,'%s\n', 'END RELPERMS');
    %check if the slash works in the simulation
    
    INITIALCOND = SearchCellArray('INITIALCOND',s);
    OWC = SearchCellArray('OWC',s);
    
    % Writing everything before OWC
    for j=INITIALCOND:OWC
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    % Writing out the OWC
    OWCValue=['-', num2str(ParameterMatrix.('OWC')(k))];
    fprintf(fileID,'%s',OWCValue);
    fprintf(fileID,'\n');
    fprintf(fileID, '/\n');
    
    EndInitialCondition=SearchCellArray('END INITIALCOND',s);
    
    % Writing out the rest
    for j=EndInitialCondition:GetNumberOfLines(BaseCaseDatPath)
        fprintf(fileID,'%c',s{j});
        fprintf(fileID,'\n');
    end
    
    fclose(fileID);
end

%% Generate PBS Scripts
Cluster=2010;
PreparePBS(NbSimu,'Compressible',TrialName,...
    '/media/Scratch2/Data/3DSLRuns/Compressible/',Cluster)

%% Copy Include Files
IncludeFolders = {'depo/','include/'};
IncludeFoldersDir = '/media/Scratch2/Data/3DSLRuns/';
OutputPath = ['/media/Scratch2/Data/3DSLRuns/Compressible/' TrialName '/'];
for i = 1:length(IncludeFolders)
    Source = [IncludeFoldersDir IncludeFolders{i}];
    Dest = [OutputPath IncludeFolders{i}];
    copyfile(Source,Dest);
end


