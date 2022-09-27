clear; close all; clc;

maxIter = 1e4; % max iteration for bootstrapping
N = 1e4;   % size of bootstrap sample set
agePartition = [0,10,20,30,40,50,60,70,200];
n = 8;
ImportFigureLegends;
ImportAgePartition;
caseContactData = readtable('contactDataProcessed.xlsx');
[n,~] = size(caseContactData);
[caseNames, uniqueIndex] = unique(caseContactData.caseNameInCaseData); caseNames(1) = [];
caseAge = caseContactData.caseAge(uniqueIndex(2:end));
contactNames = caseContactData.contactName;
contactAge = caseContactData.contactAge;
vaccines = caseContactData(:,27:end);

data = zeros(n,3);
forDiscrimination = zeros(n,3);
for i = 1:n
    
    % skip the missing contactAges
    if isnan(contactAge(i))
        continue;
    end

    % how many doses this contact has taken
    dosei = ~isempty(vaccines.time1{i}) + ~isempty(vaccines.time2{i}) + ~isempty(vaccines.time3{i});
    
    % save whether the contact is a case, for further discrimination analysis
    isCase = sum(strcmp(contactNames(i), caseNames)) >= 1;
    forDiscrimination(i,:) = [isCase, contactAge(i), dosei];

    % if the contact is a case, save to nCC
    if isCase
        data(i,1:3) = [dosei, 1, contactAge(i)]; % [contactDoses, isCase, age]
        
    else
        data(i,1:3) = [dosei, 0, contactAge(i)];
    end
end

temp1 = forDiscrimination(forDiscrimination(:,1) == 1,:);
temp0 = forDiscrimination(forDiscrimination(:,1) == 0,:);
idx = randi(size(temp1,1), [size(temp0,1) - size(temp1,1), 1]);
forDiscrimination_sampleSizeBalanced = [temp0; temp1; temp1(idx,:)];
writematrix(forDiscrimination_sampleSizeBalanced,'forDiscrimination_sampleSizeBalanced.csv');

%% bootstrapping
pBootstraped = zeros(maxIter, 4);
tarRecord = zeros(maxIter, 8, 4);
for k = 1:maxIter % repeating bootstraping
    [bootstrapData] = myBootstrap(data, N);
    bootstrapVC = zeros(8,4);
    bootstrapCC = zeros(8,4);
    for i = 1:N 

        group = whichGroup(bootstrapData(i,3), agePartition);

        for j = 0:3 % doses
            if bootstrapData(i,1) == j
                bootstrapVC(group, j+1) = bootstrapVC(group, j+1) + 1;
                if bootstrapData(i,2) == 1
                    bootstrapCC(group, j+1) = bootstrapCC(group, j+1) + 1;
                end
            end
        end

    end

    bootstrapVC;
    bootstrapCC;
    tarRecord(k, :, :) = bootstrapCC ./ bootstrapVC; 
    tarRecord(isnan(tarRecord)) = 0;
    pBootstraped(k,:) = sum(bootstrapCC, 1) ./ sum(bootstrapVC, 1);
end

writematrix(pBootstraped,'pBootstraped.xlsx');


%% 

% reform the tarRecord to a suitable format for trainning
tarForTrainning = zeros(numel(tarRecord), 3);
count = 1;
for k = 1:maxIter
    for i = 1:8
        for j = 1:4
            tarForTrainning(count, :) = [tarRecord(k,i,j), i, j];  % [value, groupi, dosej];
            count = count + 1;
        end
    end
end

writematrix(tarForTrainning,'tarForTrainning.csv')
writematrix(squeeze(sum(tarRecord,2)),'TARdistributionInDoseGroups.xlsx');
writematrix(squeeze(sum(tarRecord,3)),'TARdistributionInAgeGroups.xlsx');

%% visualizing
% boxplot: TARInDoseGroups
writetable(array2table(squeeze(sum(tarRecord,2)), 'VariableNames', doseLegend), 'forBoxplotTARInVaccinationGroups.xlsx');

% fig1 = figure;
% b1 = boxplot(squeeze(sum(tarRecord,2)),'Notch','on');
% set(gca, 'XTickLabel', doseLegend, 'FontSize', 16);
% title('Total Attack Rate In Vaccination Status');
% ylabel('TAR')
% 
% fig1 = setupCurrentFigure(fig1);
% exportgraphics(fig1, 'TARInDoseGroups.jpg','Resolution',300);  % export figure
% savefig(fig1, 'TARInDoseGroups.fig');
% 
% % boxplot: TARInAgeGroups
% writetable(array2table(squeeze(sum(tarRecord,3)), 'VariableNames', ageLegend), 'forBoxplotTARInAgeGroups.xlsx');
% fig2 = figure;
% b2 = boxplot(squeeze(sum(tarRecord,3)),'Notch','on');
% ax = gca;
% ax.TickLabelInterpreter = 'tex';
% set(gca, 'XTickLabel', ageLegend, 'FontSize', 16);
% title('Total Attack Rate In Age Groups');
% ylabel('TAR')
% 
% fig2 = setupCurrentFigure(fig2);
% exportgraphics(fig2, 'TARInAgeGroups.jpg','Resolution',300);  % export figure
% savefig(fig2, 'TARInAgeGroups.fig');

fig = figure;
tiled = tiledlayout(1,2);

ax1 = nexttile;
b1 = boxplot(squeeze(sum(tarRecord,3)),'Notch','on');
ax = gca;
ax.TickLabelInterpreter = 'tex';
set(gca, 'XTickLabel', ageLegend, 'FontSize', 16, 'FontName', 'Times New Roman');
title('Total Attack Rate In Age Groups', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('TAR', 'FontSize', 16, 'FontName', 'Times New Roman')


ax2 = nexttile;
b2 = boxplot(squeeze(sum(tarRecord,2)),'Notch','on');
set(gca, 'XTickLabel', doseLegend, 'FontSize', 16, 'FontName', 'Times New Roman');
title('Total Attack Rate In Vaccination Status', 'FontSize', 16, 'FontName', 'Times New Roman');
ylabel('TAR', 'FontSize', 16, 'FontName', 'Times New Roman')


%fig = setupCurrentFigure(fig);
exportgraphics(fig, 'TARInAgeAndVaccinationGroups.jpg','Resolution',300);  % export figure
savefig(fig, 'TARInAgeAndVaccinationGroups.fig');


% quantDistribution = [0 0.25 0.5 0.75 1];
% b3 = boxPlot3D(tarRecord);
% xticks(1:8);
% xticklabels(ageLegend);
% yticklabels(doseLegend);

figure; 
scatter3(tarForTrainning(:,2),tarForTrainning(:,3),tarForTrainning(:,1),'bo');
xticks(1:8);
xticklabels(ageLegend);
yticklabels(doseLegend);