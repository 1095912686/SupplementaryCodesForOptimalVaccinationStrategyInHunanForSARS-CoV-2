clear all; close all; clc;

% read data
ImportAgePartition;  groupCount = numel(agePartition)-1;
ImportParametersVEFIAQR;   % import parameters omega, gamma, p, kappa,...
ImportFigureLegends; % import the following cell array for figure legend:
                     % groupLegend; ageLegend;
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
C_full = readmatrix('contactMatrix.xlsx');   % the contact matrix

caseContactData = readtable('contactDataProcessed.xlsx');
caseContactData = caseContactData(~isnan(caseContactData.caseAge) & ~isnan(caseContactData.contactAge), :); %  extract data with (caseAge, contactAge) available
vaildCaseContactAgePairsCount = size(caseContactData, 1)  
[uniqueCaseName, idx] = unique(caseContactData.caseName);
uniqueCaseNameCount = numel(uniqueCaseName) 


% key variables
relation = caseContactData.relation;
setting = caseContactData.setting;

uniqueSettingCN = unique(setting); % uniqueSetting in Chinese
uniqueSettingEN = {'Empty', 'Public Transports', 'Amusement', 'School', 'Household', 'Hospital',...
    'Working Area', 'Service Industry', 'Uncategorized', 'Community', 'Catering industry'};  % uniqueSetting in English

% contacts count in each settings
contactCountInSettings = zeros(numel(uniqueSettingCN),1);
for k = 1:numel(uniqueSettingEN)
    contactCountInSettings(k) = sum( strcmp(setting, uniqueSettingCN{k}) );
end
contactCountInSettings
contactCountInNonEmptySettings = sum(contactCountInSettings(2:end)) % 


% caseCount in age groups
caseCountInGroups = findCaseCountInGroups(caseContactData(idx,:), agePartition)

%% decomposition
figureCount = 1;
C = zeros(11, groupCount, groupCount);
A = zeros(11, groupCount, groupCount);
for k = 1:numel(uniqueSettingEN)
    
    % extract rows corresponding to uniqueSetting{i}
    indices = strcmp(setting, uniqueSettingCN{k});
    caseContactDatai = caseContactData(indices,:);
    %N2(i) = numel(unique(caseContactDatai.caseName));
    
    % for caseContactData corresponding to uniqueSetting{i}, estimate the contact matrix via MLE
    [Ci, contactDataMatrix] = estimateCi(caseContactDatai, agePartition, caseCountInGroups);

    % replace nan by 0
    Ci(isnan(Ci)) = 0;
    contactDataMatrix(isnan(contactDataMatrix)) = 0;
   

%     % heatmap of contactDataMatrix
%     plotTitle = ['(', num2str(figureCount), '-', '1) ', uniqueSettingEN{i}, ', Contact Data Matrix'];
%     matrixplot(contactDataMatrix, plotTitle);  % heatmap of contactDataMatrix
%     exportgraphics(gca, ['contactDataMatrixNormalizedByAllCases', num2str(i), '.jpg'],'Resolution',300);  % export figure
% 
%     
%     % heatmap of contactMatrix
%     plotTitle = ['(', num2str(figureCount), '-', '2) ', uniqueSettingEN{i}, ', Contact Matrix Estimated via MLE'];
%     matrixplot(Ci, plotTitle);  % heatmap of contactDataMatrix
%     exportgraphics(gca, ['contactMatrixByMLENormalizedByAllCases', num2str(i), '.jpg'],'Resolution',300);  % export figure

    figureCount = figureCount + 1;
    
    C(k,:,:) = Ci;
    A(k,:,:) = contactDataMatrix;
    
    
end

C = cat(1, zeros(1,8,8), C); % adding a slice for the case without intervention

% %% visualize the recovered C and A 
% recoveredC = squeeze(sum(C,1));
% recoveredA = squeeze(sum(A,1));
% 
% figure; heatmap(recoveredC); title('Recovered Contact Matrix');
% exportgraphics(gca, 'recoveredContactMatrix.jpg','Resolution',300);  % export figure
% 
% figure; heatmap(recoveredA); title('Recovered Contact Data Matrix');
% exportgraphics(gca, 'recoveredContactDataMatrix.jpg','Resolution',300);  % export figure


%% simulate with each setting closed

% time span
tInit = 0;
tFinal = 14;

% step size
stepSize = 0.01;

% probability q of infection via a single contact
R0 = 8;
q_NGM = solveqFromR0(R0, N, C_full, 'NGM')
q_DBM = solveqFromR0(R0, N, C_full, 'DBM')


% converage and weights for simulation
coverageNameList = {'oldCoverage', 'newCoverage', 'percentageCoverage'};
weightNameList = {'uniformDelta', 'hospitalDelta', 'fatalityDelta',...
    'uniformOmicron', 'hospitalOmicron', 'fatalityOmicron'};


coverageNameRecord = {};
weightNameRecord = {};
settingNameRecord = {};
vaccinationStatusRecord = [];
ageGroupRecord = [];
weightedAccumulativeCasesRecord = [];

folderName = 'resultFigures';
mkdir(folderName)
cd(folderName)

for i = 1%:numel(coverageNameList)
    for j = 1:numel(weightNameList)

        switch coverageNameList{i}
            case 'oldCoverage'
                nVC = readmatrix('nVC.xlsx');
            case 'newCoverage'
                nVC = readmatrix('nVC_update.xlsx');
            case 'percentageCoverage'
                nVC = readmatrix('newVC2.xlsx');
        end


        switch weightNameList{j}
            case 'uniformDelta'
                weights = ones(8,4);

            case 'hospitalDelta'
                weights = hospital.Age.Delta * hospital.Vaccination.Delta.';

            case 'fatalityDelta'
                weights = fatality.Age.Delta * fatality.Vaccination.Delta.';

            case 'uniformOmicron'
                weights = ones(8,4);
                params.VE = params.VE / 2;  % VE for Omicron

            case 'hospitalOmicron'
                weights = hospital.Age.Omicron * hospital.Vaccination.Omicron.';
                params.VE = params.VE / 2;  % VE for Omicron

            case 'fatalityOmicron'
                weights = fatality.Age.Omicron * fatality.Vaccination.Omicron.';
                params.VE = params.VE / 2;  % VE for Omicron
        end


        % initial point
        feasibleVC = nVC ./ sum(nVC,2);
        xInit = [feasibleVC .* N, zeros(n,24)];
        doseCount = sum(xInit(:,1:4), 1);
        xInit(4,16) = 1;        % first case in I_{4,4}   aged 30-39, fully vaccinated
        xInit = xInit(:);



        %%% simulation
        weightedAccumulativeCases = zeros(numel(uniqueSettingEN)+1, 32); % initialize
        settingNames = ['Non-Intervention', uniqueSettingEN];
        for k = 1:numel(uniqueSettingEN)+1

            params.Beta   = (C_full - squeeze(C(k,:,:))) * q_NGM ./ (N');

            fun = @(x,t)dxdt_VEFIAQR_Vector(x,t,params);
            [x,t] = odeSolveRK4(fun, tInit, xInit, tFinal, stepSize);
            [V, E, F, I, A, Q, R] = extractVariableFromMatrix_forVEFIAQR(x);

            dailyNewCases = p(:)'.*omega(:)'.*E + (1-p(:)').*omegap(:)'.*E;
            weightedAccumulativeCases(k,:) = trapz(t, dailyNewCases) .* weights(:)'; % relatively weight, not normlized

            
            % save record
            for s = 1:8
                for t = 1:4
                    coverageNameRecord = [coverageNameRecord; coverageNameList{i}];
                    weightNameRecord = [weightNameRecord; weightNameList{j}];
                    settingNameRecord = [settingNameRecord; settingNames{k}];
                    ageGroupRecord = [ageGroupRecord; s];
                    vaccinationStatusRecord = [vaccinationStatusRecord; t];
                    weightedAccumulativeCasesRecord = [weightedAccumulativeCasesRecord; weightedAccumulativeCases(k,8*(t-1)+s)];
                end
            end

        end



        %% visualize reduction of different weight of accumulative cases in groups
        isplot = 1;
        if isplot == 1
            close all;

            % heatmap
            fig1 = figure;
            format shorte;
            xNames = doseLegend;
            yNames = [{'None'}, uniqueSettingEN];
            tempID = [1, 3:9, 11:12];
            h_DEF = heatmap(xNames, yNames(tempID), extractByVaccineDoses(weightedAccumulativeCases(tempID, :)));
            h_DEF.Title = 'Weighted Accumulative Cases in 14 Days Under different PSHM';
            h_DEF.XLabel = 'Vaccination Status';
            h_DEF.YLabel = 'Setting Closure';
            h_DEF.CellLabelFormat = '%.1e ';
            exportgraphics(fig1, ['weightedAccumulativeCasesWithin14DaysUnderDifferentPSHM_', coverageNameList{i}, '_', weightNameList{j}, '.jpg'],'Resolution',300);  % export figure


            % barplot for reduction of accumulative cases of entire population
            fig2 = figure;
            idx = [3:9, 11:12];
            y0 = sum(weightedAccumulativeCases(1,:)); % accumulative cases without PSHMs
            yValues = sum(weightedAccumulativeCases(idx,:),2); % accumulative cases with different PHSMs
            bari = bar(categorical(yNames(idx)), (y0 - yValues) ./ y0);

            xtips1 = bari(1).XEndPoints;
            ytips1 = bari(1).YEndPoints;
            labels1 = string(arrayfun(@(num)sprintf('%0.1f', num), bari(1).YData * 100, 'UniformOutput', false)) + '%';
            text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'bottom')
            ylabel('Fraction of Reduction');
            ylim([0,1]);
            title('Reduction of Weighted Accumulative Cases by Different PSHM');
            exportgraphics(fig2, ['reductionOfWeightedAccumulativeCasesWithin14DaysUnderDifferentPSHM_', coverageNameList{i}, '_', weightNameList{j},'.jpg'],'Resolution',300);  % export figure


            % barplots for reduction of accumulative cases in each age groups
            accumulativeCasesByAgeGroups = extractByAgeGroup(weightedAccumulativeCases);
            fig3 = figure;
            fig3.WindowState = 'maximized';
            tiled3 = tiledlayout("flow");

            for k = 1:size(accumulativeCasesByAgeGroups,2)
                idx = [3:9, 11:12];
                y0 = sum(accumulativeCasesByAgeGroups(1,k)); % accumulative cases without PSHMs
                yValues = sum(accumulativeCasesByAgeGroups(idx,k),2); % accumulative cases with different PHSMs

                nexttile;
                bari = bar(categorical(yNames(idx)), (y0 - yValues) ./ y0);
                xtips1 = bari(1).XEndPoints;
                ytips1 = bari(1).YEndPoints;
                labels1 = string(arrayfun(@(num)sprintf('%0.1f', num), bari(1).YData * 100, 'UniformOutput', false)) + '%';
                text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center',...
                    'VerticalAlignment', 'bottom')
                ylabel('Fraction of Reduction');
                ylim([0,1]);
                title(['age ', ageLegend{k}]);

            end
            sgtitle('Reduction of Weighted Accumulative Cases with Different PSHM by Different Age Groups')
            exportgraphics(fig3, ['weightedAccumulativeCasesWithin14DaysUnderDifferentPSHMByAgeGroup', coverageNameList{i}, '_', weightNameList{j}, '.jpg'],'Resolution',300);  % export figure


            % barplots for reduction of accumulative cases in each vaccination groups
            accumulativeCasesByVaccinationStatus = extractByVaccineDoses(weightedAccumulativeCases);
            fig4 = figure;
            fig4.WindowState = 'maximized';
            tiled4 = tiledlayout("flow");

            for k = 1:size(accumulativeCasesByVaccinationStatus,2)
                idx = [3:9, 11:12];
                y0 = sum(accumulativeCasesByVaccinationStatus(1,k)); % accumulative cases without PSHMs
                yValues = sum(accumulativeCasesByVaccinationStatus(idx,k),2); % accumulative cases with different PHSMs

                nexttile;
                bari = bar(categorical(yNames(idx)), (y0 - yValues) ./ y0);
                xtips1 = bari(1).XEndPoints;
                ytips1 = bari(1).YEndPoints;
                labels1 = string(arrayfun(@(num)sprintf('%0.1f', num), bari(1).YData * 100, 'UniformOutput', false)) + '%';
                text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center',...
                    'VerticalAlignment', 'bottom')
                ylabel('Fraction of Reduction');
                ylim([0,1]);
                title(doseLegend{k});

            end
            sgtitle('Reduction of Weighted Accumulative Cases with Different PSHM by Different Vaccination Status Groups')
            exportgraphics(fig4, ['weightedAccumulativeCasesWithin14DaysUnderDifferentPSHMByVaccinationGroup', coverageNameList{i}, '_', weightNameList{j}, '.jpg'],'Resolution',300);  % export figure
        end


    end
end


%% Make Record Table
Record = table;
Record.coverageName = coverageNameRecord;
Record.weightName = weightNameRecord;
Record.settingName = settingNameRecord;
Record.ageGroup = ageGroupRecord;
Record.vaccinationStatus = vaccinationStatusRecord;
Record.weightedAccumulativeCases = weightedAccumulativeCasesRecord;

writetable(Record, 'Record.xlsx');