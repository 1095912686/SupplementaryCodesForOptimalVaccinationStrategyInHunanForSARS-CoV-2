clear; close all; clc;

%% import data and parameters
format shortE
ImportParametersVEFIAQR;   % import parameters omega, gamma, p, kappa,...
ImportFigureLegends; % import the following cell array for figure legend: groupLegend; ageLegend;
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
C = readmatrix('contactMatrix.xlsx');   % the contact matrix

%% for different vaccine coverage and weights
coverageNameList = {'oldCoverage', 'nullCoverage'};
weightNameList = {'uniformDelta', 'hospitalDelta', 'fatalityDelta',...
    'uniformOmicron', 'hospitalOmicron', 'fatalityOmicron'};
for i = 1:numel(coverageNameList)
    for j = 1:numel(weightNameList)

        coverageName = coverageNameList{i}
        % coverageName = 'oldCoverage';
        % coverageName = 'newCoverage';
        % coverageName = 'percentageCoverage';

        weightName = weightNameList{j};
        % weightName = 'uniformDelta';
        % weightName = 'uniformOmicron';
        % weightName = 'hospitalDelta';
        % weightName = 'hospitalOmicron';
        % weightName = 'fatalityDelta';
        % weightName = 'fatalityOmicron';

        switch coverageName
            case 'oldCoverage'
                nVC = readmatrix('nVC.xlsx');
            case 'nullCoverage'
                nVC = zeros(8,4);
                nVC(:,1) = 1;
        end


        switch weightName
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





        %% setup optimization
        % time span
        tInit = 0;
        tFinal = 14;

        % step size
        stepSize = 1;

        % R0
        R0 = (1:12)';


        %% optimizing
        close all;
        for k = 1:numel(R0)
            vaccinationOptimizing(R0(k), C, nVC, weights, stepSize, tInit, tFinal, params);
        end

        %% export graphics
        cd('H:\xmuph\0 HunanCovid-19 paper\Codes\e. optimization of vaccinating process\Vaccination Optimization');
        folderName = ['../resultFigures/', coverageName, '_', weightName];
        mkdir(folderName)
        cd(folderName)
        exportgraphics(figure(1), ['EffectivenessOfVaccination', '_', coverageName, '_', weightName,'.jpg'],...
            'Resolution',300,...
            'ContentType', 'vector');  % export figure
        for k = 1:12
            exportgraphics(figure(k+1), ['OptimizedVaccinatingProcessInAgeGroups', '_', coverageName, '_', weightName, '_', 'Within14Days_R', num2str(R0(k)), '.jpg'],'Resolution',300);  % export figure
        end

    end % end weightNameList{j}

end % end coverageNameList{i}


% for legend
markers = {'o', '^', 's', 'none'};
lineStyles = {'-', '-', '-', '-.'};
colors = {'#7E2F8E', '#0072BD', '#D95319', '#000000'}; % green, blue, red, black
figure;
for j = 1:4
    plot(rand(4,1), 'LineWidth', 1, 'Color', colors{j}, 'LineStyle', lineStyles{j}, 'Marker', markers(j));
    hold on;
end
legend({'booster vaccinated', 'fully vaccinated', 'un-fully vaccinated',' population size'}, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(gca, 'legends.jpg', 'Resolution', 600);