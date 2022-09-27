clear all; close all; clc;

ImportParametersVEFIAQR;   % import parameters omega, gamma, p, kappa,...
ImportFigureLegends; % import the following cell array for figure legend:
                     % groupLegend; ageLegend;
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
C = readmatrix('contactMatrix.xlsx');   % the contact matrix

%% steup simulation

% initial point
%nVC = vaccineCoverageViaDataByIsComplete();
nVC = readmatrix('nVC.xlsx');
feasibleVC = nVC ./ sum(nVC,2);
xInit = [feasibleVC .* N, zeros(n,24)]; 
doseCount = sum(xInit(:,1:4), 1);
xInit(4,16) = 1;        % first case in I_{4,3}   aged 30-39, fully vaccinated
xInit = xInit(:);

% time span
tInit = 0;
tFinal = 100;

% step size
stepSize = 0.01;
    
% a series R0 for simulation
%R0 = [0.2, 0.5:6.5];     R0 = R0'; 
R0 = (1:8)';


%% probability vector q from R0(DBM and NGM are adopted respectively), for infection by a one-time contact
% 8 un-vaccinated age groups are used for construct the probability of infection via a single contact q from R0.
% all parameters except beta_{ij} are irrelevant to age groups, and we use

q_DBM = zeros(numel(R0),1);
q_NGM = zeros(numel(R0),1);

for i = 1:n
    for j = 1:n
        R(i,j) = C(j,i) / (p(i)*omega(i) + (1-p(i))*omegap(i)) * (kappa(i)*p(i)*omega(i)/gammap(i) + (1-p(i))*omegap(i)/gamma(i) + kappap(i)*(1-p(i))*omegap(i)/omegapp(i));
    end
end


for i = 1:numel(R0)
    temp1 = R .* (N.') ./ sum(N);
    temp2 = sum(temp1,1);
    q_DBM(i) = R0(i) / sum(temp2);
    
    temp3 = 1 / (p(i)*omega(i) + (1-p(i))*omegap(i)) * (kappa(i)*p(i)*omega(i)/gammap(i) + (1-p(i))*omegap(i)/gamma(i) + kappap(i)*(1-p(i))*omegap(i)/omegapp(i));
    q_NGM(i) = R0(i) / (max(eig(C)) * temp3);
    % R0 = max(eig(C)) * q / (p(i)*omega(i) + (1-p(i))*omegap(i)) * (kappa(i)*p(i)*omega(i)/gammap(i) + (1-p(i))*omegap(i)/gamma(i) + kappap(i)*(1-p(i))*omegap(i)/omegapp(i))
end
fprintf('[R0, q_DBM, q_NGM] = \n');
disp([R0, q_DBM, q_NGM]);



%% simulate for each self-defined base-line R0
for i = 1:numel(R0)
%     params.Beta   = C * q_DBM(i) ./ (N');
    params.Beta   = C * q_NGM(i) ./ (N');
    
    fun = @(x,t)dxdt_VEFIAQR_Vector(x,t,params);
    [x,t] = odeSolveRK4(fun, tInit, xInit, tFinal, stepSize);
    [V, E, F, I, A, Q, R] = extractVariableFromMatrix_forVEFIAQR(x);
    


    % dailyIncidence
    dailyNewSymptomaticCases_by_doses = extractByVaccineDoses(omegapp(:)'.*F) ./ doseCount;
    fig1 = figure(1);
    fig1.WindowState = 'maximized'
    subplot(2,ceil(numel(R0) / 2), i);
    plot(t, dailyNewSymptomaticCases_by_doses);
    legend(doseLegend);
    title(['R_0 = ', num2str(R0(i))]);
    sgtitle('Daily Incidence Rate by Vaccination Status with Simulated R_0', 'fontName', 'times new roman');
    xlabel('time (in days)');
    ylabel('Daily Incidence Rate (I) by Doses')
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    


    % plot(t,accumulativeNumberOfCases ./ doseCount);
    dailyNewCases = p(:)'.*omega(:)'.*E + (1-p(:)').*omegap(:)'.*E;
    accumulativeNumberOfCases = cumtrapz(t, dailyNewCases) * ones(32,1);
    TotalIncidenceRate = reshape(trapz(t, dailyNewCases,1), [8,4]) ./ (feasibleVC .* N);
    fig2 = figure(2);
    fig2.WindowState = 'maximized'
    subplot(2,ceil(numel(R0) / 2), i);
    TotalIncidenceRate(isnan(TotalIncidenceRate)) = 0;
    h_DEF = heatmap(doseLegend, ageLegend, TotalIncidenceRate);
    sgtitle('Total Attack Rate by Age and Vaccination Group with Simulated R_0', 'fontName', 'times new roman');
    title(['R_0 = ', num2str(R0(i))]);
    h_DEF.CellLabelFormat = '%.1e ';
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;

    % how to define and compute the accumulativeNumberOfCases from model directly 
    % instead of a extra numerical integration, it is still an very interesting question.
    Nt = sum(V + E + F + I + A + Q + R) / sum(N);
end

% exportgraphics(figure(1), 'simulatedIncidences14days.jpg','Resolution',300);  % export figure
% exportgraphics(figure(2), 'simulatedTARs14days.jpg','Resolution',300);  % export figure
exportgraphics(figure(1), 'simulatedIncidences100days.jpg','Resolution',300);  % export figure
exportgraphics(figure(2), 'simulatedTARs100days.jpg','Resolution',300);  % export figure
