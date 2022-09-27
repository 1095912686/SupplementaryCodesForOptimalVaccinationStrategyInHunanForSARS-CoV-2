function [] = vaccinationOptimizing(R0, C, nVC, weights, stepSize, tInit, tFinal, params)

ImportFigureLegends; % import the following cell array for figure legend: groupLegend; ageLegend;
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
n = numel(N);

%%% Initial Values for VEFIAQR model
    % initial point
    feasibleVC = nVC ./ sum(nVC,2);
    xInit = [feasibleVC .* N, zeros(n,24)];
    doseCount = sum(xInit(:,1:4), 1);
    xInit(4,16) = 1;        % first case in I_{4,3}   aged 30-39, fully vaccinated
    xInit = xInit(:);
    
    % Parameters
    q_NGM = solveqFromR0(R0, N, C, 'NGM');
    params.Beta   = C * q_NGM ./ (N');


    %%% stepwise vaccinating process
    % OD equations of this model
    fun = @(x,t)dxdt_VEFIAQR_Vector(x,t,params);

    % Init ObjFun (accumulative number of cases under init vaccine coverage)
    fInit0 = weightedAccumulativeCases(fun, params, tInit, xInit, tFinal, stepSize, weights);

    M = 1e8; % number of doses available
    m = 5e4; % stepSize for Vaccination
    record = zeros(floor(M/m), 4);
    stateRecord = zeros(floor(M/m),n,4);
    fRecord = zeros(floor(M/m)+1,1);
    fRecord(1) = fInit0;

    % stepwise vaccinating
    for k = 1 : floor(M/m)
        df = zeros(n,3); % decreasement under different vaccinating strategies
        fInit = zeros(n,3); % save the object function under the previous vaccine convergence
        m_k = zeros(n,3); % feasible stepsize for dose vaccination
        for i = 1:n     % number of age groups: n = 8
            for j = 1:3 % 0, 1, 2 doses
                [df(i,j), fInit(i,j), m_k(i,j)] = decisionEffectiveness(i, j, m, fun, params, tInit, xInit, fInit0, tFinal, stepSize, weights);
            end
        end

        % prepare for next iteration
        [temp, idx] = min(df(:));
        j = ceil(idx / 8);
        i = idx - 8*(j-1);
        dVC = feasibleDirection(i,j);
        xInit(1:32) = xInit(1:32) + dVC(:) * m_k(idx);
        fInit0 = fInit(idx);

        % disp(xInit(:,1:4))


        % save record
        fRecord(k+1) = fRecord(k) + temp * m_k(idx);
        record(k,:) = [i, j, df(idx), fInit(idx)];

        stateRecord(k,:,:) = reshape(xInit(1:32), [8,4]);
    end



    %% visualize the accumulative cases with different R0
    % initilize figure
    fig1 = figure(1)
    fig1.WindowState = 'maximized'
    if isempty(fig1.Children)
        tile1 = tiledlayout(3,4,'Padding','tight'); % weighted accumulative cases by vaccination process
    else
        tile1 = fig1.Children;
    end
    nexttile(tile1);
    plot((0:floor(M/m)) * m, fRecord, 'LineWidth', 1); hold on;
    subtitle(['R_0 = ', num2str(R0)]);
    title(tile1, ['Effectiveness of Vaccination'], 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel(tile1, 'Vaccination Process (in doses)', 'FontName', 'Times New Roman', 'FontSize', 16);
    ylabel(tile1, 'Weighted Accumulative Cases Within 14-days', 'FontName', 'Times New Roman', 'FontSize', 16);
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;



    %% cumulative summation of vaccination status
    % rearrange {'un-vaccinated', 'un-fully vaccinated', 'fully vaccinated', 'booster vaccinated'} 
    % to  {'booster vaccinated'， 'fully vaccinated'， 'un-fully vaccinated',' un-vaccinated'} 
    stateRecord = flip(stateRecord, 3); 

    % cumulative summation
    stateRecord = cumsum(stateRecord, 3);

    %% visualizing the stepwise vaccinating process for each R0
    
    % setup figure for optimal vaccinating process
    fign = figure;
    fign.WindowState = 'maximized';
    tilen = tiledlayout(2,4,'Padding','tight');

    % fign.Units = 'inches';
    % fign.OuterPosition = [0.25 0.25 3 3];

    markers = {'o', '^', 's', 'none'};
    lineStyles = {'-', '-', '-', '-.'};
    colors = {'#7E2F8E', '#0072BD', '#D95319', '#000000'}; % green, blue, red, black
%     RGB = cbrewer2('seq', 'Greens', 10, 'linear')
%     RGB = cbrewer2('seq', 'Blues', 10, 'linear')
%     RGB = cbrewer2('seq', 'Reds', 10, 'linear')
    % the 7-th rows of RGB
    faceColor = [8.0669e-01   9.2427e-01   7.8065e-01;...
                 7.9749e-01   8.7269e-01   9.4423e-01;...
                 9.9186e-01   7.6555e-01   6.7346e-01;...
                 1 1 1];




    % Plot vaccination history
    for i = 1:n

        % initialize the x-axis
        y0 = zeros(size(stateRecord, 1), 1); % initial lower bound of coloring
        x = m*(1:floor(M/m))';

        nexttile(tilen);
        hold on;

        % fill areas
        for j = 1:4
            y = stateRecord(:,i,j);

            % fill color
            p2 = fill([x; flip(x)], [y; y0], 'red');
            p2.FaceColor = faceColor(j,:);
            p2.FaceAlpha = 0.6;
            p2.EdgeColor = 'none';
            y0 = flip(y);

            % adjust fonts
            ax = gca;
            ax.FontName = 'Times New Roman';
            ax.FontSize = 14;
        end

        % plot lines
        for j = 1:4
            y = stateRecord(:,i,j);
            p1 = plot(x, y, 'LineWidth', 1, 'Color', colors{j}, 'LineStyle', lineStyles{j}, 'Marker', markers(j), 'MarkerIndices', 1:100:numel(x));
        end

        %legend(doseLegend);
        title(['Age ', ageLegend{i}]);
        ylim([0, y(end)]);
    end
  
    xlabel(tilen, 'Process (in doses)', 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(tilen, ['Population Size'], 'FontName', 'Times New Roman', 'FontSize', 18 , 'FontWeight', 'bold');
    % title(tilen, ['Optimized Vaccinating Process in Age Groups, R_0 = ', num2str(R0)], 'FontName', 'Times New Roman', 'FontSize', 18);
  




end