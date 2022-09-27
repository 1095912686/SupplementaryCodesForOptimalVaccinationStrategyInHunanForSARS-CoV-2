function [C_MLE, varargout] = estimateCi(data, agePartition, varargin)
% construct the contact matrix and contact data matrix for data with a given
% partition

    if nargin >= 3
        caseCountInGroups = varargin{1};
    end

    % caseAge, contactAge
    caseAge = data.caseAge;
    contactAge = data.contactAge;
    
    %%% Construct A, the Contact Data Matrix
    % contactDataMatrix(i,j) denotes the daily average number of close contact in group
    % j, will produced by a individual in group i.
    
    
    % age partition
    duration = 4 * ones(size(caseAge));  % case-wise duration of making contacts set as 4-day constant according to data
    if nargin < 3
        [contactDataMatrix, caseCountInGroups] = computeContactDataMatrix(data, duration, caseAge, contactAge, agePartition);
    else
        contactDataMatrix = computeContactDataMatrix(data, duration, caseAge, contactAge, agePartition, caseCountInGroups);
    end
    
    %%% read population Data
    populationData = readtable('populationData.xlsx');
    ageGroup = populationData.ageGroup;
    population = populationData.all;
    
    agePopulation = zeros(numel(agePartition)-1,1);
    for i = 1:numel(agePartition)-1
        agePopulation(i) = sum(population(ageGroup>=agePartition(i) & ageGroup < agePartition(i+1)),'omitnan');
    end
    
    % maximum likelihood estimation of contact matrix C
    C_MLE = estimateContactMatrix_MLE(contactDataMatrix,caseCountInGroups, agePopulation, 0);
    
    %exportgraphics(gca,'contactMatrixGroupedByAges_MLE.jpg','Resolution',300);
    
    % output contactDataMatrix as well
    if nargout >= 2
        varargout{1} = contactDataMatrix;
    end
    
    if nargout >= 3
        varargout{2} = caseCountInGroups;
    end

end
