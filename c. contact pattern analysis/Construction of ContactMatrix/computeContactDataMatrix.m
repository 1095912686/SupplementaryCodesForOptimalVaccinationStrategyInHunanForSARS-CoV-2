function varargout = computeContactDataMatrix(data, duration, caseAge, contactAge, agePartition, varargin)
%%% Construct A, the Contact Data Matrix
% contactDataMatrix(i,j) denotes the daily average number of close contact in group
% j, will produced by a individual in group i. 
%
%   INPUT:
%                  caseCountInGroups   : 8*1 cell array of caseCount in each groups
%                           duration   : m*1 vector, case-wise duration of making contacts
%               [caseAge, contactAge]  : m*2 matrix of n-th (caseAge, contactAge) pairs
%                       agePartition   : (n+1) * 1 vector, partitioned population by n age groups
%
%   OUTPUT: 
%                                   A  : the contact data matrix, whose ij-th entries denotes the daily average
%                                        contact in group j that a case in group i will produce
%

if nargin == 6
    caseCountInGroups = varargin{end};
else
    caseCountInGroups = findCaseCountInGroups(data, agePartition)
end

% number of groups 
groupCount = numel(agePartition) - 1; 

% initialize the contactMatrix A
A = zeros(groupCount);

%%% add (caseAge, contactAge) data pairs to the contactDataMatrix row by row 
groupOfCase = whichGroup(caseAge, agePartition);
groupOfContact = whichGroup(contactAge, agePartition);
contactRowVectorOfGroupi = zeros(1,groupCount); % initialize

for i = 1:groupCount    % inspect cases in each age group
    idx = groupOfCase == i; % index of cases in group i

    % case-wise weight
    w = 1 ./ duration(idx);

    % fill case-contact data pairs of idx in a row vector of each contact
    % groups
    for j = 1:groupCount
        contactRowVectorOfGroupi(j) = sum((groupOfContact(idx) == j) .* w);
    end

    % update the i-th row of contact data matrix A
    A(i,:) = A(i,:) + contactRowVectorOfGroupi / caseCountInGroups(i);
end


varargout{1} = A;
if nargout > 1
    varargout{2} = caseCountInGroups;
end


end