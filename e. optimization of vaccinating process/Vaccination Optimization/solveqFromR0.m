function q = solveqFromR0(R0, N, C, method)
%%% Solve the probability q of infection by a single contact From given R0 contact matrix C and other parameters in ODE model
% INPUT: 
%               R0      :   m*1 vector, contains a list of basic reproduction number for recover q
%               N       :   n*1 vector of population size of each groups
%               C       :   n*n matrix, the contact matrix of current population (n is the number of groups)
%            method     :   string of method ('NGM' for next generation method; 'DBM' for definition-based method)
% OUTPUT:
%               q       :   m*1 vector of probability of infected via single contact, under the corresponding R0
%

% parameters in ODE model
ImportParametersVEFIAQR;

% initialize
q_DBM = zeros(numel(R0),1);
q_NGM = zeros(numel(R0),1);
n = 8; % R0 dose not consider vaccination status, thus there are 8 age groups in total

% the matrix (R / q), still denote as R 
R = zeros(n); 
for i = 1:n
    for j = 1:n
        R(i,j) = C(j,i) / (p(i)*omega(i) + (1-p(i))*omegap(i)) * (kappa(i)*p(i)*omega(i)/gammap(i) + (1-p(i))*omegap(i)/gamma(i) + kappap(i)*(1-p(i))*omegap(i)/omegapp(i));
    end
end

% for all R0, compute q by two methods
for i = 1:numel(R0)
    temp1 = R .* (N.') ./ sum(N);
    temp2 = sum(temp1,1);
    q_DBM(i) = R0(i) / sum(temp2);

    temp3 = 1 / (p(i)*omega(i) + (1-p(i))*omegap(i)) * (kappa(i)*p(i)*omega(i)/gammap(i) + (1-p(i))*omegap(i)/gamma(i) + kappap(i)*(1-p(i))*omegap(i)/omegapp(i));
    q_NGM(i) = R0(i) / (max(eig(C)) * temp3);
end

% return
if strcmp(method, 'NGM')
    q = q_NGM;
elseif strcmp(method, 'DBM')
    q = q_DBM;
end

end
% fprintf('[R0, q_DBM, q_NGM] = \n');
% disp([R0, q_DBM, q_NGM]);