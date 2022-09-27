n = 8;  % number of subgroups 
params.VE = 1 - readmatrix('OR2.xlsx');
params.kappa = 0.35 * ones(n,4);
params.kappap = 0.63 * ones(n,4);
params.p      = [0.31, 0.62, 0.62, 0.62] .* ones(n,4);
params.omega  = 1/3 * ones(n,4);
params.omegap = 1/3 * ones(n,4);
params.omegapp = 1/2 * ones(n,4);
params.gamma  = 1/5 * ones(n,4);
params.gammap = 1/7 * ones(n,4);
params.gammapp = 0 * ones(n,4);
params.mu = 0 * ones(n,4);
params.f      = 0 * ones(n,4);
params.br     = 0 * ones(n,4);
params.dr     = 0 * ones(n,4);  

kappa = params.kappa;
kappap = params.kappap;
p = params.p;
omega = params.omega;
omegap = params.omegap;
omegapp = params.omegapp;
gamma = params.gamma;
gammap = params.gammap;
gammapp = params.gammapp;
mu = params.mu;
f = params.f;
br = params.br;
dr = params.dr;
p = params.p;

% https://www.cdc.gov/mmwr/volumes/71/wr/mm7104e2.htm
N = readmatrix('agePopulationVector.xlsx'); % the population vector (stratified by age)
fatality.Vaccination.Delta = [0.02344, (0.02344+0.01351)/2, 0.01351, 0.01351*(0.004/0.00684)]' * 1e-2;	
fatality.Vaccination.Omicron = [0.01337, (0.01337+0.002165)/2, 0.002165, 0.002165*(0.004/0.00684)]' * 1e-2;
fatality.Vaccination.Mixed = [0.0224, (0.0224+0.00684)/2, 0.00684, 0.004]' * 1e-2;
fatality.Age.Delta = [0.001, 0.001, 0.004, 0.04, 0.05, 0.26, 1.12, (294+502)/(5943+3165)*100]' ./ N;
fatality.Age.Omicron = [0.005,0.001,0.002,0.005,0.01,0.05,0.2, (257+725)/(31066+14165)*100]' ./ N;

hospital.Vaccination.Delta = [0.1034, (0.1034+0.02309)/2, 0.02309, 0.0166]' * 1e-2;
hospital.Vaccination.Omicron = [0.0278, (0.0278+0.01054)/2, 0.01054, 0.0043]' * 1e-2;
hospital.Age.Delta = [0.46,0.33,1.33,1.5,1.4,2.38,5.29, (795+804)/(5943+3165)*100]' ./ N;
hospital.Age.Omicron = [1.1,0.38,0.6,0.66,0.59,0.77,1.39, (1108+1574)/(31066+14165)*100]' ./ N;

% normalization
fatality.Age.Delta = fatality.Age.Delta / sum(fatality.Age.Delta);
fatality.Age.Omicron = fatality.Age.Omicron / sum(fatality.Age.Omicron);
hospital.Age.Delta = hospital.Age.Delta / sum(hospital.Age.Delta);
hospital.Age.Omicron = hospital.Age.Omicron / sum(hospital.Age.Omicron);


