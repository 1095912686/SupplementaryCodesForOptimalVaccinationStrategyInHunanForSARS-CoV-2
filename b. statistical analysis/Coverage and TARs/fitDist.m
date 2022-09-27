function [] = fitDist(data, binWidth, XLABEL, YLABEL, TITLE)
pd1_Gamma = fitdist(data,'Gamma')
pd1_Weibull = fitdist(data,'Weibull')
pd1_Lognormal = fitdist(data,'Lognormal')
a = pd1_Gamma.a;
b = pd1_Gamma.b;
A = pd1_Weibull.A;
B = pd1_Weibull.B;
mu = pd1_Lognormal.mu;
sigma = pd1_Lognormal.sigma;

x = linspace(0, max(data), 1e4);
y1_Gamma = pdf(pd1_Gamma,x);
y1_Weibull = pdf(pd1_Weibull,x);
y1_Lognormal = pdf(pd1_Lognormal,x);

figure;
plot(x,y1_Gamma,'LineWidth',2); hold on;
plot(x,y1_Weibull,'LineWidth',2); hold on;
plot(x,y1_Lognormal,'LineWidth',2); hold on;
h1 = histogram(data,'Normalization','pdf','FaceColor',[0 0.4470 0.7410],'binWidth',binWidth);

LEGEND = {['Gamma,     scale = ', num2str(b),'  shape = ', num2str(a)],...
    ['Weibull,     scale = ', num2str(A), '  shape = ',num2str(B)],...
    ['Lognormal, log-Scale = ', num2str(mu), '  log-Location = ', num2str(sigma)]};

legend(LEGEND);
xlabel(XLABEL);
ylabel(YLABEL);
title(TITLE);

end
