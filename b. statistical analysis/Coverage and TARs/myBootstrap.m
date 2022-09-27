function [bootstrapData] = myBootstrap(data,N)

[m,~] = size(data);
idx = randi(m,[N,1]);
bootstrapData = data(idx,:);
end