function I_by_doses = extractByVaccineDoses(I)


[m,~] = size(I);

I_by_doses = zeros(m,4);

for i = 1:4
    I_by_doses(:,i) = sum(I(:,(1:8)+(i-1)*8), 2);
end

end
