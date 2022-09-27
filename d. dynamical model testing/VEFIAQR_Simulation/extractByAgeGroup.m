function I_by_ages = extractByAgeGroup(I)

[m,~] = size(I);

I_by_ages = zeros(m,8);
idx = [0,1,2,3]*8 + 1;
for i = 1:8
    I_by_ages(:,i) = sum(I(:,idx+i-1), 2);
end

end