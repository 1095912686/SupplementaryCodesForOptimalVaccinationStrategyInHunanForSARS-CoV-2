function reshapedI = extractByDosesAndGroups(I)


[m,~] = size(I);

reshapedI = zeros(m,8,4);

for i = 1:m
    reshapedI(i,:,:) = reshape(I(i,:),[8,4]);
end

end