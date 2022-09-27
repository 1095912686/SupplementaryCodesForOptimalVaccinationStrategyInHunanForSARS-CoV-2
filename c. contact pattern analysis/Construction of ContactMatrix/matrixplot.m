function [] = matrixplot(matrix, TITLE, labels)
matrix(isnan(matrix)) = 0;
figure;
format shorte;
xvalue = labels;
yvalue = labels;
h_DEF = heatmap(xvalue, yvalue, matrix);
h_DEF.Title = TITLE;
h_DEF.XLabel = 'age group i';
h_DEF.YLabel = 'age group j';
h_DEF.FontName = 'times new roman';
h_DEF.FontSize = 12;
h_DEF.CellLabelFormat = '%.2f ';
end