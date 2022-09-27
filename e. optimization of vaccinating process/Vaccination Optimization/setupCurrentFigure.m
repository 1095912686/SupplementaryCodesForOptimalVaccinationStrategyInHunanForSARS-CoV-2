function fig = setupCurrentFigure(fig)
% setup figure format

axs = fig.CurrentAxes;
set(axs, 'FontSize', 16, 'FontName', 'Times New Roman');
set(get(axs, 'xlabel'), 'FontSize', 16, 'FontName', 'Times New Roman');
set(get(axs, 'ylabel'), 'FontSize', 16, 'FontName', 'Times New Roman');
set(get(axs, 'title'), 'FontSize', 16, 'FontName', 'Times New Roman');

end