function []=FormatPlot(Xlabel,Ylabel,Title)
% my format plot
% input xlabel and ylabel in braket

%legend modification:
% legend([D{1:3}],'bilayer, L=80','L=100','L=120');
% legend boxoff
% ah1=axes('position',get(gca,'position'),'visible','off');
% legend(ah1,[S{1:3}],'single layer, L=80','L=100','L=120');
% legend boxoff

set(get(gca,'XLabel'),'FontSize',12);   
set(get(gca,'YLabel'),'FontSize',12);
set(gca,'fontsize',16);
set(gca,'linewidth',2.5); 
set(get(gca,'Children'),'linewidth',2);
xlabel(Xlabel);
ylabel(Ylabel);
title(Title)
end