function yearlyPlot(tvec,data)
% Function to plot yearly data and a linear regression
% Cas Renette, Léo Martin, Oslo University 2021.

[tvec_y, data_y] = timeAverages(tvec,data,'yearly', 'mean');

plot(tvec_y,data_y,'LineWidth',2)
hold on
plot(tvec_y,data_y,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
LRout =linReg_full(tvec_y, data_y, 'SouthEast', 0);
plot(LRout.xplot,LRout.yplot,'k')
text(LRout.xpos,LRout.ypos,LRout.text_out,'FontWeight','Bold')
datetick
hold off
end