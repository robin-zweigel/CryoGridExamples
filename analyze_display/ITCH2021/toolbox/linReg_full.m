function LRout =linReg_full(x, y, position, plotbool)
% Function to do linear regression and plot relevant info.
% L.C.P. Martin, Utrecht University, 2020.

% Perform regression
LRout.p = fitlm(x,y);

% Produce points to plot
LRout.xplot=linspace(min(x),max(x),100);
LRout.yplot=LRout.p.Coefficients{1,1}+LRout.p.Coefficients{2,1}.*LRout.xplot;

% Text box
LRout.text_out=sprintf('r2 = %3.2f',LRout.p.Rsquared.Ordinary);
LRout.text_out=sprintf('%s \t delta  = %2.1e', LRout.text_out, LRout.yplot(end)-LRout.yplot(1));
LRout.text_out=sprintf('%s\na  = %3.2e', LRout.text_out, LRout.p.Coefficients{2,1});
LRout.text_out=sprintf('%s \t b  = %3.2e', LRout.text_out, LRout.p.Coefficients{1,1});

if nargin<3
    position='NorthWest';
end

if strcmp(position,'NorthEast')
    Xcoeff=0.50;
    Ycoeff=0.80;
elseif strcmp(position,'SouthEast')
    Xcoeff=0.80;
    Ycoeff=0.20;
elseif strcmp(position,'SouthWest')
    Xcoeff=0.20;
    Ycoeff=0.20;
else
    Xcoeff=0.20;
    Ycoeff=0.80;
end

LRout.xpos=min(x)+Xcoeff.*(max(x)-min(x));
LRout.ypos=min(y)+Ycoeff.*(max(y)-min(y));

if nargin<4
    plotbool=0;
end

if plotbool > 0
   % close all
   plot(x,y,'o')
   hold on
   plot(LRout.xplot,LRout.yplot)
   text(LRout.xpos,LRout.ypos,LRout.text_out,'FontWeight','Bold')
   if plotbool==2
       datetick
   end
   hold off
end

end