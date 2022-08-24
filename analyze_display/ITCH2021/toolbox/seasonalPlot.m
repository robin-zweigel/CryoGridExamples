function [outAv,outMat,years_vec] = seasonalPlot(tvec_i,data_i,period,method,disp,medianBool,cumBool)
% Function to plot a time serie as a repetition of yearly cycle based on
% daily or monthy averages.
% Author : Léo Martin, l.c.p.martin@uu.nl, October 2020, Oslo

if nargin<7
    cumBool=0;
    if nargin<6
        medianBool=0;
        if nargin<5
            disp=0;
        end
    end
end

if strcmp(period,'monthly')
    [outAv,outMat,years_vec] = monthlyPlot(tvec_i,data_i, method,disp,medianBool,cumBool);
else
    [outAv,outMat,years_vec] = dailyPlot(tvec_i,data_i, method,disp,medianBool,cumBool);
end

end

function [monthlyAv,monthlyMat,years_vec] = monthlyPlot(tvec_i,data_i, method,disp,medianBool, cumBool)
% Function to plot the seasonla signal in a variable

% Check input format
[tvec_i, data_i]=inputCheck(tvec_i,data_i);

% Calculate average
[tvec_a, data_a] = timeAverages(tvec_i,data_i,'monthly', method);

% Chop to keep full years
tochop=find(month(tvec_a)==1);
if tochop(1)>1
    tvec_a(1:tochop(1)-1)=[];
    data_a(1:tochop(1)-1)=[];
end

tochop=find(month(tvec_a)==12);
if tochop(end)<length(tvec_a)
    tvec_a(tochop(end)+1:end)=[];
    data_a(tochop(end)+1:end)=[];
end

% Reshape
nbyears=length(tvec_a)/12;
assert(nbyears==round(nbyears),'monthlyPlot : error in data partitioning')

monthlyMat=(reshape(data_a,[12, nbyears]))';
years_vec=unique(year(tvec_a));

if medianBool ==1
    monthlyAv=nanmedian(monthlyMat);
    fprintf('seasonalPlot: unsing median\n')
else
    monthlyAv=nanmean(monthlyMat);
end

if disp > 0
    cm = colormap(jet(size(monthlyMat,1)));
    for i_y=1:size(monthlyMat,1)
        if cumBool==1
            plot(1:12,cumsum(monthlyMat(i_y,:)),'Color', cm(i_y,:))
        else
            plot(1:12,monthlyMat(i_y,:),'Color', cm(i_y,:))
        end
        hold on
    end
    if cumBool==1
    plot(1:12,cumsum(monthlyAv),'k','LineWidth',3)
    plot(1:12,cumsum(monthlyAv),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
    else
    plot(1:12,monthlyAv,'k','LineWidth',3)
    plot(1:12,monthlyAv,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
    end
    xlim([1 12])
    set( gca, 'XTickLabel', {'J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'} )
    if disp ==1
        years_cell=num2cell(years_vec);
        years_cell=cellfun(@num2str,years_cell,'UniformOutput',false);
        years_cell(end+1)={'mean'};
        legend(years_cell);
    end
    hold off
end

end

function [dailyAv,dailyMat,years_vec] = dailyPlot(tvec_i,data_i, method,disp,medianBool,cumBool)
% Function to plot the seasonla signal in a variable on a daily basis

% Check input format
[tvec_i, data_i]=inputCheck(tvec_i,data_i);

% Calculate average
[tvec_a, data_a] = timeAverages(tvec_i,data_i,'daily', method);

% Chop to keep full years
tochop=find(day(tvec_a)==1 & month(tvec_a)==1);
if tochop(1)>1
    tvec_a(1:tochop(1)-1)=[];
    data_a(1:tochop(1)-1)=[];
end

tochop=find(day(tvec_a)==31 & month(tvec_a)==12);
if tochop(end)<length(tvec_a)
    tvec_a(tochop(end)+1:end)=[];
    data_a(tochop(end)+1:end)=[];
end

% remove leap day
isLeapday = month(tvec_a)==2 & day(tvec_a)==29;
tvec_a(isLeapday)=[];
data_a(isLeapday)=[];

% Reshape
nbyears=length(tvec_a)/365;
assert(nbyears==round(nbyears),'dailyPlot : error in data partitioning')

dailyMat=(reshape(data_a,[365, nbyears]))';
years_vec=unique(year(tvec_a));

if medianBool==1
    dailyAv=nanmedian(dailyMat);
    fprintf('seasonalPlot: unsing median\n')
else
    dailyAv=nanmean(dailyMat);
end

if disp > 0
    
    cm = colormap(jet(size(dailyMat,1)));
    tvec=datenum(2010,1,1)-1+(1:365);
    for i_y=1:size(dailyMat,1)
        if cumBool==1
            plot(tvec',cumsum(dailyMat(i_y,:)),'Color', cm(i_y,:))
        else
            plot(tvec',dailyMat(i_y,:),'Color', cm(i_y,:))
        end
        hold on
    end
    if cumBool==1
        plot(tvec,cumsum(dailyAv),'k','LineWidth',3)
    else
        plot(tvec,dailyAv,'k','LineWidth',3)
    end
    % plot(tvec,dailyAv,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
    xlim([datenum(2010,1,1) datenum(2010,12,31)])
    % set( gca, 'XTickLabel', {'J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'} )
    if disp ==1
        years_cell=num2cell(years_vec);
        years_cell=cellfun(@num2str,years_cell,'UniformOutput',false);
        years_cell(end+1)={'mean'};
        legend(years_cell);
    end
    datetick
    hold off
end

end