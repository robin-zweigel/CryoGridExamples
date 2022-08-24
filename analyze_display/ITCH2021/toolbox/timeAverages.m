function [tvec, averaged_data] = timeAverages(tvec_datetime,array,period, method)
% Function that computes time averages. Time vector should be in the
% datetime format. If not the function makes the conversion.
% L.C.P. Martin, Utrecht University, 2020.
%
% Use datetime_out=datetime(timevec,'ConvertFrom','datenum')
% to convert ISO proleptic time into datetime.
% Example of period: 'monthly', 'daily'
% Example of method: 'mean', 'sum'
%
% Function uses input with 1 row and many cols, so check this
% 
% Author : Léo Martin, l.c.p.martin@uu.nl, October 2020, Oslo

[array_nl,array_nc]=size(array);
[tvec_nl,tvec_nc]=size(tvec_datetime);
if array_nl>array_nc
    array=array';
end
if tvec_nl>tvec_nc
    tvec_datetime=tvec_datetime';
end

if isdatetime(tvec_datetime)==0
    tvec_datetime=datetime(tvec_datetime,'ConvertFrom','datenum');
end

% Do the calculation
[~,nbdate]=size(tvec_datetime);
[~,nc]=size(array);
assert(nc==nbdate, 'timeAverages: dimension problem');
averaged_data = array2timetable(array','RowTimes',tvec_datetime);
averaged_data = retime(averaged_data, period, method);
averaged_data = timetable2table(averaged_data);
tvec=datenum(averaged_data.Time);
averaged_data(:,1)=[];
averaged_data =table2array(averaged_data)';

% Fit input format
if array_nl>array_nc
    averaged_data=averaged_data';
end

end