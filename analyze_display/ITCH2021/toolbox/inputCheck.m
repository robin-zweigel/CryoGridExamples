function [tvec,signal] = inputCheck(tvec,signal)
% Function checking tvec and data have the same format and giving back
% column data.
% L.C.P. Martin, Utrecht University, 2021

[data_nl,data_nc]=size(signal);
[tvec_nl,tvec_nc]=size(tvec);
if data_nl<data_nc
    signal=signal';
end
if tvec_nl<tvec_nc
    tvec=tvec';
end
[nbdate,~]=size(tvec);
[nc,~]=size(signal);
assert(nc==nbdate, 'inputCheck: dimension problem');

end

