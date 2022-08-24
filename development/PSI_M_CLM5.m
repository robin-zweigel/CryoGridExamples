
function result = PSI_M_CLM5(zeta1, zeta2)
k = .4;
if zeta2 < -0.465 % very unstable
    result = real(log(abs(zeta1)) + 21./10.*k^(2/3).*(zeta1).^(1/3) - ( log(abs(zeta2)) + 21./10.*k^(2/3).*(zeta2).^(1/3)) );
elseif zeta2 < 0 & zeta2 >= -0.465 % unstable
    result = log((1-16.*zeta1).^.5+1) + 2.*log((1-16.*zeta1).^.25+1) - 2.*atan((1-16.*zeta1).^.25) - ...
        ( log((1-16.*zeta2).^.5+1) + 2.*log((1-16.*zeta2).^.25+1) - 2.*atan((1-16.*zeta2).^.25) );
elseif zeta2 >= 0 & zeta2 <= 1 % stable
    result = -5.*zeta1 + 5.*zeta2;
else % very stable
    result = -4.*log(zeta1)-zeta1 + 4.*log(zeta2) + zeta2;
end
end