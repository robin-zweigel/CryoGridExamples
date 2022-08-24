function result = PSI_H_CLM5( zeta1, zeta2)
k = 0.4;
if zeta2 < -1.574 % very unstable
    result = log(abs(zeta1)) + real(27./10.*k^(4/3).*zeta1.^(1/3)) - ( log(abs(zeta2)) + real(27./10.*k^(4/3).*zeta2.^(1/3)) );
elseif zeta2 < 0 & zeta2 >= -1.574 % unstable
    result = 2.*log((1-16.*zeta1)^.5 + 1) - 2.*log((1-16.*zeta2)^.5 + 1);
elseif zeta2 >= 0 & zeta2 <= 1 % stable
    result = -5.*zeta1 + 5.*zeta2;
else % very stable
    result = -4.*log(zeta1)-zeta1 + 4.*log(zeta2) + zeta2;
end
end
