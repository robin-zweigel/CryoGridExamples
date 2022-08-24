% Atmospheric stability function for momentum
% derived from theta functions in CLM5 documentation directly

function res = psi_M_CLM5(zeta1, zeta2)
zeta_m = -1.574;

if zeta1 <= 0 % Unstable
    
    if zeta1 < zeta_m % very unstable
        res = -log(-zeta_m) + log(-zeta1) - 1.14*( (-zeta2)^(1/3) - (-zeta_m)^(1/3) ) ...
            + (2*log(1 + (1-16*zeta_m)^.25) + log(1 + (1-16*zeta_m)^.5) -2*atan((1-16*zeta_m)^.25)) ...
            - (2*log(1 + (1-16*zeta2)^.25) + log(1 + (1-16*zeta2)^.5) -2*atan((1-16*zeta2)^.25));
    else
        res = (2*log(1 + (1-16*zeta1)^.25) + log(1 + (1-16*zeta1)^.5) -2*atan((1-16*zeta1)^.25) ) ...
            - (2*log(1 + (1-16*zeta2)^.25) + log(1 + (1-16*zeta2)^.5) -2*atan((1-16*zeta2)^.25) );
    end
else % stable
    if zeta1 > 1 % very stable
        res = -4*log(zeta1) - 5 - zeta1 + 1 + 5*zeta2;
    else
        res = 5*zeta2 - 5*zeta1;
    end
end

end