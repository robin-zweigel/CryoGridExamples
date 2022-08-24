% Atmospheric stability function for heat/vapor
% From CLM5 documentation directly

function res = psi_H_CLM5(zeta1, zeta2)
zeta_h = -0.465;

if zeta1 <= 0 % Unstable
    if zeta1 < zeta_h % very unstable
        res = -log(-zeta_h) + log(-zeta1) - 0.8*( (-zeta_h)^(-1/3) - (-zeta1)^(-1/3) ) ...
            + 2*log( 1/2 + (1-16*zeta_h)^.5 ) - 2*log( 1/2 + (1-16*zeta2)^.5 );
    else
        res = 2*log( 1/2 + (1-16*zeta1)^.5 ) - 2*log( 1/2 + (1-16*zeta2)^.5 );
    end
else % stable
    if zeta1 > 1 % very stable
        res = -4*log(zeta1) - 5 - zeta1 + 1 + 5*zeta2;
    else
        res = -5*zeta1 + 5*zeta2 ;
    end
end

end