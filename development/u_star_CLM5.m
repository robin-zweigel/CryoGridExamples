function u_star = u_star_CLM5(Lstar, uz)
            k =0.4;
            z =  10;
            z0 = 0.001;
            d = .3;
            zetam = -1.574; % zeta-boundary between unstable/very unstable, see CLM5 documentation
    
            zeta = (z-d)/Lstar;
            if zeta < zetam % very unstable
                u_star = k.*uz./( log(-zetam.*Lstar./z0) - psi_M_CLM5(zetam) + 1.14.*((-zeta)^(1/3)-(-zetam)^(1/3)) + psi_M_CLM5(z0./Lstar) );
            elseif zeta >= -1.575 && zeta < 0 % unstable
                u_star = k.*uz./( log((z-d)./z0) - psi_M_CLM5(zeta) + psi_M_CLM5(z0./Lstar) );
            elseif zeta >= 0 && zeta <= 1 % stable
                u_star = k.*uz./( log((z-d)/z0) +5*zeta - 5*z0/Lstar );
            else % very stable
                u_star = k.*uz./( log(Lstar./z0) + 5 + 5.*log(zeta) + zeta - 1 - 5.*z0./Lstar );
            end
                
end
        