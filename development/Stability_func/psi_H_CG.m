 % atmospheric stability function heat/water
% Original CryoGrid function

function res = psi_H_CG( zeta1, zeta2)
if zeta1<=0
    res= 1.9.*atanh((1 - 11.6.*zeta1).^0.5) + log(zeta1) - (1.9.*atanh((1 - 11.6.*zeta2).^0.5) + log(zeta2));
else
    res=((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta1) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta1))/2  - (((-5 + 5^0.5).*log(-3 + 5^0.5- 2.*zeta2) - (5 + 5^0.5).*log(3 + 5^0.5 + 2.*zeta2))/2);
end
end