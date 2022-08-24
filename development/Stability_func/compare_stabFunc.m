% Check the different stability functions of CryoGrid and CLM5 against
% each other, 
% R. B. Zweigel, January 2022

close all
clear all

% Parameters
z = 12;
z0 = 10^-3;
L = -100:0.01:100;

for  i = 1:length(L)
   CG_M(i) = psi_M_CG(z/L(i),z0/L(i)); 
   CG_H(i) = psi_H_CG(z/L(i),z0/L(i)); 
   CLM5_M(i) = psi_M_CLM5(z/L(i),z0/L(i)); 
   CLM5_H(i) = psi_H_CLM5(z/L(i),z0/L(i)); 
end

%%
figure
plot(L,CG_M,'.')
hold on
plot(L,CLM5_M,'.')
plot([z/-1.574 z/1],[0 0],'|k')
title('Stability function momentum')
xlabel('L_*')
ylabel('\psi_M')
legend({'CryoGrid','CLM5'})

figure
plot(L,CG_H,'.')
hold on
plot(L,CLM5_H,'.')
plot([z/-0.465 z/1],[0 0],'|k')
title('Stability function heat/vapor')
xlabel('L_*')
ylabel('\psi_H')
legend({'CryoGrid','CLM5'})




