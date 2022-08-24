% Leaf resistances
% Following Dingman (2002) and STewart (1988)
% R. B. Zweigel, December 2021

figure

%% Light
Sin = 0:5:1000;
f_Sin = 1.1046.*Sin./(Sin + 104.4);

subplot(2,2,1)
plot(Sin,f_Sin)
title('Light')
xlabel('Incoming solar radiation [W/m^2]')
ylabel('f(S_{in})')
ylim([0 1])

%% Vapor-pressure deficit
d_rho_v = 0:0.0001:0.025;
f_d_rho_v = 0.233.*ones(size(d_rho_v));
I = d_rho_v <= 0.01152;
f_d_rho_v(I) = 1-66.6.*d_rho_v(I);

subplot(2,2,2)
plot(d_rho_v,f_d_rho_v)
title('Vapor-pressure deficit')
xlabel('absolute humidity deficit [kg/m^3]')
ylabel('f(\Delta\rho_v)')
ylim([0 1])

%% Leaf temperature
T = -10:.5:50;
f_T = zeros(size(T));
I = T>=0 & T<=40;
f_T(I) = T(I).*(40-T(I)).^1.18 ./ 691;

subplot(2,2,3)
plot(T,f_T)
title('Leaf temperature')
xlabel('Leaf temperature [{\circ}C]')
ylabel('f(T_{leaf})')
ylim([0 1])

%% Leaf water content
psi = 0:-5:-1500;
psi_wp = -1500;
f_psi = 1 - 0.00119.*exp(.81.*8.4.*psi./psi_wp);

subplot(2,2,4)
plot(psi,f_psi)
set(gca,'Xdir','reverse')
title('Leaf water content')
xlabel('Soil water potential [Pa]')
ylabel('f(\Psi)')
ylim([0 1])
% xlim([0 -1500])
