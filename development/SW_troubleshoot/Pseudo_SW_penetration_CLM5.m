% Pseudocode of Penetrate_SW_CLM5 for troubleshooting
figure
hold on

nir_fraction = .3;

for sun_angle = 0:.05:90
    L = 1;
    S = 0;
    alpha_leaf_nir = .35;
    alpha_leaf_vis = .07;
    alpha_stem_nir = .39;
    alpha_stem_vis = .16;
    tau_leaf_nir = .10;
    tau_leaf_vis = .05;
    tau_stem_nir = .001;
    tau_stem_vis =.001;
    Khi_L = .01;
    alpha_g = .15;
    
    spectral_weights = [nir_fraction 1-nir_fraction ];
    zenith = min(89,90 - sun_angle); % sun angles close to 90 produce NaNs
    
    alpha_leaf = [alpha_leaf_nir alpha_leaf_vis]; % leaf reflectances
    alpha_stem = [alpha_stem_nir alpha_stem_vis]; % stem reflectances
    tau_leaf = [tau_leaf_nir tau_leaf_vis]; % leaf transmittances
    tau_stem = [tau_stem_nir tau_stem_vis]; % stem transmittances
    
    w_leaf = L/(L+S); % leaf weighting
    w_stem = S/(L+S); % stem weighting
    alpha = alpha_leaf.*w_leaf + alpha_stem.*w_stem; % canopy reflectance, Eq. 3.11
    tau = tau_leaf.*w_leaf + tau_stem.*w_stem; % canopy transmittance, Eq.
    phi1 = .5 - .633.*Khi_L - .33.*Khi_L.^2; % for -.4 <= Khi_L <=.6
    phi2 = .877 .* (1-2.*phi1);
    my = cosd(zenith);
    G = phi1 + phi2.*my; % Relative projecter area of canopy (leaf and stem), Eq. 3.3
    K = G./my; % optical depth
    my_bar = 1./phi2*( 1- phi1./phi2 * log( (phi1+phi2)./phi1) ); % average inverse diffuse optical depth per unit leaf and stem area, Eq. 3.4
    omega = alpha + tau; % scattering coefficient
    cos_theta = (1+Khi_L)./2; % theta = mean leaf inclination angle relative to the horizontal plane, Eq. 3.14
    omega_beta = .5.*(alpha + tau + (alpha - tau).*cos_theta.^2); % upscatter for diffuse radiation, Eq. 3.13
    as = omega./2 .* G./max(my.*phi2+G,1e-6) .* (1 - my.*phi1./max(my.*phi2+G,1e-6) .* log((my.*phi1+max(my.*phi2+G,1e-6))./(my.*phi1)) ); % single scatter albedo, Eq. 3.16
    omega_beta_0 = (1+my_bar.*K)./(my_bar.*K).*as; % upscatter for direct beam radiation, Eq. 3.15
    beta_0 = omega_beta_0./omega;
    
    b = 1 - omega + omega_beta; % Eq. 3.31
    c = omega_beta; % Eq. 3.32
    d = omega.*my_bar.*K.*beta_0; % Eq.3.33
    f = omega.*my_bar.*K.*(1-beta_0); % Eq. 3.34
    h = sqrt(b.^2-c.^2)./my_bar; % Eq. 3.35
    sigma = (my_bar.*K).^2 + c.^2 - b.^2; % Eq. 3.36
    u1 = b - c./alpha_g; % Eq. 3.37
    u2 = b - c.*alpha_g; % Eq. 3.38
    u3 = f + c.*alpha_g; % Eq. 3.39
    s1 = exp(-min(h.*(L+S),40)); % Eq. 3.40
    s2 = exp(-min(K*(L+S),40)); % Eq. 3.41
    p1 = b + my_bar.*h; % Eq. 3.42
    p2 = b - my_bar.*h; % Eq. 3.43
    p3 = b + my_bar*K; % Eq. 3.44
    p4 = b - my_bar*K; % Eq. 3.45
    d1 = p1.*(u1-my_bar*h)./s1 - p2.*(u1+my_bar*h).*s1; % Eq. 3.46
    d2 = (u2 + my_bar.*h)./s1 - (u2 - my_bar*h).*s1; % Eq. 3.47
    h1 = -d.*p4 - c.*f; % Eq. 3.48
    h2 = 1./d1.*( (d-h1./sigma.*p3).*(u1-my_bar.*h)./s1 - p2.*(d-c-h1./sigma.*(u1+my_bar.*K)).*s2 ); % Eq. 3.50
    h3 = -1./d1.*( (d-h1./sigma.*p3).*(u1+my_bar.*h).*s1 - p1.*(d-c-h1./sigma.*(u1+my_bar.*K)).*s2 ); % Eq. 3.50
    h4 = -f.*p3 - c.*d; % Eq. 3.51
    h5 = -1./d2.*( h4.*(u2+my_bar.*h)./(sigma.*s1) + (u3-h4./sigma.*(u2-my_bar*K)).*s2 ); % Eq. 3.52
    h6 = 1./d2.*( h4./sigma.*(u2-my_bar.*h).*s1 + (u3-h4./sigma.*(u2 - my_bar.*K)).*s2 ); % Eq. 3.53
    h7 = c.*(u1-my_bar.*h)./(d1.*s1); % Eq. 3.54
    h8 = -c.*(u1+my_bar*h).*s1./d1; % Eq. 3.55
    h9 = (u2 + my_bar.*h)./(d2.*s1); % Eq. 3.56
    h10 = -s1.*(u2 - my_bar*h)./d2; % Eq. 3.57
    
    % Downwelling shortwave fluxes
    I_down_from_dir = (h4./sigma.*exp(-K*(L+S)) + h5.*s1 + h6./s1);%*spectral_weights'; % Eq. 3.19
    I_down_from_dif = (h9.*s1 + h10./s1)*spectral_weights'; % Eq. 3.20
    I_transmitted = exp(-K*(L+S));
    
%     S_down =  Sin_dir*(I_down_from_dir + I_transmitted) + Sin_dif*I_down_from_dif;
    
%     [seb.NEXT, S_up] = penetrate_SW(seb.NEXT, S_down);
    
    % Upwelling shortwave fluxes
    I_out_from_dir = (h1./sigma + h2 + h3);%*spectral_weights'; % Eq. 3.17
%     if sun_angle > 25 & sun_angle < 45
%         I_out_from_dir(2) = 0.06;
%     end
%   	if sun_angle > 30 & sun_angle < 50
%         I_out_from_dir(1) = 0.13;
%     end
    I_out_from_dif = (h7 + h8)*spectral_weights'; % Eq. 3.18
    
    plot(my,h1(1)/sigma(1),'.r')
    plot(my,h2(1),'.b')
    plot(my,h3(1),'g.')
    plot(my,I_out_from_dir(1),'k.')
end

grid on