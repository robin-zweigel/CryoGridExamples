%Test albedo code based on full CTSM code:
%https://github.com/ESCOMP/CTSM/blob/master/src/biogeophys/SurfaceAlbedoMod.F90
%subroutine TwoStream
%Kjetil Aas
%07.12.2021

close all;
clearvars;

%Input
S=0.0;L=1.0;
p=1; %Dummy, values are now hardcoded to NET
c=1; %Dummy
nlevcan=1;

albtestval=0.15;
albsod=[albtestval albtestval]; %Test values!
albsoi=[albtestval albtestval]; %Test values!
albgrd=[albtestval albtestval]; %Test values!
albgri=[albtestval albtestval]; %Test values!

lSFonly = false; 
snowveg_affects_radiation=false;
fwet=0.0; fcansno=0.0; 
omegas=[1e10 1e10]; betads=1e10; betais=1e10; %Only for fwet>0.0

%NET Temperate/NET Boreal/NDT Boreal. [vis nir]
chiL=0.01;	alphaLeaf=[0.07 0.35]; alphaStem=[0.16 0.39]; tauLeaf=[0.05 0.10]; tauStem=[0.001 0.001];

%Arrays (mostly for saving output)
cosz_arr=[0.001 0.001:0.0001:1]; 

sigma_arr=zeros(2,length(cosz_arr));c1_arr=zeros(2,length(cosz_arr)); b_arr=zeros(2,length(cosz_arr));

albd_arr=zeros(2,length(cosz_arr));
ftid_arr=zeros(2,length(cosz_arr));
ftdd_arr=zeros(2,length(cosz_arr));
fabd_arr=zeros(2,length(cosz_arr));
albi_arr=zeros(2,length(cosz_arr));
ftii_arr=zeros(2,length(cosz_arr));
fabi_arr=zeros(2,length(cosz_arr));

albdSF_arr=zeros(2,length(cosz_arr));albiSF_arr=zeros(2,length(cosz_arr)); 
h1_arr=zeros(2,length(cosz_arr)); h2_arr=zeros(2,length(cosz_arr)); h3_arr=zeros(2,length(cosz_arr));
twostext_arr=zeros(2,length(cosz_arr));avmu_arr=zeros(2,length(cosz_arr));omegal_arr=zeros(2,length(cosz_arr));


%Initialize with nans
albdSF=[nan nan]; albiSF=[nan nan];
rho=[nan nan]; tau=[nan nan];

%Calculate tau and rho
elai=L; esai=S; mpe=0.00000001;
wl(p) = elai(p) / max( elai(p)+esai(p), mpe );
ws(p) = esai(p) / max( elai(p)+esai(p), mpe );
for ib = [1 2]
  rho(p,ib) = max( alphaLeaf(p,ib)*wl(p) + alphaStem(p,ib)*ws(p), mpe );
  tau(p,ib) = max( tauLeaf(p,ib)*wl(p) + tauStem(p,ib)*ws(p), mpe );
end 

%Independent of wave band
chil = chiL; 

%----Start of CTSM code (Twostream in SurfaceAlbedoMod.F90)----

for i=1:length(cosz_arr)
    cosz=cosz_arr(i);     
    %cosz=0.6104;
    phi1 = 0.5 - 0.633*chil - 0.330*chil*chil;
    phi2 = 0.877 * (1.-2.*phi1);
    gdir(p) = phi1 + phi2*cosz;
    twostext(p) = gdir(p)/cosz;
    avmu(p) = ( 1 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2;

    temp0(p) = max(gdir(p) + phi2*cosz,1.e-6);
    temp1 = phi1*cosz;
    temp2(p) = ( 1. - temp1/temp0(p) * log((temp1+temp0(p))/temp1) );
 
   % Loop over all wavebands to calculate for the full canopy the scattered fluxes
   % reflected upward and transmitted downward by the canopy and the flux absorbed by the
   % canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given
   % an underlying surface of known albedo.
   %
   % Output:
   % ------------------
   % Direct beam fluxes
   % ------------------
   % albd       - Upward scattered flux above canopy (per unit direct beam flux)
   % ftid       - Downward scattered flux below canopy (per unit direct beam flux)
   % ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
   % fabd       - Flux absorbed by canopy (per unit direct beam flux)
   % fabd_sun   - Sunlit portion of fabd
   % fabd_sha   - Shaded portion of fabd
   % fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
   % fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
   % ------------------
   % Diffuse fluxes
   % ------------------
   % albi       - Upward scattered flux above canopy (per unit diffuse flux)
   % ftii       - Downward scattered flux below canopy (per unit diffuse flux)
   % fabi       - Flux absorbed by canopy (per unit diffuse flux)
   % fabi_sun   - Sunlit portion of fabi
   % fabi_sha   - Shaded portion of fabi
   % fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
   % fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer


    for ib=[1 2]

          % Calculate two-stream parameters omega, betad, and betai.
          % Omega, betad, betai are adjusted for snow. Values for omega*betad
          % and omega*betai are calculated and then divided by the new omega
          % because the product omega*betai, omega*betad is used in solution.
          % Also, the transmittances and reflectances (tau, rho) are linear
          % weights of leaf and stem values.

        omegal = rho(p,ib) + tau(p,ib);
        asu = 0.5*omegal*gdir(p)/temp0(p) *temp2(p);
        betadl = (1.+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu;
        betail = 0.5 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) ...
             * ((1.+chil(p))/2.).^2) / omegal;

          if ( lSFonly )
             % Keep omega, betad, and betai as they are (for Snow free case or
             % when there is no snow
             tmp0 = omegal;
             tmp1 = betadl;
             tmp2 = betail;
          else
             % Adjust omega, betad, and betai for intercepted snow
             if (snowveg_affects_radiation) 
                tmp0 =   (1.-fcansno(p))*omegal        + fcansno(p)*omegas(ib);
                tmp1 = ( (1.-fcansno(p))*omegal*betadl + fcansno(p)*omegas(ib)*betads ) / tmp0;
                tmp2 = ( (1.-fcansno(p))*omegal*betail + fcansno(p)*omegas(ib)*betais ) / tmp0;
             else
                tmp0 =   (1.-fwet(p))*omegal        + fwet(p)*omegas(ib);
                tmp1 = ( (1.-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0;
                tmp2 = ( (1.-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0;
             end 
          end   % end Snow free

          omega(p,ib) = tmp0;
          betad = tmp1;
          betai = tmp2;

          % Common terms

          b = 1. - omega(p,ib) + omega(p,ib)*betai;
          c1 = omega(p,ib)*betai;
          tmp0 = avmu(p)*twostext(p);
          d = tmp0 * omega(p,ib)*betad;
          f = tmp0 * omega(p,ib)*(1.-betad);
          tmp1 = b*b - c1*c1;
          h = sqrt(tmp1) / avmu(p);
          sigma = tmp0*tmp0 - tmp1;
          p1 = b + avmu(p)*h;
          p2 = b - avmu(p)*h;
          p3 = b + tmp0;
          p4 = b - tmp0;

          % Absorbed, reflected, transmitted fluxes per unit incoming radiation
          % for full canopy

          t1 = min(h*(elai(p)+esai(p)), 40.);
          s1 = exp(-t1);
          t1 = min(twostext(p)*(elai(p)+esai(p)), 40.);
          s2 = exp(-t1);

          % Direct beam
          if ( ~ lSFonly )
             u1 = b - c1/albgrd(c,ib);
             u2 = b - c1*albgrd(c,ib);
             u3 = f + c1*albgrd(c,ib);
          else
             % Snow Free (SF) only 
             % albsod instead of albgrd here:
             u1 = b - c1/albsod(c,ib);
             u2 = b - c1*albsod(c,ib);
             u3 = f + c1*albsod(c,ib);
          end
          tmp2 = u1 - avmu(p)*h;
          tmp3 = u1 + avmu(p)*h;
          d1 = p1*tmp2/s1 - p2*tmp3*s1;
          tmp4 = u2 + avmu(p)*h;
          tmp5 = u2 - avmu(p)*h;
          d2 = tmp4/s1 - tmp5*s1;
          h1 = -d*p4 - c1*f;
          tmp6 = d - h1*p3/sigma;
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2;
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1;
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1;
          h4 = -f*p3 - c1*d;
          tmp8 = h4/sigma;
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2;
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2;
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2;
          if ( ~lSFonly )
            albd(p,ib) = h1/sigma + h2 + h3;
            ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1;
            ftdd(p,ib) = s2;
            fabd(p,ib) = 1. - albd(p,ib) - (1.-albgrd(c,ib))*ftdd(p,ib) - (1.-albgri(c,ib))*ftid(p,ib);
          else
            albdSF(p,ib) = h1/sigma + h2 + h3;
          end
          

          a1 = h1 / sigma * (1. - s2*s2) / (2. * twostext(p)) ...
             + h2         * (1. - s2*s1) / (twostext(p) + h) ...
             + h3         * (1. - s2/s1) / (twostext(p) - h);

          a2 = h4 / sigma * (1. - s2*s2) / (2. * twostext(p)) ...
             + h5         * (1. - s2*s1) / (twostext(p) + h) ...
             + h6         * (1. - s2/s1) / (twostext(p) - h);
          if ( ~ lSFonly )
            fabd_sun(p,ib) = (1. - omega(p,ib)) * ( 1. - s2 + 1. / avmu(p) * (a1 + a2) );
            fabd_sha(p,ib) = fabd(p,ib) - fabd_sun(p,ib);
          end

          % Diffuse
          if ( ~ lSFonly )
            u1 = b - c1/albgri(c,ib);
            u2 = b - c1*albgri(c,ib);
          else
             % Snow Free (SF) only 
             % albsoi instead of albgri here:
            u1 = b - c1/albsoi(c,ib);
            u2 = b - c1*albsoi(c,ib);
          end
          tmp2 = u1 - avmu(p)*h;
          tmp3 = u1 + avmu(p)*h;
          d1 = p1*tmp2/s1 - p2*tmp3*s1;
          tmp4 = u2 + avmu(p)*h;
          tmp5 = u2 - avmu(p)*h;
          d2 = tmp4/s1 - tmp5*s1;
          h7 = (c1*tmp2) / (d1*s1);
          h8 = (-c1*tmp3*s1) / d1;
          h9 = tmp4 / (d2*s1);
          h10 = (-tmp5*s1) / d2;

  
          % Final Snow Free albedo
          if ( lSFonly )
            albiSF(p,ib) = h7 + h8;
          else
            % For non snow Free case, adjustments continue
            albi(p,ib) = h7 + h8;
            ftii(p,ib) = h9*s1 + h10/s1;
            fabi(p,ib) = 1. - albi(p,ib) - (1.-albgri(c,ib))*ftii(p,ib);

            a1 = h7 * (1. - s2*s1) / (twostext(p) + h) +  h8 * (1. - s2/s1) / (twostext(p) - h);
            a2 = h9 * (1. - s2*s1) / (twostext(p) + h) + h10 * (1. - s2/s1) / (twostext(p) - h);

            fabi_sun(p,ib) = (1. - omega(p,ib)) / avmu(p) * (a1 + a2);
            fabi_sha(p,ib) = fabi(p,ib) - fabi_sun(p,ib);
  
            % Repeat two-stream calculations for each canopy layer to calculate derivatives.
            % tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are
            % calculated at the center of the layer. Derivatives are needed only for the
            % visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
            % Derivatives are calculated first per unit lai+sai and then normalized for sunlit
            % or shaded fraction of canopy layer.
  
            % Sun/shade big leaf code uses only one layer, with canopy integrated values from above
            % and also canopy-integrated scaling coefficients
  
            if (ib == 1) 
               if (nlevcan == 1) 
  
                  % sunlit fraction of canopy
                  fsun_z(p,1) = (1. - s2) / t1;
  
                  % absorbed PAR (per unit sun/shade lai+sai)
                  laisum = elai(p)+esai(p);
                  fabd_sun_z(p,1) = fabd_sun(p,ib) / (fsun_z(p,1)*laisum);
                  fabi_sun_z(p,1) = fabi_sun(p,ib) / (fsun_z(p,1)*laisum);
                  fabd_sha_z(p,1) = fabd_sha(p,ib) / ((1. - fsun_z(p,1))*laisum);
                  fabi_sha_z(p,1) = fabi_sha(p,ib) / ((1. - fsun_z(p,1))*laisum);
  
                  % leaf to canopy scaling coefficients
                  extkn = 0.30;
                  extkb = twostext(p);
                  vcmaxcintsun(p) = (1. - exp(-(extkn+extkb)*elai(p))) / (extkn + extkb);
                  vcmaxcintsha(p) = (1. - exp(-extkn*elai(p))) / extkn - vcmaxcintsun(p);
                  if (elai(p)  >  0.) 
                    vcmaxcintsun(p) = vcmaxcintsun(p) / (fsun_z(p,1)*elai(p));
                    vcmaxcintsha(p) = vcmaxcintsha(p) / ((1. - fsun_z(p,1))*elai(p));
                  else
                    vcmaxcintsun(p) = 0.;
                    vcmaxcintsha(p) = 0.;
                  end
  
               elseif (nlevcan > 1)
                  for iv = 1, nrad(p)
  
                     % Cumulative lai+sai at center of layer
  
                     if (iv == 1) 
                        laisum = 0.5 * (tlai_z(p,iv)+tsai_z(p,iv));
                     else
                        laisum = laisum + 0.5 * ((tlai_z(p,iv-1)+tsai_z(p,iv-1))+(tlai_z(p,iv)+tsai_z(p,iv)));
                     end
  
                     % Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction
     
                     t1 = min(h*laisum, 40.);
                     s1 = exp(-t1);
                     t1 = min(twostext(p)*laisum, 40.);
                     s2 = exp(-t1);
                     fsun_z(p,iv) = s2;
  
                     % ===============
                     % Direct beam
                     % ===============
  
                     % Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai
  
                     u1 = b - c1/albgrd(c,ib);
                     u2 = b - c1*albgrd(c,ib);
                     u3 = f + c1*albgrd(c,ib);
  
                     % Derivatives for h2, h3, h5, h6 and a1, a2
  
                     v = d1;
                     dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1;
  
                     u = tmp6 * tmp2 / s1 - p2 * tmp7;
                     du = h * tmp6 * tmp2 / s1 + twostext(p) * p2 * tmp7;
                     dh2 = (v * du - u * dv) / (v * v);
  
                     u = -tmp6 * tmp3 * s1 + p1 * tmp7;
                     du = h * tmp6 * tmp3 * s1 - twostext(p) * p1 * tmp7;
                     dh3 = (v * du - u * dv) / (v * v);
  
                     v = d2;
                     dv = h * tmp4 / s1 + h * tmp5 * s1;
     
                     u = -h4/sigma * tmp4 / s1 - tmp9;
                     du = -h * h4/sigma * tmp4 / s1 + twostext(p) * tmp9;
                     dh5 = (v * du - u * dv) / (v * v);
  
                     u = h4/sigma * tmp5 * s1 + tmp9;
                     du = -h * h4/sigma * tmp5 * s1 - twostext(p) * tmp9;
                     dh6 = (v * du - u * dv) / (v * v);
  
                     da1 = h1/sigma * s2*s2 + h2 * s2*s1 + h3 * s2/s1 ...
                         + (1. - s2*s1) / (twostext(p) + h) * dh2 ...
                         + (1. - s2/s1) / (twostext(p) - h) * dh3;
                     da2 = h4/sigma * s2*s2 + h5 * s2*s1 + h6 * s2/s1 ...
                         + (1. - s2*s1) / (twostext(p) + h) * dh5 ...
                         + (1. - s2/s1) / (twostext(p) - h) * dh6;
  
                     % Flux derivatives
     
                     d_ftid = -twostext(p)*h4/sigma*s2 - h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1;
                     d_fabd = -(dh2+dh3) + (1.-albgrd(c,ib))*twostext(p)*s2 - (1.-albgri(c,ib))*d_ftid;
                     d_fabd_sun = (1. - omega(p,ib)) * (twostext(p)*s2 + 1. / avmu(p) * (da1 + da2));
                     d_fabd_sha = d_fabd - d_fabd_sun;
  
                     fabd_sun_z(p,iv) = max(d_fabd_sun, 0.);
                     fabd_sha_z(p,iv) = max(d_fabd_sha, 0.);
  
                     % Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                     % to normalize derivatives by sunlit or shaded fraction to get
                     % APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
  
                     fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / fsun_z(p,iv);
                     fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / (1. - fsun_z(p,iv));
  
                     % ===============
                     % Diffuse
                     % ===============
  
                     % Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai
  
                     u1 = b - c1/albgri(c,ib);
                     u2 = b - c1*albgri(c,ib);

                     a1 = h7 * (1. - s2*s1) / (twostext(p) + h) +  h8 * (1. - s2/s1) / (twostext(p) - h);
                     a2 = h9 * (1. - s2*s1) / (twostext(p) + h) + h10 * (1. - s2/s1) / (twostext(p) - h);
     
                     % Derivatives for h7, h8, h9, h10 and a1, a2
  
                     v = d1;
                     dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1;
     
                     u = c1 * tmp2 / s1;
                     du = h * c1 * tmp2 / s1;
                     dh7 = (v * du - u * dv) / (v * v);
  
                     u = -c1 * tmp3 * s1;
                     du = h * c1 * tmp3 * s1;
                     dh8 = (v * du - u * dv) / (v * v);
  
                     v = d2;
                     dv = h * tmp4 / s1 + h * tmp5 * s1;
  
                     u = tmp4 / s1;
                     du = h * tmp4 / s1;
                     dh9 = (v * du - u * dv) / (v * v);
  
                     u = -tmp5 * s1;
                     du = h * tmp5 * s1;
                     dh10 = (v * du - u * dv) / (v * v);
  
                     da1 = h7*s2*s1 +  h8*s2/s1 + (1.-s2*s1)/(twostext(p)+h)*dh7 + (1.-s2/s1)/(twostext(p)-h)*dh8;
                     da2 = h9*s2*s1 + h10*s2/s1 + (1.-s2*s1)/(twostext(p)+h)*dh9 + (1.-s2/s1)/(twostext(p)-h)*dh10;
  
                     % Flux derivatives
  
                     d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1;
                     d_fabi = -(dh7+dh8) - (1.-albgri(c,ib))*d_ftii;
                     d_fabi_sun = (1. - omega(p,ib)) / avmu(p) * (da1 + da2);
                     d_fabi_sha = d_fabi - d_fabi_sun;
  
                     fabi_sun_z(p,iv) = max(d_fabi_sun, 0.);
                     fabi_sha_z(p,iv) = max(d_fabi_sha, 0.);
  
                     % Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                     % to normalize derivatives by sunlit or shaded fraction to get
                     % APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha
  
                     fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / fsun_z(p,iv);
                     fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / (1. - fsun_z(p,iv));
  
                    
                   end % end of iv loop
               end % nlevcan
            end   % first band
          end  % NOT lSFonly   
          
    sigma_arr(ib,i)=sigma; 
    
    albd_arr(ib,i)=albd(ib);
    ftid_arr(ib,i)=ftid(ib);
    ftdd_arr(ib,i)=ftdd(ib);
    fabd_arr(ib,i)=fabd(ib);
    albi_arr(ib,i)=albi(ib);
    ftii_arr(ib,i)=ftii(ib);
    fabi_arr(ib,i)=fabi(ib);
    
    albdSF_arr(ib,i)=albdSF(ib);
    albiSF_arr(ib,i)=albiSF(ib);
    h1_arr(ib,i)=h1;h2_arr(ib,i)=h2;h3_arr(ib,i)=h3;
    twostext_arr(ib,i)=twostext;avmu_arr(ib,i)=avmu;omegal_arr(ib,i)=omegal;
    if abs(sigma)<0.00000001
        disp('cosz, ib, sigma, h1, h4, h1*p3/sigma, h1/sigma, h4/sigma'); 
        disp([cosz,ib,sigma,h1,h4,h1*p3/sigma,h1/sigma,h4/sigma]); 
        return
    end
    
    end       % end of radiation band loop (ib)
end %cosz loop

%Plot results
ib=2; %Vis
figure;
plot(cosz_arr,h1_arr(ib,:)./sigma_arr(ib,:)); grid on; box on;hold on;
plot(cosz_arr,h2_arr(ib,:));plot(cosz_arr,h3_arr(ib,:));
plot(cosz_arr,h1_arr(ib,:)./sigma_arr(ib,:)+h2_arr(ib,:)+h3_arr(ib,:),'k')
legend('h1/sigma','h2','h3')

figure;hold on; box on; grid on;
plot(cosz_arr,albd_arr(ib,:));
plot(cosz_arr,ftid_arr(ib,:));
plot(cosz_arr,ftdd_arr(ib,:));
plot(cosz_arr,fabd_arr(ib,:));
plot(cosz_arr,albi_arr(ib,:));
plot(cosz_arr,ftii_arr(ib,:));
plot(cosz_arr,fabi_arr(ib,:));
legend('Upward scattered dir.','Downward scattered dir.','Transmitted direct',...
    'Flux absorbed dir.','Upward scattered dif.','Downward scattered dif.','Flux absorbed dif.')
title(['Lai=' num2str(L), ', Sai=', num2str(S), ', ib=', num2str(ib)])
xlabel('Zenit angle')

