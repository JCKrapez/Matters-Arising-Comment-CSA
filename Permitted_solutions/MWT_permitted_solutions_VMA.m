% MWT_Permitted_solutions_VMA.m

% Developed in MATLAB R2021a
% Author: J.-C. Krapez
% Version 1.0: 2026-02-25
% DOTA, ONERA, 13300 Salon de Provence, France
% jean-claude.krapez@onera.fr

% Calculation of the permitted solutions for the Multi-Wavelength Thermometry (MWT) underdetermined problem

clear all;
Size_label=12;
Size_ticks=10;
format longG
close all
for ifig=1:3
    figure(ifig);clf;
    h=gcf;set(h,'Position',[450*(ifig-1)+10 670 550 300]);
end


h=6.62607015e-34;  % Js
c=299792458 ; % m/s
C1=2*h*c^2;
C2=1.438776877e-2;

disp([' Calculation of the permitted solutions for the MWT (underdetermined) problem ']);
disp(['  ']);

disp([' To perform the analysis, several options are available and many parameters can be changed at will']);
disp([' The selection by the user is made by answering to a series of questions; suggested values are in parentheses, the default value is in square brackets']);
disp([' Just press ENTER to select the default value ']);

disp(['  ']);
disp([' Planck''s law or Wien approximation can be selected to calculate the "measured" radiance used as input for the analysis of the permitted solutions']);
disp([' Planck''s law was selected in Matters Arising (core paper)']);
avis_Planck_direct = input(' Select Wien (0) or Planck (1) for the calculation of the "measured" radiance  (0-1) [1] ?  ');
if isempty(avis_Planck_direct) avis_Planck_direct = 1;end
disp([' Planck''s law or Wien approximation can be selected to calculate the radiance used to infer the permitted emissivity spectrum from the "measured" radiance']);
disp([' Planck''s law was selected in Matters Arising (core paper)']);
avis_Planck_inv = input(' Select Wien (0) or Planck (1) for the calculation of the radiance for the inversion (0-1) [1] ?  ');
if isempty(avis_Planck_inv) avis_Planck_inv = 1;end


disp(['  ']);
disp([' Scenarios established from the paper by He et al.  Nat. Comm. 2025,16:10090  ']);
disp([' Spectral band: [0.45 µm; 0.65 µm]']);
disp([' Spectral resolution: 0.01 µm, i.e. 21 wavelengths from 0.45 to 0.65 µm']);
disp(['      - Scenario 1 : blackbody furnace  ']);
disp(['      - Scenario 2 : ceramic Al2O3 ']);
avis_case = input(' Select a scenario (1-2) [1] ?  ');
if isempty(avis_case) avis_case = 1;end

Temper = input('Enter the assumed TRUE temperature (K) (1850 was used for scenario 1, 1350 was used for scenario 2) [1850]  ');
if isempty(Temper) Temper = 1850;end

lambda_inf = 0.45;
lambda_sup = 0.65;
nlambda = 21;

v_lambda=linspace(lambda_inf,lambda_sup,nlambda)'*1e-6;
v_xticks=[0.45 0.5 0.55 0.6 0.65];
v_xticklabels={'0.45' '0.5' '0.55' '0.6' '0.65'};
avis_Marker_emiss = 1;
avis_Marker_rad = 1;
avis_Marker_sol = 3;
avis_log_rad = 0;
rad_exc=1;


nref=21;
v_iref=round([1:((nlambda-1)/(nref-1)):nlambda]);
iaff=v_iref;

if(avis_Marker_emiss==1)
    Marker_emiss='o';
elseif(avis_Marker_emiss==2)
    Marker_emiss='-';
elseif(avis_Marker_emiss==3)
    Marker_emiss='-o';
end
if(avis_Marker_rad==1)
    Marker_rad='o';
elseif(avis_Marker_rad==2)
    Marker_rad='-';
elseif(avis_Marker_rad==3)
    Marker_rad='-o';
end
if(avis_Marker_sol==1)
    Marker_sol='o';
elseif(avis_Marker_sol==2)
    Marker_sol='-';
elseif(avis_Marker_sol==3)
    Marker_sol='-o';
end

disp([' Selection of the emissivity allowable interval (according to Coates, 1981) ']);
emiss_max = input(' Maximum allowable emissivity (1 was used for both scenarios) ? [1] :  ');
if isempty(emiss_max) emiss_max = 1;end
emiss_min = input(' Minimum allowable emissivity (0.7 was used for scenario 1, 0.05 for scenario 2) ? [0.7] :  ');
if isempty(emiss_min) emiss_min = 0.7;end

function_LCN = @(v_lambda,T) (C1./(v_lambda.^5.*(exp(C2./(v_lambda*T))-1*avis_Planck_direct)));
function_inv_LCN = @(v_lambda,T) (C1./(v_lambda.^5.*(exp(C2./(v_lambda*T))-1*avis_Planck_inv)));

if(avis_case==1)
    disp([' For scenario 1, the ''true'' emissivity was assumed to be a linear function from 0.95 to 0.97 over the spectral band ']);
    v_emiss=linspace(0.95,0.97,nlambda)';
    avis_log_emiss = 0;         % Emissivity in lin scale
    v_yticks=[0.7 0.8 0.9 1];
    v_yticklabels={'0.7' '0.8' '0.9' '1'};
elseif(avis_case==2)
    disp([' For scenario 2, the ''true'' emissivity was taken from ']);
    disp([' Lisitsyn, A. V. et al. Near-infrared optical properties of a porous alumina ']);
    disp([' ceramics produced by hydrothermal oxidation of aluminum. Infrared ']);
    disp([' Physics & Technology 77, 162{170 (2016).  ']);
  
    v_emiss= [8.35 7.95 7.82 7.8 7.82 7.85 7.98 8.13 8.24 8.4 8.48 8.47 8.45 8.3 8.13 7.84 7.65 7.47 7.43 7.38 7.34]'*0.2/10.65;
    avis_log_emiss = 1;         % Emissivity in log scale
    v_yticks=[0.05 0.1 0.2 0.3 0.5 0.7 1];
    v_yticklabels={'0.05' '0.1' '0.2' '0.3' '0.5' '0.7' '1'};
end

v_CN= function_LCN(v_lambda,Temper);
v_L=v_emiss.*v_CN;

figure(1);clf;
plot(v_lambda(v_iref)/1e-6,v_emiss(v_iref),Marker_emiss,'Color','k','MarkerSize',5);hold on; 
set(gca,'XColor','k','Ygrid','on','YColor','k','Xgrid','on','XColor','k',...
    'Xlabel',text('String','Wavelength (µm)','Fontsize',Size_label,'FontName','Times'),...
    'Ylabel',text('String','Emissivity','Fontsize',Size_label,'FontName','Times'),...
    'Xlim',[v_xticks(1) v_xticks(end)],...
    'Ylim',[emiss_min emiss_max],...
    'Fontsize',Size_ticks);
xticks(v_xticks);
xticklabels(v_xticklabels);
yticks(v_yticks);
yticklabels(v_yticklabels);
set(gca,'XTickLabel', get(gca,'XTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
set(gca,'YTickLabel', get(gca,'YTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
figure(1);


figure(2);
if(avis_log_rad==0)
    plot(v_lambda(v_iref)/1e-6,rad_exc*v_L(v_iref)/1e9,Marker_rad,'Color','k','MarkerSize',5);hold on;
elseif(avis_log_rad==1)
    semilogy(v_lambda(v_iref)/1e-6,rad_exc*v_L(v_iref)/1e9,Marker_rad,'Color','k','MarkerSize',5);hold on;
end
set(gca,'XColor','k','Ygrid','on','YColor','k','Xgrid','on','XColor','k',...
    'Xlabel',text('String','Wavelength (µm)','Fontsize',Size_label,'FontName','Times'),...
    'Xlim',[v_xticks(1) v_xticks(end)],...
    'Fontsize',Size_ticks);
if(avis_log_rad==1)
    toto=get(gca,'Ylim');
    set(gca,'Ylim', [10^floor(log10(toto(1))) 10^ceil(log10(toto(2)))]),
end
set(gca,'Ylabel',text('String','Radiance (kW/m^{2}/µm/sr)','Fontsize',Size_label,'FontName','Times'));   % (10^{9} W/m^{3}/sr)

xticks(v_xticks);
xticklabels(v_xticklabels);
set(gca,'XTickLabel', get(gca,'XTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
set(gca,'YTickLabel', get(gca,'YTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
figure(2);

T_min=max(C2./v_lambda(v_iref)./log(emiss_max*C1./v_lambda(v_iref).^5./v_L(v_iref)+1*avis_Planck_inv));   
T_max=min(C2./v_lambda(v_iref)./log(emiss_min*C1./v_lambda(v_iref).^5./v_L(v_iref)+1*avis_Planck_inv));  
disp([' The lower limit for the permitted temperature is  : ' num2str(T_min)  ' K']);
disp([' The upper limit for the permitted temperature is  : ' num2str(T_max)  ' K']);


disp([' The permitted temperatures will be plotted by steps starting from a plot-reference temperature (not necessarily equal to the TRUE temperature)']);
Temper_aff = input(' Enter the plot-reference temperature (K) [T_true] :  ');
if isempty(Temper_aff) Temper_aff = Temper;end
T_step = input(' Enter the temperature step (K) (5 for scenario 1, 10 for scenario 2) [5] :  ');
if isempty(T_step) T_step = 5;end

nstep=floor((T_max-Temper_aff)/T_step);
T_test_sup=[Temper_aff:T_step:(Temper_aff+nstep*T_step) T_max];
n_test_sup=numel(T_test_sup);
nstep=floor((Temper_aff-T_min)/T_step);
T_test_inf=[Temper_aff:-T_step:(Temper_aff-nstep*T_step) T_min];
n_test_inf=numel(T_test_inf);

figure(3);clf;


for i_test_inf=1:n_test_inf-1
    emiss_test=(v_L)./(function_inv_LCN(v_lambda,T_test_inf(i_test_inf)));
    plot(v_lambda(iaff)/1e-6,emiss_test(iaff),Marker_sol,'Color','k','MarkerSize',3);hold on;
end
i_test_inf=n_test_inf;   % Corresponds to Tmin
emiss_test=(v_L)./(function_inv_LCN(v_lambda,T_test_inf(i_test_inf)));
color_T=[0 0. 1.];
plot(v_lambda(iaff)/1e-6,emiss_test(iaff),Marker_sol,'Color',color_T,'MarkerSize',5,'MarkerEdgeColor',color_T,...
    'MarkerFaceColor',color_T,'Linewidth',2);hold on;

for i_test_sup=2:n_test_sup-1
    emiss_test=(v_L)./(function_inv_LCN(v_lambda,T_test_sup(i_test_sup)));
    figure(3);
    plot(v_lambda(iaff)/1e-6,emiss_test(iaff),Marker_sol,'Color','k','MarkerSize',3);hold on;  % ,'Color',colorc(color_sel(profile+3),:),'Linestyle','-'
end

i_test_sup=n_test_sup;    % Corresponds to Tmax
emiss_test=(v_L)./(function_inv_LCN(v_lambda,T_test_sup(i_test_sup)));
color_T=[1 0. 0.];
plot(v_lambda(iaff)/1e-6,emiss_test(iaff),Marker_sol,'Color',color_T,'MarkerSize',5,'MarkerEdgeColor',color_T,...
    'MarkerFaceColor',color_T,'Linewidth',2);hold on;

set(gca,'XColor','k','Ygrid','on','YColor','k','Xgrid','on','XColor','k',...
    'Xlabel',text('String','Wavelength (µm)','Fontsize',Size_label,'FontName','Times'),...
    'Ylabel',text('String','Emissivity','Fontsize',Size_label,'FontName','Times'),...
    'Xlim',[v_xticks(1) v_xticks(end)],...
    'Fontsize',Size_ticks);
if(avis_log_emiss==0)
    set(gca,'Ylim',[emiss_min 1]);
else
    set(gca,'Ylim',[emiss_min 1],'YScale','log');
end

xticks(v_xticks);
xticklabels(v_xticklabels);
yticks(v_yticks);
yticklabels(v_yticklabels);
set(gca,'XTickLabel', get(gca,'XTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
set(gca,'YTickLabel', get(gca,'YTickLabel'), 'FontName', 'Times', 'Fontsize', Size_ticks);
text(lambda_sup+(lambda_sup-lambda_inf)*0.015,1,'$\hat{T}(K)$','Interpreter','Latex','Fontsize',Size_label,'FontName','Times');
text(lambda_sup-(lambda_sup-lambda_inf)*0.2,emiss_max*1.02,['$\hat{T}_{L}=$' ' ' num2str(round(T_min,1))],'Interpreter','Latex','Fontsize',Size_ticks,'FontName','Times');
text(lambda_sup-(lambda_sup-lambda_inf)*0.2,emiss_min*1.02,['$\hat{T}_{U}=$' ' ' num2str(round(T_max,1))],'Interpreter','Latex','Fontsize',Size_ticks,'FontName','Times');
for i_test_inf=1:2:n_test_inf-1
    text(lambda_sup+(lambda_sup-lambda_inf)*0.0,(v_L(end))./(function_inv_LCN(v_lambda(end),T_test_inf(i_test_inf))),['<' num2str(T_test_inf(i_test_inf))],'Fontsize',Size_ticks,'FontName','Times');
end
for i_test_sup=(1+3):3:n_test_sup-1
    text(lambda_sup+(lambda_sup-lambda_inf)*0.0,(v_L(end))./(function_inv_LCN(v_lambda(end),T_test_sup(i_test_sup))),['<' num2str(T_test_sup(i_test_sup))],'Fontsize',Size_ticks,'FontName','Times');
end


