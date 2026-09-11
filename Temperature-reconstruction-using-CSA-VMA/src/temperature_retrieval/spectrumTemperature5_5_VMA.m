% spectrumTemperature5_5_VMA.m

% The initial version of this function is
% spectrumTemperature5_5.m
% It was obtained from https://github.com/KirinShi/Temperature-reconstruction-using-CSA
% He, Y., Chen, M. K., Huang, M. et al. Dispersive Meta-lens Thermometry for High-temperature Measurements.
% Nature Communications 16, 10090 (2025). https://doi.org/10.1038/s41467-025-65171-7

% The analysis of the paper, next the code and data made available by the
% authors, led to the submission to Nature Communications of a document "Matters Arising"

% The changes made to the code were made by:
% J.-C. Krapez
% 2026-09-07
% DOTA, ONERA, 13300 Salon de Provence, France
% jean-claude.krapez@onera.fr

% VMA is the flag used to highlight the changes (VMA for Version Matters
% Arising):
% - modifications to enable the display of intermediate and additional results
% - introduction of global variables to share with other
% functions and scripts
% - correction of a bug in Chameleon.m  -> Chameleon_VMA.m

function[T_fmincon,E]=spectrumTemperature5_5_VMA(V_b,V_lamda_i,T_b,lambdalist,exposuretime_i,exposuretime_b,nop,iter,ub,lb)
global ii jj T_unknown iTC % VMA
global stdTr  % added variable containing the best fit, i.e. lowest std of Tr (fmin0) which was surprisingly overwritten by the last instruction in the He version of Chameleon.
global avis_permitted avis_synth avis_disp_emiss avis_bb_ceram
global min_allowed_EMR max_allowed_EMR
global  T_perm_min T_perm_max
global t_threshold  % introduction of a threshold to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
global avis_thalweg i1_thalweg i2_thalweg
Lamda_i = lambdalist;
C2=1.4338*10^(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VMA Part added to plot the cost function surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(avis_thalweg ==1)  % VMA plot of the cost function surface
    npts=100;
    D=C2./(Lamda_i([i1_thalweg i2_thalweg]).*T_b)-log((V_lamda_i([i1_thalweg i2_thalweg])/exposuretime_i)./(V_b([i1_thalweg i2_thalweg])/exposuretime_b));
    [X,Y] = meshgrid(logspace(log10(min_allowed_EMR),log10(max_allowed_EMR),npts),logspace(log10(min_allowed_EMR),log10(max_allowed_EMR),npts));
    A=C2 ./Lamda_i(i1_thalweg) ./ (X(:)'+D(1));
    A=[A;C2 ./Lamda_i(i2_thalweg) ./ (Y(:)'+D(2))];
    J=reshape(std(A),[npts npts]);
    figure(70);clf;
    surf(X,Y,J,'EdgeColor','none','FaceColor','interp');
    xlabel('Emissivity at \lambda_1');
    ylabel('Emissivity at \lambda_2');
    xlim([min_allowed_EMR max_allowed_EMR]);
    ylim([min_allowed_EMR max_allowed_EMR]);
    colorbar;
    view(2);
    
    
    T_perm_min=max(C2./Lamda_i([i1_thalweg i2_thalweg])./(log(max_allowed_EMR)+D));
    T_perm_max=min(C2./Lamda_i([i1_thalweg i2_thalweg])./(log(min_allowed_EMR)+D));
    disp([' Minimum permitted temperature : ' num2str(T_perm_min) ' K']);
    disp([' Maximum permitted temperature : ' num2str(T_perm_max) ' K']);
%         disp([' Type ENTER to resume ' ]);
%         pause
    T_fmincon=0;  % not used
    E=0;  % not used
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  VMA Part added to plot the signals, their ratio and the emissivity deduced
%  therefrom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(avis_disp_emiss>0)
    % VMA : figures added to draw the signals and their ratio
    if(avis_bb_ceram==1)
        if(iTC==1)
            siTC='-';
        elseif(iTC==2)
            siTC='--';
        elseif(iTC==3)
            siTC='-.';
        elseif(iTC==4)
            siTC=':';
        end
    end
    
    
    figure(70);%clf;
    if(avis_bb_ceram==0)
        plot(lambdalist,V_b,'k','Linewidth',2.5);% VMA
        plot(lambdalist,V_lamda_i,'k');hold on;% VMA
    elseif(avis_bb_ceram==1)
        semilogy(lambdalist,V_b,'k','Linewidth',2.5);% VMA
        semilogy(lambdalist,V_lamda_i,'LineStyle',siTC,'Color','k');hold on;% VMA
        set(gca, 'YScale', 'log')
    end
    
    grid on; xlim([0.45 0.65]);
    %VMA: title({'Signals (BB in black)'},'interpreter','latex'); % VMA
    xlabel('Wavelength (µm)'); % VMA
    ylabel('V (a.u.)');  % VMA
    
    figure(71);%clf;
    if(avis_bb_ceram==0)
    plot(lambdalist,V_lamda_i./V_b,'k');hold on;% VMA
    elseif(avis_bb_ceram==1)
    plot(lambdalist,V_lamda_i./V_b,'LineStyle',siTC,'Color','k');hold on;% VMA
    end
    %VMA: title({'Signal ratio'},'interpreter','latex'); % VMA
    xlabel('Wavelength (µm)'); % VMA
    ylabel('V/Vb (-)');  % VMA
    % pause(0.3) %VMA
    
    % VMA: Part added to draw the expected emissivity from V/Vb and the true temperature T_unknown
    if(avis_disp_emiss>=1)
        figure(73);
        emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_unknown-1./T_b));
        emiss_perm_min=min(emiss_perm);
        emiss_perm_max=max(emiss_perm);
        if(avis_bb_ceram==0)
            plot(lambdalist,emiss_perm,'-k');hold on;% VMA
        elseif(avis_bb_ceram==1)
            plot(lambdalist,emiss_perm,'LineStyle',siTC,'Color','k');hold on;% VMA
        end
        if(avis_disp_emiss==2)
            emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./(T_unknown+20)-1./T_b));
            emiss_perm_min=min(emiss_perm_min,min(emiss_perm));
            emiss_perm_max=max(emiss_perm_max,max(emiss_perm));
            if(avis_bb_ceram==0)
                plot(lambdalist,emiss_perm,'r');hold on;% VMA
            elseif(avis_bb_ceram==1)
                plot(lambdalist,emiss_perm,'LineStyle',siTC,'Color','r');hold on;% VMA
            end
            emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./(T_unknown-20)-1./T_b));
            emiss_perm_min=min(emiss_perm_min,min(emiss_perm));
            emiss_perm_max=max(emiss_perm_max,max(emiss_perm));
            if(avis_bb_ceram==0)
                plot(lambdalist,emiss_perm,'b');hold on;% VMA
            elseif(avis_bb_ceram==1)
                plot(lambdalist,emiss_perm,'LineStyle',siTC,'Color','b');hold on;% VMA
            end
        end
        grid on; xlim([0.45 0.65]);
        %        ylim([max(0,emiss_perm_min-0.12*(emiss_perm_max-emiss_perm_min)) emiss_perm_max+0.15*(emiss_perm_max-emiss_perm_min)]);
        xlabel('Wavelength (µm)'); % VMA
        if(avis_bb_ceram==0)
            ylabel('Emissivity ratio (-)');  % VMA
            ylim([0.7 1.3]);
        elseif(avis_bb_ceram==1)
            ylabel('Emissivity (-)');  % VMA
            ylim([0 3]);
        end
    end
    %     disp([' Type ENTER to resume ' ]);
    %     pause
    T_fmincon=0;  % not used
    E=0;  % not used
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VMA Part added to test a hypothetic case with constant emissivity (GREY MATERIAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(avis_synth ==1)
    disp([' Test performed by introducing a realistic and constant value for the emissivity (GREY MATERIAL)']); % VMA
    emiss_grey = input('Enter the constant value for the emissivity (in Matters Arising: 0.2 for alumina, 0.5 for ZrO2, 0.9 for SiC) ? [0.2] : ');
    if isempty(emiss_grey) emiss_grey = 0.2;end
    T_grey = input('Enter the true temperature of this grey material (in Matters Arising: 1380K for alumina and Zr02, 1300K for SiC) (K)? [1380] : ');
    if isempty(T_grey) T_grey = 1380;end
    V_lamda_i=V_b./exp(-C2./(Lamda_i.*T_b)).*exp(-C2./(Lamda_i.*T_grey))*emiss_grey;
end  % VMA end of part added for avis_synth ==1


D=C2./(Lamda_i.*T_b)-log((V_lamda_i/exposuretime_i)./(V_b/exposuretime_b));
m =C2 ./ Lamda_i;
TT = @(x) C2 ./Lamda_i ./ (x + D);
fun = @(x) std(C2 ./Lamda_i ./ (x + D)); % object function
dim = length(lambdalist);

v_ub = log(ub) * ones(1, dim);
v_lb = log(lb) * ones(1, dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VMA Part added to draw a selection of permitted solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(avis_permitted ==1)
    figure(72);
    %clf;
           
    T_perm_min=max(C2./Lamda_i./(log(max_allowed_EMR)+D));
    emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_perm_min-1./T_b));
    plot(lambdalist,emiss_perm,'b','Linewidth',2);hold on;
    Mlambda=find(lambdalist==max(lambdalist));
    text(lambdalist(Mlambda), emiss_perm(Mlambda), sprintf(' %.1f K', T_perm_min),'FontName','Times','fontsize',9);  
    
    T_perm_max=min(C2./Lamda_i./(log(min_allowed_EMR)+D));
    emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_perm_max-1./T_b));
    plot(lambdalist,emiss_perm,'r','Linewidth',2);hold on;
    text(lambdalist(Mlambda), emiss_perm(Mlambda), sprintf(' %.1f K', T_perm_max),'FontName','Times','fontsize',9);  
    grid on; xlim([0.45 0.65]);
    if(avis_bb_ceram==0)
        ymin=max(0,min_allowed_EMR-0.12*(max_allowed_EMR-min_allowed_EMR));
        ymax=max_allowed_EMR+0.15*(max_allowed_EMR-min_allowed_EMR);
    else
        ymin=0;
        ymax=2.2;
    end
    ylim([ymin ymax]);
    xlabel('Wavelength (µm)'); % VMA
    
    if(avis_bb_ceram==0)
        ylabel('Emissivity ratio (-)');  % VMA
    elseif(avis_bb_ceram==1)
        ylabel('Emissivity (-)');  % VMA
    end
    %text(0.47, 0.5,...
    %sprintf('Pixel : ii= %.0f /41, jj= %.0f /37', ii,jj),'FontName','Times','fontsize',10)
    %     text(0.47, 0.5,...
    %     sprintf('Pixel : ii=%.0f,  jj=%.0f ', ii,jj),'FontName','Times','fontsize',10)
    
    % Plot of the curve corresponding to the true temperature (T_grey in the
    % case of synthetic grey data, T_unknown in the case of experimental
    % Al2O3)
    if(avis_synth==1)
        emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_grey-1./T_b));
    else
        emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_unknown-1./T_b));
    end
    plot(lambdalist,emiss_perm,'Color', [0, 0.8, 0] ,'Linewidth',2);hold on;  % green
    
    
    if(T_perm_max>T_perm_min+10)  % intermediate Permitted solutions will be drawn only in this case
        vsteps=[10 25 50 100];  % Possible values of the temperature steps
        vp=find(round(((T_perm_max-T_perm_min)/5)./vsteps)==1);   % /3  or /5
        step=vsteps(vp(1)); % VMA : if iq=1,2,3,4, the steps will be of 10, 25, 50, 100K
        pT_min=(round(T_perm_min/step)+1)*step;
        if(pT_min<=T_perm_min)
            pT_min=(round(T_perm_min/step)+1)*step;
        end
        
        pT_max=(round(T_perm_max/step)-1)*step;
        if(pT_max>=pT_min)
            for pT=pT_min:step:pT_max
                emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./pT-1./T_b));
                plot(lambdalist,emiss_perm,'k');hold on;
                    text(lambdalist(Mlambda), emiss_perm(Mlambda), sprintf(' %.1f K', pT),'FontName','Times','fontsize',9);   % VMA: for BB
            end
        end
    end
    if(avis_synth==1)
        text(0.52, ymin+(ymax-ymin)*0.11, sprintf('"True" temperature: %.1f K', T_grey),'FontName','Times','fontsize',10);
    else
        text(0.52, ymin+(ymax-ymin)*0.11, sprintf('"True" temperature: %.1f K', T_unknown),'FontName','Times','fontsize',10);
    end
    text(0.52, ymin+(ymax-ymin)*0.05, sprintf('Emissivity(-ratio) limits: [%.2f, %.1f] ', min_allowed_EMR,max_allowed_EMR),'FontName','Times','fontsize',10);
    
    disp(['Pixel : ' num2str(ii)  ', ' num2str(jj) ]); % VMA
    
    disp(['  T_perm such that max(emiss) = ub : ' num2str(T_perm_min)  'K, spectrum in blue' ]);
    disp(['  T_perm such that min(emiss) = lb : ' num2str(T_perm_max)  'K, spectrum in red'  ]);
    
    if(T_perm_min>T_perm_max)
        disp([' There is no possibility to have simultaneously max(emiss) < ub and min(emiss) > lb ']);
        disp([' There is no solution to the problem !!! ']);
        %pause(0.5);
    end
    
    if(avis_bb_ceram==1)
        figure(72);
        set(gca, 'YScale', 'log')
        figure(72);
    end
%     avis_linlog = input('Logarithmic scale for emissivity (0-1) ? [0] : ');
%     if isempty(avis_linlog) avis_linlog = 0;end
%     if(avis_linlog==1)
%         set(gca, 'YScale', 'log')
%         figure(72);
%     end
%     disp([' Press ENTER to resume ']);
%     pause
    
end  % VMA End of part added to draw the permitted solutions


% CSA parameters
noP = nop;
maxIter = iter;

[xmin, bestPosition, CSAConvCurve] =Chameleon_VMA(noP,maxIter,v_lb,v_ub,dim,fun);

% VMA: BEWARE, xmin, which contained the best fit (minimum std of the radiance temperatures)
% in the Braik version, now
% contains the best agent, i.e. best emissivity spectrum.
% The best fit is now in stdTr

% VMA Original He version in which Trv, the vector of radiance temperatures, was added to be plotted next
xmin = double(xmin);
for uu = 1:dim
    eval(['y' num2str(uu) '= xmin(' num2str(uu) ')',';']);
end
Tr = 0;
for uu = 1:dim
    eval(['f' num2str(uu) ' = m(uu) ./ (y' num2str(uu) '+ D(uu))',';']);
    Tr = Tr + eval(['f' num2str(uu)]);
    Trv(uu)=eval(['f' num2str(uu)]) ; % VMA : vector of radiance temperatures to be plotted next
end
Tr = Tr/dim;
T_fmincon=Tr;
E = exp(xmin);


% Tr = 0;
% for uu = 1:dim
%     Trv(uu)=m(uu)/(bestPosition(uu)+ D(uu));
%     Tr = Tr + Trv(uu);
% end
% Tr = Tr/dim;
% T_fmincon=Tr;
% E = exp(bestPosition);
% stdTr=xmin;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% VMA Plot of the results
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


if(avis_permitted ==1)
    figure(72);
    emiss_perm=V_lamda_i./V_b.*exp(C2./Lamda_i.*(1./T_fmincon-1./T_b));
    if(t_threshold>1.0001) % CSA-Corrected
        if(max_allowed_EMR<1.01)
            plot(lambdalist,emiss_perm, '-o', 'Color', [0.6, 0.6, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA  olive
        else
            plot(lambdalist,emiss_perm, '-o', 'Color', [0., 0.9, 0.9], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA cyan foncé
        end
    else  % CSA-buggy
        if(max_allowed_EMR<1.01)
            plot(lambdalist,emiss_perm, '-x', 'Color', [1, 0.647, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA orange
        else
            plot(lambdalist,emiss_perm, '-x', 'Color', [0.7, 0, 0.7], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA entre violet et magenta
        end
    end
end


disp([' Final radiance temperatures : ' ])
disp([' Mean radiance temperature : ' num2str(T_fmincon) ' K'])  % VMA;
disp([' stdTr : ' num2str(stdTr) ' K'])  % VMA;

figure(80);clf; % VMA
plot(lambdalist,Trv,'k');hold on; % VMA

    if(t_threshold>1.0001) % CSA-Corrected
        if(max_allowed_EMR<1.01)
            plot(lambdalist,T_fmincon*ones(1,dim), '-o', 'Color', [0.6, 0.6, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA  olive
        else
            plot(lambdalist,T_fmincon*ones(1,dim), '-o', 'Color', [0., 0.9, 0.9], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA cyan foncé
        end
    else  % CSA-buggy
        if(max_allowed_EMR<1.01)
            plot(lambdalist,T_fmincon*ones(1,dim), '-x', 'Color', [1, 0.647, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA orange
        else
            plot(lambdalist,T_fmincon*ones(1,dim), '-x', 'Color', [0.7, 0, 0.7], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA entre violet et magenta
        end
    end

% VMA: title({'Final radiance temperatures'},'interpreter','latex'); % VMA
xlabel('Wavelength (µm)'); % VMA
ylabel('Temperature (K)');  % VMA
if(avis_synth==1)
    text(0.52, (max(Trv)+min(Trv))/2+(max(Trv)-min(Trv))*0.35, sprintf('"True" temper.: %.1f K', T_grey),'FontName','Times','fontsize',10)
    rmse=sqrt(stdTr^2+(T_fmincon-T_grey)^2);
else
    text(0.52, (max(Trv)+min(Trv))/2+(max(Trv)-min(Trv))*0.35, sprintf('"True" temper.: %.1f K', T_unknown),'FontName','Times','fontsize',10)
    rmse=sqrt(stdTr^2+(T_fmincon-T_unknown)^2);
end
disp([' rmse : ' num2str(rmse) ' K'])  % VMA;
text(0.52, (max(Trv)+min(Trv))/2+(max(Trv)-min(Trv))*0.25, sprintf('Mean radiance temper.: %.1f K', T_fmincon),'FontName','Times','fontsize',10)
text(0.52, (max(Trv)+min(Trv))/2+(max(Trv)-min(Trv))*0.15, sprintf('St. dev.: %.2f K', stdTr),'FontName','Times','fontsize',10)
text(0.52, (max(Trv)+min(Trv))/2+(max(Trv)-min(Trv))*0.05, sprintf('RMSE: %.2f K', rmse),'FontName','Times','fontsize',10)
% text(0.5, T_fmincon+(max(Trv)-min(Trv))*0.15,...
%     sprintf('Pixel : ii= %.0f /41, jj= %.0f /37', ii,jj),'FontName','Times','fontsize',10)
% text(0.5, mean(Trv)+(max(Trv)-min(Trv))*0.15,...
%     sprintf('Pixel : ii=%.0f,  jj=%.0f ', ii,jj),'FontName','Times','fontsize',10)
grid on;
figure(80);

disp([' Calculated emissivity (ratio)  : ' ])  % VMA;
%E  % VMA  % Emissivity as obtained from Chameleon
E_from_T_fmincon=exp(C2/T_fmincon./lambdalist-D); % VMA : Emissivity as obtained from the Mean radiance temperature
%[E E_from_T_fmincon]
figure(81);clf; % VMA
plot(lambdalist,E,'k');hold on; % VMA
    if(t_threshold>1.0001) % CSA-Corrected
        if(max_allowed_EMR<1.01)
            plot(lambdalist,E_from_T_fmincon, '-o', 'Color', [0.6, 0.6, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA  olive
        else
            plot(lambdalist,E_from_T_fmincon, '-o', 'Color', [0., 0.9, 0.9], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA cyan foncé
        end
    else  % CSA-buggy
        if(max_allowed_EMR<1.01)
            plot(lambdalist,E_from_T_fmincon, '-x', 'Color', [1, 0.647, 0], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA orange
        else
            plot(lambdalist,E_from_T_fmincon, '-x', 'Color', [0.7, 0, 0.7], 'LineWidth', 1.5, 'MarkerSize', 7);hold on;  % VMA entre violet et magenta
        end
    end

% VMA title({'Final emissivity (ratio) spectrum'},'interpreter','latex'); % VMA  (ratio) was added
xlabel('Wavelength (µm)'); % VMA
if(avis_bb_ceram==0)
    ylabel('Emissivity ratio (-)');  % VMA
    text(0.52, (max(E_from_T_fmincon)+min(E_from_T_fmincon))/2+(max(E_from_T_fmincon)-min(E_from_T_fmincon))*0.2, sprintf('Emissivity-ratio limits : [%.2f, %.1f] ', min_allowed_EMR,max_allowed_EMR),'FontName','Times','fontsize',10)
elseif(avis_bb_ceram==1)
    ylabel('Emissivity (-)');  % VMA
    text(0.52, (max(E_from_T_fmincon)+min(E_from_T_fmincon))/2+(max(E_from_T_fmincon)-min(E_from_T_fmincon))*0.2, sprintf('Emissivity limits : [%.2f, %.1f] ', min_allowed_EMR,max_allowed_EMR),'FontName','Times','fontsize',10)
end
% text(0.5, mean(E_from_T_fmincon)+(max(E_from_T_fmincon)-min(E_from_T_fmincon))*0.15,...
%     sprintf('Pixel : ii=%.0f,  jj=%.0f ', ii,jj),'FontName','Times','fontsize',10)
grid on;
figure(81);

E_mean=exp(mean(log(E)));
E_from_T_fmincon_mean=exp(mean(log(E_from_T_fmincon)));
allowed_EMR_mean=sqrt(min_allowed_EMR*max_allowed_EMR);
differ_E_pc=(E_mean/allowed_EMR_mean-1)*100;
differ_E_from_T_fmincon_pc=(E_from_T_fmincon_mean/allowed_EMR_mean-1)*100;

disp(['Comparison with the mean value sqrt(min_allowed_EMR*max_allowed_EMR :  ' num2str(allowed_EMR_mean)])
disp(['of the mean value of emiss ratio given by CSA : ' num2str(E_mean)  ',  ' num2str(differ_E_pc) ' % '])
disp(['of the mean value of emiss ratio obtained by mean(Tr) given by CSA : ' num2str(E_from_T_fmincon_mean)  ',  ' num2str(differ_E_from_T_fmincon_pc) ' % '])
%pause % VMA

