% temperature_reconstruction_blackbody_cali_VMA.m

% The initial version of this script is
% temperature_reconstruction_blackbody_cali.m
% It was obtained from https://github.com/KirinShi/Temperature-reconstruction-using-CSA
% He, Y., Chen, M. K., Huang, M. et al. Dispersive Meta-lens Thermometry for High-temperature Measurements.
% Nature Communications 16, 10090 (2025). https://doi.org/10.1038/s41467-025-65171-7

% The analysis of the paper, next the code and data made available by the
% authors, led to the submission to Nature Communications of a document "Matters Arising"
% The results obtained with the present version of the code are in the Supplementary material

% The changes made to the code were made by:
% J.-C. Krapez
% 2026-09-07
% DOTA, ONERA, 13300 Salon de Provence, France
% jean-claude.krapez@onera.fr

% VMA is the flag used to highlight the changes (VMA for Version Matters
% Arising):
% - modifications to enable display of intermediate or additional results
% - introduction of global variables to share them with other
% functions and scripts
% - correction of a bug in Chameleon.m  -> Chameleon_VMA.m

% REMARK 1:
% In temperature_reconstruction_blackbody_cali.m the wavelengths
% are decreasing whereas in temperature_reconstruction_materials.m
% they are increasing.
% However, Chameleon.m depends on the order.
% For consistency, it would have been better to use always the same
% order. The results in Matters Arising were obtained by using the original arrangements (i.e. wavelengths in opposite order).
% The dependency on the order in Chameleon_VMA.m was detected
% afterwards by comparing the results for simulated GREY materials
% when using
% either temperature_reconstruction_blackbody_cali_VMA.m or
% temperature_reconstruction_materials_VMA.m. The temperature results
% are slightly different (by 1 to 1.5K). The comparison between
% temperature_reconstruction_blackbody_cali_VMA.m and
% temperature_reconstruction_materials_VMA.m is relevant only for
% simulated GREY materials. Therefore, for consistency with the
% results published in Matters Arising for GREY materials (which were obtained with
% temperature_reconstruction_materials_VMA.m), it was decided to
% reverse the order of wavelenghts when
% spectrumTemperature5_5_VMA.m is called by temperature_reconstruction_blackbody_cali_VMA.m, but only in the case of
% simulation of GREY materials
% (this text is reproduced below, when spectrumTemperature5_5_VMA is
% called)


% clc  % VMA
clear all
% close all
format long % VMA

% VMA Introduction of global variables that are set here and used in other
% functions or scripts
global ii jj T_unknown
global avis_plot
global avis_permitted avis_synth avis_disp_emiss avis_bb_ceram
global Nagents Nbiter min_allowed_EMR max_allowed_EMR
global stdTr  % VMA added variable containing the estimation error (fmin0) which was surprisingly overwritten by the last instruction in the He version of Chameleon.
global avis_seed seed % VMA seed used to initiate the random numbers
global t_threshold  % VMA  introduction of a threshold to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
global  T_perm_min T_perm_max
global avis_thalweg i1_thalweg i2_thalweg

disp([' To perform the analysis, several options are available and many parameters can be changed at will.']);
disp([' The selection by the user is made by answering to a series of questions; suggested values are in parentheses, the default value is in square brackets.']);
disp([' Just press ENTER to select the default value. ']);

% VMA avis_bb_ceram: global flag to separate the case BLACKBODY from the case CERAMIC in spectrumTemperature5_5_VMA
% avis_bb_ceram = input('Enter 0 if blackbody, 1 if ceramic (0-1) ? [1] : ');
% if isempty(avis_bb_ceram) avis_bb_ceram = 1;end
avis_bb_ceram=0;  %  0 if blackbody, 1 if ceramic

% formal assignements (not used) to allow calling the function
% spectrumTemperature5_5_VMA when no use of these parameter is made therein
Nagents=5;
Nbiter=10;
max_allowed_EMR=1;
min_allowed_EMR=0.1;



disp('ADDED ANALYSIS: ');
disp('Display of the surface of the cost function ');

avis_thalweg = input('Select 1(YES) or 0(NO) to plot the surface of the cost function and its thalweg (fig. 1 in Matters Arising)? [0] : ');
if isempty(avis_thalweg) avis_thalweg = 0;end
if(avis_thalweg==1)
    i1_thalweg = input('Select the number of the 1st wavelength ? [1] : ');
    if isempty(i1_thalweg) i1_thalweg = 1;end
    i2_thalweg = input('Select the number of the 2nd wavelength ? [21] : ');
    if isempty(i2_thalweg) i2_thalweg = 21;end
    
elseif(avis_thalweg==0)
    
    disp('ADDED ANALYSIS: ');
    disp('Possible analysis of the raw signals V and Vb, their ratio V/Vb, and the emissivity derived thereof, given a temperature (fig. 2 in Matters Arising)');
    avis_disp_emiss = input('Enter 0 (NO), or 1 (YES, with the assumed true temperature Ttrue), 2 (YES, with Ttrue and Ttrue +/- 20K) ? [0] : ');
    if isempty(avis_disp_emiss) avis_disp_emiss = 0;end
    
    % VMA Preparation of some figures
    if(avis_disp_emiss>0)
        if(ishandle(73))
            figure(73);
            % VMA Preparation of the figures
            avis_erase = input('Allow the following results to be superimposed on existing Figures 73 (0) or clear the previous results from this figure (1)?  [0] : ');
            if isempty(avis_erase) avis_erase = 0;end
        end
        if((ishandle(73) && avis_erase==1) || ~ishandle(73))
            for ifig=[73]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-70)+10 730-340 500 270]);
            end
        end
        
    elseif(avis_disp_emiss==0)
        if(ishandle(62))
            figure(62);
            figure(72);
            avis_erase = input('Allow the following results to be superimposed on existing Figures 62 and 72 (0) or clear the previous results from these figures (1)?  [0] : ');
            if isempty(avis_erase) avis_erase = 0;end
        end
        if((ishandle(62) && avis_erase==1) || ~ishandle(62))
            for ifig=[62]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-60)+10 730 500 270]);  % left bottom width height  [450*(ifig-1)+10 670 scrsz(3)/3.5 300]
            end
            for ifig=[72]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-70)+10 730-340 500 270]);  % left bottom width height  [450*(ifig-1)+10 670 scrsz(3)/3.5 300]
            end
        end
        
        % VMA Two types of analyses with CSA
        disp('There are two options for applying CSA. The input signal ratio V/Vb can come from two sources: ');
        disp('1- either that simulated from a hypothetic material with constant emissivity (GREY MATERIAL) at a given temperature');
        disp('2- or that measured experimentally by He et al. on a blackbody at one or several temperatures between 1673K and 1873K');
        avis_synth = input('Enter the selection (1-2) ? [1] : ');
        if isempty(avis_synth) avis_synth = 1;end
        
        disp('ADDED ANALYSIS: ');
        avis_permitted = input('Select 1(YES) or 0(NO) for a computation and display of the permitted solutions ? [1] : ');
        if isempty(avis_permitted) avis_permitted = 1;end
        
        
        % VMA introduction of a seed for the random number generation
        disp('CSA uses random numbers');
        disp('To ensure reproducibility for the random number generation, a seed is required (no seed is used in the original code by He at al.)');
        disp('There are three options: ');
        disp(' 1- no seed is used (as in the original code by He at al.)');
        disp(' 2- the same seed is used for each CSA call');
        disp(' 3- the seed is incremented by 1 before each CSA call (as in Matters Arising fig. 3,4,5,6,C1,C2,C3,C4 starting from seed=2) ');
        avis_seed = input('Enter the seed option  ? [3] : ');
        if isempty(avis_seed) avis_seed = 3;end
        if(avis_seed==2)
            seed = input('Enter the constant value of the seed (integer number) ? [2] : ');
            if isempty(seed) seed = 2;end
        elseif(avis_seed==3)
            seed = input('Enter the initial value of the seed (integer number) ? [2] : ');
            if isempty(seed) seed = 2;end
            seed=seed-1;  % later it will be incremented by 1 before each CSA call
        end
        
        % VMA bug fix or not in Chameleon.m ->  Chameleon_VMA.m
        disp(' The velocity updates need the accelation at the denominator, which depends on the interation number t. The denominator is 0 for t=1 (in Braik and He versions) ');
        disp(' To solve this bug, t has been replaced by max(t,t_threshold) in the formula for the acceleration.');
        disp(' Select 1 for t_threshold to perform the (buggy) calculations (as with the Braik and He codes)');
        disp(' Select a value higher than 1 for t_threshold to fix the bug (1.1 was used in Matters Arising; it should be not less than 1.01) ');
        t_threshold = input('Enter a value for t_threshold ? [1.1] : ');
        if isempty(t_threshold) t_threshold = 1.1;end
        
        
        
        % VMA Selection of the number of agents (chameleons)
        disp(' The number of agents (chameleons) was: ');
        disp(' - 5 in the original code by He at al. (orange and violet crosses in fig. 3,C1,C2,C3,C4 in Matters Arising) ');
        disp(' - 30 (blue circles in fig. 3,C2,C4; curves with crosses and circles in fig. 4,5,6)');
        Nagents = input('Enter the number of agents ? [5] : ');
        if isempty(Nagents) Nagents = 5;end
        
        % VMA Selection of the number of iterations
        disp(' The number of iterations was: ');
        disp(' - 100 in the original code by He at al. (orange and violet crosses in fig. 3,C1,C2,C3,C4 in Matters Arising)');
        disp(' - 300 (blue circles in fig. 3,C2,C4; curves with crosses and circles in fig. 4,5,6) ');
        Nbiter = input('Enter the number of iterations ? [100] : ');
        if isempty(Nbiter) Nbiter = 100;end
        
        % VMA Selection of the minimum allowed emissivity (ratio)
        min_allowed_EMR = input('Enter the minimum allowed emissivity ratio (0.1 in the original code and in Matters Arising) ? [0.1] : ');
        if isempty(min_allowed_EMR) min_allowed_EMR = 0.1;end
        
        % VMA Selection of the maximum allowed emissivity (ratio)
        disp(' The maximum allowed emissivity ratio was: ');
        disp(' - 0.99999 in the original code by He at al. (orange and violet crosses in fig. 3,C1,C2,C3,C4 in Matters Arising)');
        disp(' - 1.2 in Matters Arising for the blackbody tests (violet crosses and blue circles in fig. 3)');
        disp(' - 1 in Matters Arising for the tests on simulated GREY materials (curves with crosses and circles in fig. 4,5,6)');
        
        max_allowed_EMR = input('Enter the maximum allowed emissivity ratio ? [0.99999] : ');
        if isempty(max_allowed_EMR) max_allowed_EMR = 0.99999;end
        
        
        % VMA How many intermediate plots of the Chameleon coordinates (i.e.
        % emissivity) should be displayed ?
        disp(' Intermediate results of the "Chameleon emissivity" can be displayed with a varying levels of detail, depending on a display index:');
        disp(' 1 for no plots ');
        disp(' 2 for a plot after initialization');
        disp(' 3 for a plot after exploration (6)');
        disp(' 5 for a plot after Rotation (30)');
        disp(' 7 for a plot after velocity update (210)');
        disp(' 11 for a plot after handling of boundary violations (2310)');
        avis_plot = input('Enter the display index (product of the first prime numbers; 2310 for all intermediate results)? [1] : ');
        if isempty(avis_plot) avis_plot = 1;end
    end
end

if(~isempty(avis_synth) && avis_synth==1)
    % By construction, the "synthetic" case does not depend on the
    % blackbody temperature and pixel. Therefore the following specifications are arbitrary:
    avis_case = 1;
    avis_pixel = 173;
else
    % VMA Selection of the Blackbody temperature or set of temperatures
    disp(' There are several options for performing the analysis in temperature');
    disp(' Notice that Case 11 corresponds to the reference (1773 K), Case 1 is for the highest temperature, 1873K, Case 21 is for the lowest temperature, 1673K');
    if(avis_thalweg==1 || avis_disp_emiss>=1)
        disp(' - option X (X from 1 to 21) : the analysis is performed only for Case X');
        avis_case = input('Enter the option to consider ? [1] : ');
        if isempty(avis_case) avis_case = 1;end
    else
        disp(' - option 0 : the analysis is performed for all 21 Cases, from 1873 K to 1673 K (default)');
        disp(' - option X (X from 1 to 21) : the analysis is performed only for Case X');
        disp(' - option -1 : the analysis is performed only over the two limiting Cases, 1 and 21');
        avis_case = input('Enter the option to consider ? [0] : ');
        if isempty(avis_case) avis_case = 0;end
    end
    
    % VMA Which pixel should be considered in the image or should the whole image be taken into account
    disp(' There are several options for performing the analysis in space');
    %disp(' global_num is equal to the product of maximum values of global_x=23 and global_y=15');
    if(avis_thalweg==1 || avis_disp_emiss>=1)
        disp(' - option N (N from 1 to 345) : the analysis is performed only for pixel of coordinate x (to right) ,y (y upwards) with N=(x-1)*15+(17-y)=(C-1)*15+L  ');
        disp(' - option 173 for the central pixel ');
        disp(' - option 79 for the central pixel in the top left quadrant ');
        disp(' - option 244 for the central pixel in the top right quadrant ');
        disp(' - option 87 for the central pixel in the bottom left quadrant ');
        disp(' - option 252 for the central pixel in the bottom right quadrant ');
        avis_pixel = input(' Select the option to consider ? [173] : ');
        if isempty(avis_pixel) avis_pixel = 173;end
    else
        disp(' - option 0 : the analysis is performed for all 15 lines X 23 columns =345 pixels (nominal)');
        disp(' - option N (N from 1 to 345) : the analysis is performed only for pixel of coordinate x (to right) ,y (y upwards) with N=(x-1)*15+(17-y)=(C-1)*15+L  ');
        disp(' - option 173 for the central pixel ');
        disp(' - option 79 for the central pixel in the top left quadrant ');
        disp(' - option 244 for the central pixel in the top right quadrant ');
        disp(' - option 87 for the central pixel in the bottom left quadrant ');
        disp(' - option 252 for the central pixel in the bottom right quadrant ');
        avis_pixel = input(' Select the option to consider ? [173] : ');
        if isempty(avis_pixel) avis_pixel = 173;end
    end
end


% constant
e = 2.718281828;
C2 = 1.4388*10^4;    % um*K

% VMA: initial settings for ub and lb are now selected interactively by the user :
% ub = 0.99999;
% lb = 0.1;
ub = max_allowed_EMR;
lb = min_allowed_EMR;

addpath '..\temperature_retrieval'
addpath '..\temperature_retrieval\CSA'

global_x = 23;
global_y = 15;
global_num = global_x * global_y;
case_num = 21;
case_ref = 11;
T_ref = 1773;
T_step = -10;
T_start = T_ref - (case_ref - 1) * T_step;
T_end = T_ref + (case_num - case_ref) * T_step;
Temperature_faceresults = zeros(global_num, case_num);
Error_faceresults = zeros(global_num, case_num);
stdTr_faceresults = zeros(global_num, case_num);
stdTr_pc_faceresults = zeros(global_num, case_num);
Bias_faceresults = zeros(global_num, case_num);
Bias_pc_faceresults = zeros(global_num, case_num);
RMSE_faceresults = zeros(global_num, case_num);
RMSE_pc_faceresults = zeros(global_num, case_num);

load("..\..\data\blackbody\lambdalistt.mat")

if(avis_pixel==0)
    domain=1:global_num;
elseif(avis_pixel>0)
    domain=avis_pixel;
end

if(avis_case==0)
    case_test=1:1:case_num;  %VMA  Nominal use
elseif(avis_case>0 && avis_case<=21)
    case_test=avis_case;
elseif(avis_case==-1)
    case_test=[1 case_num] ;
end


for ifig=[60 61 63]
    figure(ifig);clf;
    h=gcf;set(h,'Position',[500*(ifig-60)+10 730 500 270]);
end

for ifig=[70 71]
    figure(ifig);clf;
    h=gcf;set(h,'Position',[500*(ifig-70)+10 730-340 500 270]);
end
for ifig=[80 81 82]
    figure(ifig);clf;
    h=gcf;set(h,'Position',[500*(ifig-80)+10 730-680 500 270]);
end



%for uu = 1:1:case_num  % VMA: first loop in the He version (in true
%temperature)
for uu = case_test  % VMA: new first loop (in BB temperature)
    %for globali = 1:global_num  % VMA: second loop in the He version (in
    %space)
    for globali = domain  % VMA: new second loop (in space)
        jj=floor(globali / global_y)+1;    % VMA : used for plots in functions
        ii = mod(globali, global_y);       % VMA : used for plots in functions
        disp(num2str(uu)+"  "+num2str(globali));
        load("..\..\data\blackbody\" + ...
            case_ref + "\" + globali + ...
            "\calibratedfinal_" + num2str(T_ref) + "_nor.mat");
        cdataRep = calibratedfinal_unknown_nor'; % reference
        
        T_unknown = T_start + (uu - 1) * T_step;
        load("..\..\data\blackbody\" + ...
            uu + "\" + globali + ...
            "\calibratedfinal_" + num2str(T_unknown) + "_nor.mat");
        cdataRep_1 = calibratedfinal_unknown_nor'; % target
        
        if(avis_seed==3)  % VMA: seed incremented
            seed=seed+1;
        end
        
        if(avis_synth==1)
            [R_temperature,x] = spectrumTemperature5_5_VMA(cdataRep(1,end:-1:1), cdataRep_1(1,end:-1:1),T_ref,lambdalistt(end:-1:1),1,1,Nagents,Nbiter,ub,lb);
            % In temperature_reconstruction_blackbody_cali.m the wavelengths
            % are decreasing whereas in temperature_reconstruction_materials.m
            % they are increasing.
            % However, Chameleon.m depends on the order (probably because of the "rotations").
            % For consistency, it would have been better to use always the same
            % order. The results in Matters Arising were obtained by using the original arrangements (i.e. wavelengths in opposite order).
            % The dependency on the order in Chameleon_VMA.m was detected
            % afterwards by comparing the results for simulated GREY materials
            % when using
            % either temperature_reconstruction_blackbody_cali_VMA.m or
            % temperature_reconstruction_materials_VMA.m. The temperature results
            % are slightly different (by 1 to 1.5K). The comparison between
            % temperature_reconstruction_blackbody_cali_VMA.m and
            % temperature_reconstruction_materials_VMA.m is relevant only for
            % simulated GREY materials. Therefore, for consistency with the
            % results published in Matters Arising for GREY materials (which were obtained with
            % temperature_reconstruction_materials_VMA.m), it was decided to
            % reverse the order of wavelenghts when
            % spectrumTemperature5_5_VMA.m is called by temperature_reconstruction_blackbody_cali_VMA.m, but only in the case of
            % simulation of GREY materials
            %
        else
            %[R_temperature,x] =
            %spectrumTemperature5_5(cdataRep(1,:),cdataRep_1(1,:),T_ref,lambdalistt,1,1,5,100,ub,lb);   % VMA original call of spectrumTemperature5_5 in the code of He et al.
            [R_temperature,x] = spectrumTemperature5_5_VMA(cdataRep(1,:), cdataRep_1(1,:),T_ref,lambdalistt,1,1,Nagents,Nbiter,ub,lb);
        end
        if(avis_thalweg ==1 || avis_disp_emiss>=1)  % VMA plot of the cost function surface or the raw signals and inferred emissivity
            return
        end
        Error = abs(R_temperature-T_unknown)/(T_unknown) * 100;
        % VMA: He et al. version used the absolute value of the relative difference (in %) with respect to the true temperature
        Bias = R_temperature-T_unknown;  % VMA: Bias with respect to the true temperature
        Bias_pc = Bias/T_unknown*100  ;  % VMA: Bias (in %) with respect to the true temperature
        
        RMSE = sqrt(stdTr^2+Bias^2);     % VMA: the RMSE was added to the results
        RMSE_pc = RMSE/T_unknown*100;     % VMA:  RMSE (in %)
        Temperature_faceresults(globali,uu) = R_temperature;
        Error_faceresults(globali,uu) = Error;
        stdTr_faceresults(globali,uu) = stdTr;
        stdTr_pc_faceresults(globali,uu) = stdTr/T_unknown*100;     % VMA:  Std of radiance temperature (in %)
        Bias_faceresults(globali,uu) = Bias;
        Bias_pc_faceresults(globali,uu) = Bias_pc;
        RMSE_faceresults(globali,uu) = RMSE;
        RMSE_pc_faceresults(globali,uu) = RMSE_pc;
        
    end
    %pause
end



if(avis_pixel==0) % VMA: Plot of the distributions if the analysis was made for all pixels
    
    % Temperature
    % for ti = 1:1:case_num  % VMA:  loop in the He version
    for ei = case_test  % VMA:  new loop
        aaa = reshape(Temperature_faceresults(:,ei),[global_y,global_x]);
        figure(ei)
        imagesc(aaa(:,1:end-1));
        if(t_threshold<1.001 && max_allowed_EMR<1.001)
            caxis([1873-20, 1873+20]-(ei-1)*10)
        elseif(t_threshold<1.001 && max_allowed_EMR>1.001)
            caxis([1873-30, 1873+30]-20-(ei-1)*10)
        elseif(t_threshold>1.001)
            caxis([1873-30, 1873+30]+130-(ei-1)*10)
        end
        colormap(flipud(hot));
        colorbar
        axis image;
        x_ticks = linspace(0.5, 22.5, 3);
        y_ticks = linspace(0.5, 15.5, 3);
        set(gca, 'XTick', x_ticks, 'YTick', y_ticks);
        axis image;
        impixelinfo;
    end
    %end
    
    
    % Error (i.e. absolute bias in % as in the paper by He et al.)
    % for ei = 1:1:case_num   % VMA:  loop in the He version
    for ei = case_test   % VMA:  new loop
        aaa = reshape(Error_faceresults(:,ei),[global_y,global_x]);
        figure(ei+case_num)
        imagesc(aaa(:,1:end-1));
        caxis([0, 1])
        x_ticks = linspace(0.5, 22.5, 3);
        y_ticks = linspace(0.5, 15.5, 3 );
        colormap jet
        colorbar
        set(gca, 'XTick', x_ticks, 'YTick', y_ticks);
        axis image;
        impixelinfo;
    end
    
    % RMSE_pc  RMSE in %
    % for ei = 1:1:case_num   % VMA:  loop in the He version
    for ei = case_test   % VMA:  new loop
        aaa = reshape(RMSE_pc_faceresults(:,ei),[global_y,global_x]);
        figure(ei+2*case_num)
        imagesc(aaa(:,1:end-1));
        caxis([0, 1])
        x_ticks = linspace(0.5, 22.5, 3);
        y_ticks = linspace(0.5, 15.5, 3 );
        colormap jet
        colorbar
        set(gca, 'XTick', x_ticks, 'YTick', y_ticks);
        axis image;
        impixelinfo;
    end
end

avgtem = mean(Temperature_faceresults(domain,:),1); % VMA row vector where each element is the average of the corresponding column (average in space)
avgerr = mean(Error_faceresults(domain,:),1)';
avgstdTr = mean(stdTr_faceresults(domain,:),1)';
avgstdTr_pc = mean(stdTr_pc_faceresults(domain,:),1)';
avgbias = mean(Bias_faceresults(domain,:),1)';
avgbias_pc = mean(Bias_pc_faceresults(domain,:),1)';
avgRMSE = mean(RMSE_faceresults(domain,:),1)';
avgRMSE_pc = mean(RMSE_pc_faceresults(domain,:),1)';



T_truth= [T_start:T_step:T_start+T_step*(case_num-1)];

%figure  % VMA
if(avis_synth==2)
    figure(62);
    %clf; % VMA
    % plot(T_truth(1:10),T_truth(1:10),'--k');hold on
    % errorbar(T_truth(1:10), avgtem(1:10), avgstdTr(1:10), 'go', 'LineWidth', 1.); hold on
    % plot(T_truth(12:end),T_truth(12:end),'--k');hold on
    % errorbar(T_truth(12:end), avgtem(12:end), avgstdTr(12:end), 'go', 'LineWidth', 1.); hold on
    
    plot(T_truth(1:10),T_truth(1:10),'--k');hold on
    plot(T_truth(12:end),T_truth(12:end),'--k');hold on
    
    if(t_threshold>1.0001) % CSA-Corrected
        if(max_allowed_EMR<1.01)
            errorbar(T_truth(1:10), avgtem(1:10), avgstdTr(1:10), 'o', 'Color', [0.6, 0.6, 0], 'LineWidth', 1.);hold on;  % VMA olive
            errorbar(T_truth(12:end), avgtem(12:end), avgstdTr(12:end), 'o', 'Color', [0.6, 0.6, 0], 'LineWidth', 1.);hold on;  % VMA olive
        else
            errorbar(T_truth(1:10), avgtem(1:10), avgstdTr(1:10), 'o', 'Color', [0., 0.9, 0.9], 'LineWidth', 1.);hold on;  % VMA cyan foncé
            errorbar(T_truth(12:end), avgtem(12:end), avgstdTr(12:end), 'o', 'Color', [0., 0.9, 0.9], 'LineWidth', 1.);hold on;  % VMA cyan foncé
        end
    else  % CSA-buggy
        if(max_allowed_EMR<1.01)
            errorbar(T_truth(1:10), avgtem(1:10), avgstdTr(1:10), 'x', 'Color', [1, 0.647, 0], 'LineWidth', 1., 'MarkerSize', 7);hold on;   % VMA orange
            errorbar(T_truth(12:end), avgtem(12:end), avgstdTr(12:end), 'x', 'Color', [1, 0.647, 0], 'LineWidth', 1., 'MarkerSize', 7);hold on;  % VMA orange
        else
            errorbar(T_truth(1:10), avgtem(1:10), avgstdTr(1:10), 'x', 'Color', [0.7, 0, 0.7], 'LineWidth', 1., 'MarkerSize', 7);hold on;   % VMA entre violet et magenta
            errorbar(T_truth(12:end), avgtem(12:end), avgstdTr(12:end), 'x', 'Color', [0.7, 0, 0.7], 'LineWidth', 1., 'MarkerSize', 7);hold on;  % VMA entre violet et magenta
        end
    end
    
    grid on; xlabel('True temperature (K)'); ylabel('Temperature (K)');
    xlim([1650 1900]);
    ylim([1650 2050]);
    %legend('True T', 'Estimated T');
    
    
    %figure  % VMA
    figure(63);
    %clf; % VMA
    %plot(avgerr);ylim([0 2])
    %plot(T_truth,avgerr,'om');hold on;
    plot(T_truth(1:10),avgbias_pc(1:10),'xk');hold on;
    plot(T_truth(1:10),avgstdTr_pc(1:10),'ob');hold on;
    plot(T_truth(1:10),avgRMSE_pc(1:10),'or');hold on;
    plot(T_truth(12:end),avgbias_pc(12:end),'xk');hold on;
    plot(T_truth(12:end),avgstdTr_pc(12:end),'ob');hold on;
    plot(T_truth(12:end),avgRMSE_pc(12:end),'or');hold on;
    xlim([1650 1900]);
    ue=max([max(avgbias_pc) max(avgstdTr_pc) max(avgRMSE_pc)]);
    le=min([min(avgbias_pc) min(avgstdTr_pc) min(avgRMSE_pc)]);
    ylim([le-(ue-le)*0.1 ue+(ue-le)*0.1 ]);
    grid on; xlabel('True temperature (K)'); ylabel('Error (%)');
    %legend('abs (bias)', 'signed bias', 'Std(Tr)', 'RMSE')
    legend('signed bias', 'Std(Tr)', 'RMSE')
    disp(' Results regarding the AVERAGE errors over all  temperatures ');
    disp([ 'Bias (K)     : ' num2str(mean(avgbias)) ' +/- ' num2str(std(avgbias)) ])
    disp([ 'std(Tr) (K)  : ' num2str(mean(avgstdTr)) ' +/- ' num2str(std(avgstdTr)) ])
    disp([ 'RMSE (K)     : ' num2str(mean(avgRMSE)) ' +/- ' num2str(std(avgRMSE)) ])
    disp([ 'Bias (%)     : ' num2str(mean(avgbias_pc)) ' +/- ' num2str(std(avgbias_pc)) ])
    disp([ 'std(Tr) (%)  : ' num2str(mean(avgstdTr_pc)) ' +/- ' num2str(std(avgstdTr_pc)) ])
    disp([ 'RMSE (%)     : ' num2str(mean(avgRMSE_pc)) ' +/- ' num2str(std(avgRMSE_pc)) ])
    
    
    disp(' Results regarding the MAXIMUM errors over all  temperatures ');
    disp([ 'Bias (K)     : ' num2str(max(abs(avgbias)))  ])
    disp([ 'std(Tr) (K)  : ' num2str(max(avgstdTr))  ])
    disp([ 'RMSE (K)     : ' num2str(max(avgRMSE))  ])
    disp([ 'Bias (%)     : ' num2str(max(abs(avgbias_pc)))  ])
    disp([ 'std(Tr) (%)  : ' num2str(max(avgstdTr_pc))  ])
    disp([ 'RMSE (%)     : ' num2str(max(avgRMSE_pc))  ])
    
end
