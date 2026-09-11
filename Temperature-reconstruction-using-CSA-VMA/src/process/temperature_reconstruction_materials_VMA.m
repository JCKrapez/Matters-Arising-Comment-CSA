% temperature_reconstruction_materials_VMA.m

% The initial version of this script is
% temperature_reconstruction_materials.m
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
% - modifications to enable the display of intermediate and additional results
% - introduction of global variables to share with other
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

% REMARK 2:
% There is a flipud when displaying the image of the calculated temperature
% zone 1 center:[row 6, col 19](flame), which was at top center, is displayed at bottom center.
% zone 2 center:[row 20, col 6], which was at mid-height left, keeps roughly the same position.
% zone 3 center:[row 21, col 32], which was at mid-height right, keeps roughly the same position.
% zone 4 center:[row 37, col 18] (far from flame), which was at bottom center, is displayed at top center.
% There is an  LEFT-RIGHT FLIP ERROR with respect to fig. 6b in the original
% paper by He et al.
% In the code zone 2 is kept left, while it is represented at right in fig
% 6b.


% clc  % VMA
clear all
% close all
format long   % VMA

% VMA: Introduction of global variables that are set here and used in other
% functions or scripts
global ii jj T_unknown iTC
global avis_plot
global avis_permitted avis_synth avis_disp_emiss avis_bb_ceram
global Nagents Nbiter min_allowed_EMR max_allowed_EMR
global stdTr  % VMA added variable containing the estimation error (fmin0) which was surprisingly overwritten by the last instruction in the He version of Chameleon.
global avis_seed seed %  VMA seed used to initiate the random numbers
global t_threshold  % VMA  introduction of a threshold to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
global  T_perm_min T_perm_max
global avis_thalweg i1_thalweg i2_thalweg


disp([' To perform the analysis, several options are available and many parameters can be changed at will.']);
disp([' The selection by the user is made by answering to a series of questions; suggested values are in parentheses, the default value is in square brackets.']);
disp([' Just press ENTER to select the default value. ']);

% VMA avis_bb_ceram: global flag to separate the case BLACKBODY from the case CERAMIC in spectrumTemperature5_5_VMA
% avis_bb_ceram = input('Enter 0 if blackbody, 1 if ceramic (0-1) ? [1] : ');
% if isempty(avis_bb_ceram) avis_bb_ceram = 1;end
avis_bb_ceram=1;  %  0 if blackbody, 1 if ceramic

% formal assignements (not used) to allow calling the function
% spectrumTemperature5_5_VMA when no use of these parameter is made therein
Nagents=5;
Nbiter=10;
max_allowed_EMR=1;
min_allowed_EMR=0.1;
inter=2;

disp('ADDED ANALYSIS: ');
disp('Display of the surface of the cost function ');

avis_thalweg = input('Select 1(YES) or 0(NO) to plot the surface of the cost function and its thalweg ? [0] : ');
if isempty(avis_thalweg) avis_thalweg = 0;end
if(avis_thalweg==1)
    i1_thalweg = input('Select the number of the 1st wavelength ? [1] : ');
    if isempty(i1_thalweg) i1_thalweg = 1;end
    i2_thalweg = input('Select the number of the 2nd wavelength ? [21] : ');
    if isempty(i2_thalweg) i2_thalweg = 21;end
    
elseif(avis_thalweg==0)
    
    disp('ADDED ANALYSIS: ');
    disp('Possible analysis of the raw signals V and Vb (fig. 9 in Matters Arising), their ratio V/Vb, and the emissivity derived thereof, given a temperature (fig. 10 in Matters Arising) ');
    avis_disp_emiss = input('Enter 0 (NO), or 1 (YES, with the assumed true temperature Ttrue), 2 (YES, with Ttrue and Ttrue +/- 20K) ? [0] : ');
    if isempty(avis_disp_emiss) avis_disp_emiss = 0;end
    
    % VMA Preparation of some figures
    if(avis_disp_emiss>0)
        if(ishandle(73))
            figure(70);figure(71);figure(73);
            avis_erase = input('Allow the following results to be superimposed on existing Figures 70,71,73 (0) or clear the previous results from this figure (1)?  [0] : ');
            if isempty(avis_erase) avis_erase = 0;end
        end
        if((ishandle(73) && avis_erase==1) || ~ishandle(73))
            for ifig=[70 71 73]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-70)+10 730-340 500 270]);
            end
        end
        
    elseif(avis_disp_emiss==0)
        if(ishandle(62))
            figure(62);figure(63);
            avis_erase = input('Allow the following results to be superimposed on existing Figures 62,63 (0) or clear the previous results from these figures (1)?  [0] : ');
            if isempty(avis_erase) avis_erase = 0;end
        end
        if((ishandle(62) && avis_erase==1) || ~ishandle(62))
            for ifig=[62 63]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-60)+10 730 500 270]);  % left bottom width height  [450*(ifig-1)+10 670 scrsz(3)/3.5 300]
            end
        end
        
        
        if(ishandle(72))
            figure(72);
            avis_erase = input('Allow the following results to be superimposed on existing Figure 72 (0) or clear the previous results from these figures (1)?  [0] : ');
            if isempty(avis_erase) avis_erase = 0;end
        end
        if((ishandle(72) && avis_erase==1) || ~ishandle(72))
            for ifig=[72]
                figure(ifig);clf;
                h=gcf;set(h,'Position',[500*(ifig-70)+10 730-340 500 270]);  % left bottom width height  [450*(ifig-1)+10 670 scrsz(3)/3.5 300]
            end
        end
        
        % VMA Two types of analyses with CSA
        disp('There are two options for applying CSA. The input signal ratio V/Vb can come from two sources: ');
        disp('1- either that simulated from a hypothetic material with constant emissivity (GREY MATERIAL) at a given temperature');
        disp('2- or that measured experimentally by He et al. on an alumina plate (fig. 9-18 in Matters Arising)');
        avis_synth = input('Enter the selection (1-2) ? [2] : ');
        if isempty(avis_synth) avis_synth = 2;end
        
        disp('ADDED ANALYSIS: ');
        avis_permitted = input('Select 1(YES) or 0(NO) for a computation and display of the permitted solutions (fig. 13, 17 in Matters Arising) ? [1] : ');
        if isempty(avis_permitted) avis_permitted = 1;end
        
        
        % VMA introduction of a seed for the random number generation
        disp('CSA uses random numbers');
        disp('To ensure reproducibility for the random number generation, a seed is required (no seed is used in the original code by He at al.)');
        disp('There are three options: ');
        disp(' 1- no seed is used (as in the original code by He at al.)');
        disp(' 2- the same seed is used for each CSA call (fig. 11-18 in Matters Arising)');
        disp(' 3- the seed is incremented by 1 before each CSA call ');
        avis_seed = input('Enter the seed option  ? [2] : ');
        if isempty(avis_seed) avis_seed = 2;end
        if(avis_seed==2)
            disp(' The constant value of the seed was 2 for the results in fig. 11-13, 15-17');
            disp(' For the results related to the four zones in fig. 14 and 18, it was first uniformy set to 1, then to 2, and finally to 3');
            
            seed = input('Enter the constant value of the seed (integer number) [2] : ');
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
        disp(' - 15 in the original code by He at al. (fig. 11-14 in Matters Arising) ');
        disp(' - 30 (fig. 15-18 in Matters Arising)');
        Nagents = input('Enter the number of agents ? [15] : ');
        if isempty(Nagents) Nagents = 15;end
        
        % VMA Selection of the number of iterations
        disp(' The number of iterations was: ');
        disp(' - 30 in the original code by He at al. (fig. 11-14 in Matters Arising)');
        disp(' - 300 (fig. 15-18 in Matters Arising) ');
        Nbiter = input('Enter the number of iterations ? [30] : ');
        if isempty(Nbiter) Nbiter = 30;end
        
        % VMA Selection of the minimum allowed emissivity (ratio)
        disp(' The minimum allowed emissivity (ratio) was: ');
        disp(' - 0.1 in the original code by He at al. (fig. 11-14 in Matters Arising)');
        disp(' - 0.01 (fig. 15-18 in Matters Arising) ');
        min_allowed_EMR = input('Enter the minimum allowed emissivity ratio ? [0.1] : ');
        if isempty(min_allowed_EMR) min_allowed_EMR = 0.1;end
        
        % VMA Selection of the maximum allowed emissivity (ratio)
        disp(' The maximum allowed emissivity (ratio) was: ');
        disp(' - 0.99999 in the original code by He at al. and in Matters Arising (for all tests excluding those on a blackbody)');
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
    % Al2O3 experimental Case and pixel. Therefore the following specifications are arbitrary:
    case_num = 2;
    avis_option_zone = -1;
else
    
    % VMA Selection of the experimental case
    disp(' There are 3 cases for Al2O3: ');
    disp(' - Case 1 was only presented in the supplement (Case 1 was selected in the distributed version of the script temperature_reconstruction_materials.m)  ');
    disp(' - Case 2 was described in the main paper ');
    disp(' - Case 3 was only presented in the supplement ');
    case_num = input(' Select a case number ? [2] : ');
    if isempty(case_num) case_num = 2;end
    
    % VMA Which zone should be considered (e.g. full zone or only the
    % central pixel) ?
    disp(' There are several options for performing the analysis in space');
    if(avis_thalweg==1 || avis_disp_emiss>=1)
        disp(' - option -X (X=1,2,3 or 4) : the analysis is performed only at the center of zone X');
        avis_option_zone = input('Enter the option to consider ? [-1] : ');
        if isempty(avis_option_zone) avis_option_zone = -1;end
    else
        disp(' - option 0 : the analysis is performed for all pixels (nominal)');
        disp(' - option X (X=1,2,3 or 4) : the analysis is performed only over the zone X and the mean value is calculated for it (zone 1 is the closest to the flame)');
        disp(' - option 1234 : the analysis is performed only over the four zones 1,2,3,4 and the mean value is calculated for each');
        disp(' - option -X (X=1,2,3 or 4) : the analysis is performed only at the center of zone X');
        disp(' - option -1234 : the analysis is performed only at the center of all four zones 1,2,3,4');
        avis_option_zone = input('Enter the option to consider ? [-1] : ');
        if isempty(avis_option_zone) avis_option_zone = -1;end
        if(avis_option_zone>0)
            disp(' The default number of pixels on each side of the center point in each zone is 2');
            inter = input('Enter the number of pixels on each side of the center point to consider (0-1-2) [2] : ');
            if isempty(inter) inter = 2;end
        end
    end
end
% VMA Preparation of the remaining figures
for ifig=[60 61]
    figure(ifig);clf;
    h=gcf;set(h,'Position',[500*(ifig-60)+10 730 500 270]);
end

for ifig=[80 81 82]
    figure(ifig);clf;
    h=gcf;set(h,'Position',[500*(ifig-80)+10 730-680 500 270]);
end


temref = 1300; % blackbody reference temperature in Celsius

load("..\..\data\materials\al2o3\thermocouple.mat") % Temperature from the thermocouples
temperature_true_matrix = fliplr(T_truecase);  % VMA : in T_truecase, the max temperature is in the 1st column
% VMA : in T_truecase the TC readings 1, 2, 3, 4 are in columns 1, 2, 3, 4  (the max TC temperature is in the 1st column)
% VMA : in temperature_true_matrix, the TC readings 1, 2, 3, 4 are in columns 4, 3, 2, 1 (the max TC temperature is in the 4th
% column)
% VMA: the variable temperature is not used, hence  temperature_true_matrix and hole_index and pixel_belong_mask are useless


load_path = "..\..\data\materials\al2o3\Intensity_case"+num2str(case_num)+".mat";
load(load_path);

Tdistri = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
StdTrdistri = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
StdTrdistri_pc = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
Bias = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
Bias_pc = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
RMSE = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37
RMSE_pc = zeros(size(calibrated_inten_matrix, 1), size(calibrated_inten_matrix, 2));   % VMA:  size : 41x37


%%
addpath '..\temperature_retrieval'
addpath '..\temperature_retrieval\CSA'

T_true  = T_truecase(case_num, :) + 273;  % VMA This instruction has been moved upwards for convenience
% VMA : in T_truecase the TC readings 1, 2, 3, 4 are in columns 1, 2, 3, 4  (the max TC temperature is in the 1st column)


% VMA
% Next lines are used to define the set of scanned pixels of the image
if(avis_option_zone==-1234)
    indices= fliplr(center_int);   % center_int contains the coordinates (column, line) of the center of the four zones, which explains the use of fliplr
elseif(avis_option_zone<0 && avis_option_zone >=-4)
    indices= fliplr(center_int(-avis_option_zone,:));
elseif(avis_option_zone==0)
    % [I, J] = ndgrid(1:avis_option_zone:size(Tdistri, 1), 1:avis_option_zone:size(Tdistri, 2));
    [I, J] = ndgrid(1:size(Tdistri, 1), 1:size(Tdistri, 2));
    indices = [I(:), J(:)];
elseif(avis_option_zone>0 && avis_option_zone <=4)
    [I, J] = ndgrid(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    indices = [I(:), J(:)];
elseif(avis_option_zone==1234)
    indices=[];
    for sami=1:4
        [I, J] = ndgrid(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        indices = [indices;[I(:), J(:)]];
    end
end


% VMA
% for ii = 1:size(Tdistri, 1)  % VMA first loop for ii (He version)
%     for jj = 1:size(Tdistri, 2)  % VMA second loop for jj (He version)

% VMA
% New version with a single loop for the coordinates of the considered pixels
for k = 1:size(indices, 1)
    ii = indices(k, 1);
    jj = indices(k, 2);
    
    % Fig. 6 has:
    % index 1 at bottom center
    % index 2 at mid-height right
    % index 3 at mid-height left
    % index 4 at top center
    % VMA : description of the content of pixel_belong_mask:
    % index 4 is in 5x5 cells at bottom center
    % index 2 is in 5x5 cells at mid-height left
    % index 3 is both in 5x5 cells at mid-height right and in the 5x5 cells at top center
    % There is essentially a flipud and fliplr between Fig. 6 and
    % pixel_belong_mask (if we take into account the error that in the zone at top center of pixel_belong_mask the index should have been 1, not 3).
    % Recall that in temperature_true_matrix, the TC readings 1, 2, 3, 4 are in columns 4, 3, 2, 1
    % The following 6 lines make that:
    % where hole_index = 4 (5x5 cells at bottom center of pixel_belong_mask) the matrix "temperature" contains TC 1
    % where hole_index = 2 (5x5 cells at mid-height left of pixel_belong_mask) the matrix "temperature" contains TC 3
    % where hole_index = 3 (5x5 cells at mid-height right and in the 5x5 cells at top center of pixel_belong_mask) the matrix "temperature" contains TC 2
    % where hole_index = 0 (everywhere else)  the matrix "temperature" contains TC 4
    % Which is false.
    % Hopefully, there is no consequence since the matrix "temperature" is not used !...
    
    hole_index = pixel_belong_mask(ii, jj); % VMA temperature is not used; hence hole_index is useless
    if hole_index > 0
        temperature = temperature_true_matrix(case_num, hole_index); % VMA temperature is not used; hence hole_index is useless
    else
        temperature = temperature_true_matrix(case_num, 1);  % VMA temperature is useless; hence hole_index is useless
    end
    % VMA: temperature is not used, hence  temperature_true_matrix and hole_index and   pixel_belong_mask are useless
    
    intencasecali  = calibrated_inten_matrix{ii, jj}.intencasecali;
    intenblackcali = calibrated_inten_matrix{ii, jj}.intenblackcali;
    
    if(avis_seed==3)  % VMA: seed incremented
        seed=seed+1;
    end
    
    if( avis_option_zone~=0)  % VMA: Search for the TC number for a pixel in zones 1 to 4 (not available for avis_option_zone=0)
        iTC=find(abs(ii-center_int(:,2))<=inter & abs(jj-center_int(:,1))<=inter);  % VMA : number of closest TC
        T_unknown=T_true(iTC);  % VMA : T_unknown is the "true" unknown temperature used in spectrumTemperature5-5_VMA for plotting the "true" emissivity and "true" permitted solution
        % Do not make the confusion with  Tunknown which will be the solution given
        % by CSA in two lines, after run temrecon_singlepoint_VMA.m.
    end
    
    
    %run temrecon_singlepoint.m   % VMA Original call of script
    run temrecon_singlepoint_VMA.m   % VMA call of modified script
    % VMA temrecon_singlepoint  % VMA Different call
    
    if(avis_thalweg ==1 || avis_disp_emiss>0)
        return  % VMA The code stops after these accessory analyses (plot of the cost function surface or the raw signals and inferred emissivity)
    end
    
    Tdistri(ii, jj) = Tunknown;
    
    StdTrdistri(ii, jj) = stdTr;  % VMA New matrix with the distribution of the error (fit residual, i.e. std of the set of radiance temperatures)
    
    disp(['Pixel : ' num2str(ii) '/' num2str(size(Tdistri, 1))  ', ' num2str(jj) '/' num2str(size(Tdistri, 2)) ': ' num2str(Tunknown)  ' +/- ' num2str(stdTr)  ]); % VMA
    if( avis_option_zone~=0)  % VMA: Calculation of the difference with the corresponding TC reading
        StdTrdistri_pc(ii, jj) = stdTr/T_true(iTC)*100;
        Bias(ii, jj)=Tunknown-T_true(iTC);
        Bias_pc(ii, jj)=Bias(ii, jj)/T_true(iTC)*100;
        RMSE(ii, jj)=sqrt(stdTr^2+Bias(ii, jj)^2);
        RMSE_pc(ii, jj)=RMSE(ii, jj)/T_true(iTC)*100;
        disp(['The difference with the TC temperature ' num2str(T_true(iTC)) ' is ' num2str(Bias(ii, jj)) ' i.e. ' num2str(Bias_pc(ii, jj))  ' %']); % VMA
        disp(['The RMSE is ' num2str(RMSE(ii, jj)) ' i.e.  ' num2str(RMSE_pc(ii, jj)) ' %']); % VMA
        %pause % VMA
    end
end   % VMA  Now one single "for" loop for k
% end  In the He version there were two "for" loops, i.e. for ii and jj  %VMA

Tsample = zeros(1, 4);
Std_sample = zeros(1, 4);
Std_sample_pc = zeros(1, 4);
Bias_sample = zeros(1, 4);
Bias_sample_pc = zeros(1, 4);
RMSE_sample = zeros(1, 4);
RMSE_sample_pc = zeros(1, 4);

if( avis_option_zone == 1234) % VMA: He version where Tsample is obtained by calculating the mean in each four zones
    for sami = 1:4
        a1 = Tdistri(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a2 = StdTrdistri(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a3 = StdTrdistri_pc(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a4 = Bias(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a5 = Bias_pc(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a6 = RMSE(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        a7 = RMSE_pc(center_int(sami, 2)-inter:center_int(sami, 2)+inter, center_int(sami, 1)-inter:center_int(sami, 1)+inter);
        Tsample(sami) = mean(a1(:));
        Std_sample(sami)=sqrt(mean(a2(:).*a2(:)+(a1(:)-Tsample(sami)).*(a1(:)-Tsample(sami))));
        Std_sample_pc(sami)=sqrt(mean(a3(:).*a3(:)+(a1(:)-Tsample(sami)).*(a1(:)-Tsample(sami))/T_true(iTC)/T_true(iTC)));
        Bias_sample(sami)=mean(a4(:));
        Bias_sample_pc(sami)=mean(a5(:));
        RMSE_sample(sami)=sqrt(mean(a6(:).*a6(:)));
        RMSE_sample_pc(sami)=sqrt(mean(a7(:).*a7(:)));
    end
elseif( avis_option_zone == -1234) % VMA: Only at the center of all 4 zones
    for sami = 1:4
        Tsample(sami)=Tdistri(center_int(sami, 2),center_int(sami, 1));
        Std_sample(sami) = StdTrdistri(center_int(sami, 2),center_int(sami, 1));
        Std_sample_pc(sami) = StdTrdistri_pc(center_int(sami, 2),center_int(sami, 1));
        Bias_sample(sami) = Bias(center_int(sami, 2),center_int(sami, 1));
        Bias_sample_pc(sami) = Bias_pc(center_int(sami, 2),center_int(sami, 1));
        RMSE_sample(sami) = RMSE(center_int(sami, 2),center_int(sami, 1));
        RMSE_sample_pc(sami) = RMSE_pc(center_int(sami, 2),center_int(sami, 1));
    end
elseif( avis_option_zone<0 && avis_option_zone >=-4)  % VMA: Only at the center of zone number -avis_option_zone
    Tsample(-avis_option_zone)=Tdistri(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    Std_sample(-avis_option_zone)=StdTrdistri(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    Std_sample_pc(-avis_option_zone)=StdTrdistri_pc(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    Bias_sample(-avis_option_zone)=Bias(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    Bias_sample_pc(-avis_option_zone)=Bias_pc(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    RMSE_sample(-avis_option_zone)=RMSE(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
    RMSE_sample_pc(-avis_option_zone)=RMSE_pc(center_int(-avis_option_zone, 2),center_int(-avis_option_zone, 1));
elseif( avis_option_zone>0 && avis_option_zone <=4)  % VMA: The mean will be calculated for the zone number avis_option_zone
    a1 = Tdistri(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a2 = StdTrdistri(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a3 = StdTrdistri_pc(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a4 = Bias(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a5 = Bias_pc(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a6 = RMSE(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    a7 = RMSE_pc(center_int(avis_option_zone, 2)-inter:center_int(avis_option_zone, 2)+inter, center_int(avis_option_zone, 1)-inter:center_int(avis_option_zone, 1)+inter);
    Tsample(avis_option_zone) = mean(a1(:));
    Std_sample(avis_option_zone)=sqrt(mean(a2(:).*a2(:)+(a1(:)-Tsample(avis_option_zone)).*(a1(:)-Tsample(avis_option_zone))));
    Std_sample_pc(avis_option_zone)=sqrt(mean(a3(:).*a3(:)+(a1(:)-Tsample(avis_option_zone)).*(a1(:)-Tsample(avis_option_zone))/T_true(iTC)/T_true(iTC)));
    Bias_sample(avis_option_zone)=mean(a4(:));
    Bias_sample_pc(avis_option_zone)=mean(a5(:));
    RMSE_sample(avis_option_zone)=sqrt(mean(a6(:).*a6(:)));
    RMSE_sample_pc(avis_option_zone)=sqrt(mean(a7(:).*a7(:)));
end

if( avis_option_zone == 0) % VMA: Start of 2D plot of He version where all pixels are considered
    disp('Different options for padding/filtering:');
    disp(' 1- "replicate" padding by 1 then median filtering on a [7,7] neighborhood (original method by He et al. Padding by 1 is not enough');
    disp(' 2- "symmetric" padding by 3 then median filtering on a [7,7] neighborhood ');
    disp(' 3- median filtering on a [7,7] neighborhood by specifying "symmetric" padding (equivalent to option 2) ');
    avis_pad = input('Enter the option (1-2-3)? [1] : ');
    if isempty(avis_pad) avis_pad = 1;end
    
    if(avis_pad==1)
        Tdistrimed_prov = medfilt2(padarray(Tdistri, [1, 1], "replicate"), [7, 7]);
        Tdistrimed = Tdistrimed_prov(2:end-1, 2:end-1);
    elseif(avis_pad==2)
        % Padding by 1 on each side is not enough. We must pad by 3.
        Tdistrimed_prov = medfilt2(padarray(Tdistri, [3, 3], "symmetric"), [7, 7]);
        Tdistrimed = Tdistrimed_prov(4:end-3, 4:end-3);
    elseif(avis_pad==3)
        % Or apply the 'symmetric' padding options (default padding is with 0s)
        Tdistrimed = medfilt2(Tdistri, [7, 7],'symmetric');
    end
    figure(1);clf;
    imagesc(flipud(Tdistrimed)); colorbar; axis equal tight;
    colormap(flipud(hot));
    % clim([1220 1440]);  % VMA does not work with R2021a version
    caxis([1230 1430]);   % VMA does work with R2021a version
    if(t_threshold==1)
        caxis([1250 1500]);   % VMA does work with R2021a version
    elseif(t_threshold>1)
        caxis([1400 1650]);   % VMA does work with R2021a version
    end
    xticks([]); yticks([]);
    title(['Al2O3 Case ' num2str(case_num) ' Temperature Distribution (filtered) (K)']);
    impixelinfo
    
    % VMA We also plot the raw results (not only 7x7 mean filtering as in the original paper) :
    %   figure('Name', ['AlO Case ' num2str(case_num) 'raw results'])
    figure(2);clf;
    imagesc(flipud(Tdistri(1:end,:))); colorbar; axis equal tight;
    colormap(flipud(hot));
    % clim([1220 1440]);  % VMA does not work with R2021a version
    if(t_threshold==1)
        caxis([1250 1500]);   % VMA does work with R2021a version
    elseif(t_threshold>1)
        caxis([1400 1650]);   % VMA does work with R2021a version
    end
    xticks([]); yticks([]);
    title(['Al2O3 Case ' num2str(case_num) ' Temperature Distribution (raw) (K)']);
    impixelinfo
    
    % VMA Plot of the 2D distribution of the temperature std
    figure(3);clf;
    imagesc(flipud(StdTrdistri(1:end,:))); colorbar; axis equal tight;
    colormap(flipud(hot));
    % clim([1220 1440]);  % VMA does not work with R2021a version
    if(t_threshold==1)
        caxis([0 80]);   % VMA does work with R2021a version
    elseif(t_threshold>1)
        caxis([0 10]);   % VMA does work with R2021a version
    end
    xticks([]); yticks([]);
    title(['Al2O3 Case ' num2str(case_num) ' Std Distribution (K)']);
    impixelinfo
    
end  % VMA: End of 2D plot of He version where all pixels are considered

if( avis_option_zone ~= 0 && avis_synth==2)
    figure(62);%clf;
    avis_true = input('Plot of the true temperature (thermocouple reading) (0-1) ? [0] : ');
    if isempty(avis_true) avis_true = 0;end
    if(avis_true==1)
        plot(1:4, T_true,  's', 'Color', [0, 0.8, 0] ,'LineWidth', 2,'MarkerSize', 8); hold on
    end
    if(t_threshold==1)
        symbc='x';
    elseif(t_threshold>1)
        symbc='o';
    end
    avis_color = input('Choose the color (1-orange, 2-violet, 3-dark blue, 4-olive green)  ? [1] : ');
    if isempty(avis_color) avis_color = 1;end
    if(avis_color==1)
        V_Color=[1, 0.647, 0]; %orange
    elseif(avis_color==2)
        V_Color=[0.7, 0, 0.7]; %violet
    elseif(avis_color==3)
        V_Color=[0, 0, 0.8]; %dark blue
    elseif(avis_color==4)
        V_Color=[0.6, 0.6, 0]; %olive green
    end
    
    if(avis_option_zone == -1234)
        errorbar([1 2 3 4]+0.07*seed, Tsample, Std_sample, symbc, 'Color', V_Color , 'LineWidth', 1.5,'MarkerSize', 8); hold on
    else
        errorbar([1 2 3 4], Tsample, Std_sample, symbc, 'Color', V_Color , 'LineWidth', 1.5,'MarkerSize', 8); hold on
    end
    grid on; xticks([1 2 3 4]); xlabel('Position');
    ylabel('Temperature (K)');
    %legend('Thermocouple','Reconstruction')
    if(t_threshold==1)
        ylim([1250 1450])
    elseif(t_threshold>1)
        ylim([1350 1600])
    end
    figure(62);
    
    figure(63);%clf;
    %plot(1:4, Tsample, 'ro', 'LineWidth', 2); hold on
    %plot([1 2 3 4]+0.07*seed, Bias_sample_pc, 'k', symbc, 'MarkerSize', 8); hold on
    %plot([1 2 3 4]+0.07*seed, Std_sample_pc, 'r', symbc,'MarkerSize', 8); hold on
    %plot([1 2 3 4]+0.07*seed, RMSE_sample_pc, 'm',symbc,'MarkerSize', 8); hold on
    if(avis_option_zone == -1234)
        plot([1 2 3 4]+0.07*seed, RMSE_sample_pc, symbc,'Color', V_Color ,'LineWidth', 1.5,'MarkerSize', 8); hold on
    else
        plot([1 2 3 4], RMSE_sample_pc, symbc,'Color', V_Color ,'LineWidth', 1.5,'MarkerSize', 8); hold on
    end
    grid on; xticks([1 2 3 4]); xlabel('Position'); ylabel('RMS Error (%)');
    %legend('Mean bias', 'Mean std', 'RMSE')
    %     ymin=min([min(Bias_sample_pc) min(RMSE_sample_pc)]);
    %     ymax=max([max(Bias_sample_pc) max(RMSE_sample_pc)]);
    %     ylim([ymin-(ymax-ymin)*0.1 ymax+(ymax-ymin)*0.1]);
    if(t_threshold==1)
        ylim([0 6]);yticks(0:1:6);
    elseif(t_threshold>1)
        ylim([0 16]);yticks(0:2:16);
    end
    figure(63);
    
    figure(62);
    figure(70);
    figure(71);
    figure(72);
    figure(73);
    
end