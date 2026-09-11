% temrecon_singlepoint_VMA.m

% The initial version of this function is
% temrecon_singlepoint.m
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

% VMA Global variables added
global Nagents Nbiter min_allowed_EMR max_allowed_EMR  % VMA

intenref = intenblackcali;
intenunknown = intencasecali;
lambdalist = (0.45:0.01:0.65);
Ttrue = temperature+273;  % VMA temperature and Ttrue are actually not used
C2=1.4338*10^(4);
Tunknown = 0;

% VMA Original instruction :
%  [Tunknown,emissivity] = spectrumTemperature5_5(intenref,intenunknown,temref+273, ...
%      lambdalist,1,1,15,500,0.99999,0.1);   % VMA Original  instruction

% VMA Modified instruction with global variables:
[Tunknown,emissivity] = spectrumTemperature5_5_VMA(intenref,intenunknown,temref+273, ...
   lambdalist,1,1,Nagents,Nbiter,max_allowed_EMR,min_allowed_EMR);   % VMA New  instruction


errorlist = abs(Tunknown - Ttrue)/Ttrue*100;  % VMA: Ttrue and errorlist are actually not used