% --------------------------------------------------------------
% SSTDR Full System Simulation Using a Systematic Solution Procedure
% By: Hunter Ellis, Naveen K. T. Jayakumar, Mashad Uddin Saleh, Joel B. Harley, and
% Cynthia Furse
%
% This code provides a method of simulating the SSTDR response of systems
% with lengths of transmission lines separating loads that are in series and
% parallel
%
% For more information on the concepts behind this code, please see: 
% A Model for SSTDR signal propagation Through PV Strings
% Hunter Ellis, Naveen K. T. Jayakumar, Samuel Kingston, Evan Benoit, Mashhad Uddin Saleh, Michael A. Scapula, Joel B. Harley, Cynthia Furse
% IEEE Sensors Journal
%
% This material is based upon work supported by the U.S. Department
% of Energy’s Office of Energy Efficiency and Renewable Energy (EERE)
% under Solar Energy Technologies Office (SETO) Agreement Number
% DE-EE0008169. Also, we want to express gratitude to Gardner Energy
% (West Haven, UT) for their support with testing on their PV modules. 
% --------------------------------------------------------------
 
% Get the directory of this script file
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '../data/element_data'));

% Ensure results directory exists relative to this script's location
results_dir = fullfile(script_dir, '..', 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_dir = [results_dir filesep]; % Add trailing separator for convenience

%==========================================================================
% First simulation: single 100 ohm resistor with RG58 transmission line
%==========================================================================
num_z = 1;                      %number of impednaces to simulate 
imped = ones(40960, 1).*100;    %impedance of elements (100 Ohms)
len = [14];               %distance between each elements
s_p = [];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[single_resistor_sig, single_resistor_dist] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation
 
%==========================================================================
% second simulation: two 100 ohm resistors in parallel with RG58 transmission line
%==========================================================================
num_z = 2;                      %number of impednaces to simulate 
imped = ones(40960, 2).*100;    %impedance of elements (100 Ohms)
len = [14,1.83];               %distance between each elements
s_p = ['p'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[two_resistor_sig_parallel, two_resistor_dist_parallel] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation
 
%==========================================================================
% third simulation: three 100 ohm resistors in parallel with RG58 transmission line
%==========================================================================
num_z = 3;                      %number of impednaces to simulate 
imped = ones(40960, 3).*100;    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43];               %distance between each elements
s_p = ['p','p'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[three_resistor_sig_parallel, three_resistor_dist_parallel] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

%==========================================================================
% fourth simulation: four 100 ohm resistors in parallel with RG58 transmission line
%==========================================================================
num_z = 4;                      %number of impednaces to simulate 
imped = ones(40960, 4).*100;    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['p','p','p'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_resistor_sig_parallel, four_resistor_dist_parallel] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

 %==========================================================================
% plot parallel resistor simulations
%==========================================================================
 figure(1)
 
 subplot(2,2,1)
 plot(single_resistor_dist,single_resistor_sig)
 xlim([8,23])
 
  subplot(2,2,2)
 plot(two_resistor_dist_parallel,two_resistor_sig_parallel)
 xlim([8,23])
 
   subplot(2,2,3)
 plot(three_resistor_dist_parallel,three_resistor_sig_parallel)
 xlim([8,23])
 
   subplot(2,2,4)
 plot(four_resistor_dist_parallel,four_resistor_sig_parallel)
 xlim([8,23])
 
 saveas(gcf, [results_dir 'Parallel_Resistor.png'])
 
 
 %==========================================================================
% simulation: two 100 ohm resistors in series with RG58 transmission line
%==========================================================================
num_z = 2;                      %number of impednaces to simulate 
imped = ones(40960, 2).*100;    %impedance of elements (100 Ohms)
len = [14,1.83];               %distance between each elements
s_p = ['s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[two_resistor_sig_series, two_resistor_dist_series] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation
 
  %==========================================================================
% third simulation: three 100 ohm resistors in series with RG58 transmission line
%==========================================================================
num_z = 3;                      %number of impednaces to simulate 
imped = ones(40960, 3).*100;    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43];               %distance between each elements
s_p = ['s','s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[three_resistor_sig_series, three_resistor_dist_series] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

 %==========================================================================
% fourth simulation: four 100 ohm resistors in series with RG58 transmission line
%==========================================================================
num_z = 4;                      %number of impednaces to simulate 
imped = ones(40960, 4).*100;    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['s','s','s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_resistor_sig_series, four_resistor_dist_series] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

 %==========================================================================
% plot series resistor simulations
%==========================================================================
 figure(2)
 
 subplot(2,2,1)
 plot(single_resistor_dist,single_resistor_sig)
 xlim([8,23])
 
  subplot(2,2,2)
 plot(two_resistor_dist_series,two_resistor_sig_series)
 xlim([8,23])
 
   subplot(2,2,3)
 plot(three_resistor_dist_series,three_resistor_sig_series)
 xlim([8,23])
 
   subplot(2,2,4)
 plot(four_resistor_dist_series,four_resistor_sig_series)
 xlim([8,23])
 
 saveas(gcf, [results_dir 'Series_Resistor.png'])
 
%==========================================================================
% simulation: four 5 pF Capacitors in parallel with RG58 transmission line
%==========================================================================

% define impedance of capacitors
Fsampl = 4*24e6;        %sampling rate at 48 MHz
Fmodu  = Fsampl/4;      % Modulation rate
Fchip  = Fsampl/4;      % Chip rate
PNL    = 10240;         % Number of chips
Q = PNL*Fsampl/(Fchip);         % Length of simulated signal
n = 1:Q;                        % Sample axis
t = n/Fsampl;                   % Time axis
r = floor(Q/2)+1;               % Center sample in frequency
f = ifftshift((n-r)/Q)*Fsampl;  % Frequency axis
f(1) = 1;                       % avoid discontinuity at 0 Hz
imp_cap = (1./(1j*2*pi.*f.*5*(10^-12))).'; %impedance of a 5 pf capacitor

num_z = 4;                      %number of impednaces to simulate 
imped = [imp_cap, imp_cap, imp_cap, imp_cap];    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['p','p','p'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_cap_sig_parallel, four_cap_dist_parallel] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation
 
  %==========================================================================
% simulation: four 5 pF Capacitors in series with RG58 transmission line
%==========================================================================
num_z = 4;                      %number of impednaces to simulate 
imped = [imp_cap, imp_cap, imp_cap, imp_cap];    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['s','s','s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_cap_sig_series, four_cap_dist_series] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

 %==========================================================================
% plot 5 pF capacitor simulations
%==========================================================================
 figure(3)
 
 subplot(2,1,1)
 plot(four_cap_dist_parallel,four_cap_sig_parallel)
  xlim([16,28])
 
  subplot(2,1,2)
 plot(four_cap_dist_series,four_cap_sig_series)
xlim([10,20])
 
 saveas(gcf, [results_dir '5_pF_capacitors.png'])
 
%==========================================================================
% simulation: four 220 pF Capacitors in parallel with RG58 transmission line
%==========================================================================
% determin impedance of capacitor from measured data
file = importdata('220PF.TXT');
imp_cap = impedance(Q,f,file);

num_z = 4;                      %number of impednaces to simulate 
imped = [imp_cap, imp_cap, imp_cap, imp_cap];    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['p','p','p'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_cap220_sig_parallel, four_cap220_dist_parallel] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation
 
  %==========================================================================
% simulation: four 220 pF Capacitors in series with RG58 transmission line
%==========================================================================
num_z = 4;                      %number of impednaces to simulate 
imped = [imp_cap, imp_cap, imp_cap, imp_cap];    %impedance of elements (100 Ohms)
len = [14,1.83, 2.43, 3.05];     %distance between each elements
s_p = ['s','s','s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 
[four_cap220_sig_series, four_cap220_dist_series] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns); %run simulation

 %==========================================================================
% plot 220 pF capacitor simulations plots
%==========================================================================
 figure(4)
 
 subplot(2,1,1)
 plot(four_cap220_dist_parallel,four_cap220_sig_parallel)
 xlim([9,33])
 
  subplot(2,1,2)
 plot(four_cap220_dist_series,four_cap220_sig_series)
 xlim([9,33])
 
 saveas(gcf, [results_dir '220_pF_capacitors.png'])
 
%==========================================================================
% simulate system with PV cable as the transmission line
%==========================================================================
% determin impedance of cells
file = importdata('CELL_GAT.TXT');
imp_cap = impedance(Q,f,file);

num_z = 4;                      %number of impednaces to simulate 
imped = [imp_cap, imp_cap, imp_cap, imp_cap];    %impedance of elements (100 Ohms)
len = [15.24,.30, .30, .30];     %distance between each elements %% some distances can lead to an artifact in the data
s_p = ['s','s','s'];            % configuration of each element (series, parallel, series)
Ns = 200;                       %number of points to output 

[ref_sig_cell, dist_m_cell] = SSTDR_SSP_PV_cable(num_z, imped, len, s_p, Ns); 

%==========================================================================
% plot 220 pF capacitor simulations plots
%==========================================================================
figure(5)
 plot(dist_m_cell, ref_sig_cell);
  xlim([10,35])
 saveas(gcf, [results_dir 'cells.png'])



