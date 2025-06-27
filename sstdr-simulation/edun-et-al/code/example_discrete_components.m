% --------------------------------------------------------------
% SSTDR LUMPED ELEMENT REFLECTION DEMONSTRATION
% By: Ayobami S. Edun and Joel B. Harley
%
% This code demonstrates the theory for the reflection coefficient of a
% lumped element in symmetric and asymmetric transmission lines with an
% spread spectrum time domain reflectometry (SSTDR) signal. The theory is
% compared with experimentally acquired data. 
%
% For more information on the concepts behind this code, please see: 
% Spread Spectrum Time Domain Reflectometry with Lumped Elements on Asymmetric Transmission Lines
% Ayobami S. Edun, Naveen K. T. Jayakumar, Samuel Kingston, Cynthia Furse, Michael Scarpulla, Joel B. Harley
% IEEE Sensors Journal
%
% This material is based upon work supported by the U.S. Department
% of Energy’s Office of Energy Efficiency and Renewable Energy (EERE)
% under Solar Energy Technologies Office (SETO) Agreement Number
% DE-EE0008169. Also, we want to express gratitude to Gardner Energy
% (West Haven, UT) for their support with testing on their PV modules. 
% --------------------------------------------------------------
clear;

% Get the directory of this script file
script_dir = fileparts(mfilename('fullpath'));

% Add paths relative to the script location
addpath(fullfile(script_dir, '../data/capacitor_data/asymmetric'));
addpath(fullfile(script_dir, '../data/capacitor_data/symmetric'));
addpath(fullfile(script_dir, '../data/Parallel measurement/parallel_c_measurments'));
addpath(fullfile(script_dir, '../data/Parallel measurement/parallel_resistors'));
addpath(fullfile(script_dir, '../data/pv_cable'));
addpath(fullfile(script_dir, 'functions'));

% Check if data directory exists
data_dir = fullfile(script_dir, '../data');
if ~exist(data_dir, 'dir')
    error(['Data directory not found: ' data_dir newline ...
           'This script expects experimental data files to be present in the data directory.' newline ...
           'Please see REPRODUCING.md for instructions on how to run this code with the required data.']);
end

% Ensure results directory exists
results_dir = fullfile(script_dir, '../results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% SSTDR PARAMETERS
Fsampl = 192e6;             % Sampling rate for Wilma
Fmodu  = Fsampl/4;          % Modulation rate
Fchip  = Fsampl/4;          % Chip rate
PNL    = 10240;             % Number of chips
vop    = 0.7;               % Velocity of propgation (relative to c)

% ANALYSIS PARAMETERS
sindx = 10240/2;            % Frequency index to perform analysis with

% PLOTTING PARAMETERS
fontsize = 9;               % Font size of plots

% -------------------------------------------------
% DEFINE IMPEDANCE PROPERTIES
% -------------------------------------------------

% TRANSMISSION LINE PARAMETERS
R0 = 0.08;                  % Transmission line resistance
C0 = 3.4546e-11;            % Transmission line capacitance 
G0 = 4.5602e-14;            % Transmission line conductance
L0 = 8.0518e-7;             % Transmission line inductance      

% -------------------------------------------------
% INITIALIZE SIMULATION AND EXPERIMENT VARIABLES
% -------------------------------------------------
files       = {};       % Data files
folders     = {};       % Data folder
baselines   = {};       % Baseline file names
data_center = [];       % Isolates the first reflection
names       = {};       % Data names
types       = {};       % Data type
n = 0;                  % Initialize

% -------------------------------------------------
% DEFINE SIMULATION AND EXPERIMENTAL PARAMETERS
% -------------------------------------------------
n = n + 1;
% ===========================
% SYMMETRIC RESISTORS (SIMULATED INTERFACE IMPEDANCES)
RI{n} = 0:1:2000;                       % Series interface resistance
RLI{n} = 1e15*ones(1,length(RI{n}));    % Loss interface resistance
CI{n} = 1e5*ones(1,length(RI{n}));      % Series capacitance
LI{n} = zeros(1,length(RI{n}));         % Series inductance
Zf{n} = 1;                              % Symmetric = 1, Asymmetric = 1/2

% SYMMETRIC RESISTORS (EXPERIMENTAL DATA)
names{n} = ' Sym. Res.';
types{n} = 'series';
baselines{n} = '91ft_dr_harley_resistor_tests_Za_SC_Zb_SC_Zc_OC.lws';
data_center(n,:) = 35:49; % Isolates the first reflection. location range of first reflection 
folders{n} = fullfile(script_dir, '../data/pv_cable/');
files{n} = { ...
         '91ft_dr_harley_resistor_tests_Za_10ohm_Zb_10ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_20ohm_Zb_20ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_30ohm_Zb_30ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_56ohm_Zb_56ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_130ohm_Zb_130ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_503ohm_Zb_503ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_1kohm_Zb_1kohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_2kohm_Zb_2kohm_Zc_OC.lws', ...
         };
RE{n} = [10 20 30 56 130 500 1000 2000];
CE{n} = 1e120*ones(1,8);            
LE{n} = zeros(1,8);                     
% ===========================

n = n + 1;
% ===========================
% ASYMMETRIC RESISTORS (SIMULATED INTERFACE IMPEDANCES)
RI{n} = 0:1:2000;                       % Series resistance
RLI{n} = 1e15*ones(1,length(RI{n}));    % Loss resistance
CI{n} = 1e5*ones(1,length(RI{n}));      % Series capacitance
LI{n} = zeros(1,length(RI{n}));         % Series inductance
Zf{n} = 1/2;                            % Symmetric = 1, Asymmetric = 1/2

% ASYMMETRIC RESISTORS (EXPERIMENTAL DATA)
names{n} = 'Asym. Res.';
types{n} = 'series';
baselines{n} = '91ft_dr_harley_resistor_tests_Za_SC_Zb_SC_Zc_OC.lws';
data_center(n,:) = 35:49;
folders{n} = fullfile(script_dir, '../data/pv_cable/');
files{n} = {...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_10ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_20ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_30ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_56ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_130ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_503ohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_1kohm_Zc_OC.lws', ...
         '91ft_dr_harley_resistor_tests_Za_SC_Zb_2kohm_Zc_OC.lws', ...
         };  
RE{n} = [10 20 30 56 130 500 1000 2000];
CE{n} = 1e120*ones(1,8);
LE{n} = zeros(1,8);
% ===========================

n = n + 1;
% ===========================
% PARALLEL RESISTORS (SIMULATED INTERFACE IMPEDANCES)
RI{n} =  0:1:2000;                      % Parallel resistance
RLI{n} = 1e15*ones(1,length(RI{n}));    % Loss resistance
CI{n} = 1e5*ones(1,length(RI{n}));      % Series capacitance
LI{n} = zeros(1,length(RI{n}));         % Series inductance
Zf{n} = 1;                              % Symmetric = 1, Asymmetric = 1/2

% PARALLEL RESISTORS (EXPERIMENTAL DATA)
names{n} = ' Par. Res.';
types{n} = 'parallel';
baselines{n} = 'parallel_c_Za_OC_Zc_OC.lws';
data_center(n,:) = 35:49;
folders{n} = fullfile(script_dir, '../data/Parallel measurement/parallel_resistors/');
files{n} = {...
         'parallel_R_measurements_Za_10ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_20ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_30ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_56ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_130ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_503ohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_1kohm_zc_OC.lws', ...
         'parallel_R_measurements_Za_2kohm_zc_OC.lws', ...
         };  
RE{n} = [10 20 30 56 130 503 1000 2000];
CE{n} = 1e120*ones(1,8);
LE{n} = zeros(1,8);
% ===========================

n = n + 1;
% ===========================
% SYMMETRIC CAPACITORS (SIMULATED INTERFACE IMPEDANCES)
CI{n} = 1e-9:-1e-12:1e-12;              % Series capacitance
RI{n} = 5e-3*ones(1,length(CI{n}));     % Series resistance
RLI{n} = 1e4*ones(1,length(CI{n}));     % Loss resistance
LI{n} = 1e-9*ones(1,length(CI{n}));     % Series inductance (parasitic inductance)
Zf{n} = 1;                              % Symmetric = 1, Asymmetric = 1/2

% SYMMETRIC CAPACITORS (EXPERIMENTAL DATA) 
names{n} = ' Sym. Cap.';
types{n} = 'series';
baselines{n} = '91ft_PV_cable_cap_tests_Za_SC_Zb_SC_Zc_OC.lws';
data_center(n,:) = 35:49;
folders{n} = fullfile(script_dir, '../data/capacitor_data/symmetric/');
files{n} = { ...
        '91ft_PV_cable_cap_tests_Za_1pf_Zb_1pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_2pf_Zb_2pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_5pf_Zb_5pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_8pf_Zb_8pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_15pf_Zb_15pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_20pf_Zb_20pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_25pf_Zb_25pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_47pf_Zb_47pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_100pf_Zb_100pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_330pf_Zb_330pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_470pf_Zb_470pf_Zc_OC.lws', ...
         };
RE{n} = zeros(1,11);
CE{n} = [ 1 2 5 8 15 20 25 47 100 330 470]*1e-12;
LE{n} = zeros(1,11);
% ===========================

n = n + 1;
% ===========================
% ASYMMETRIC CAPACITORS (SIMULATED INTERFACE IMPEDANCES) 
CI{n} = 1e-9:-1e-12:1e-12;              % Series capacitance
RI{n} = 5e-3*ones(1,length(CI{n}));     % Series resistance
RLI{n} = 1e4*ones(1,length(CI{n}));     % Loss resistance
LI{n} = 1e-9*ones(1,length(CI{n}));     % Series inductance (parasitic inductance)
Zf{n} = 1/2;                            % Symmetric = 1, Asymmetric = 1/2

% ASYMMETRIC CAPACITORS (EXPERIMENTAL DATA) 
names{n} = 'Asym. Cap.';
types{n} = 'series';
baselines{n} = '91ft_PV_cable_cap_tests_Za_SC_Zb_SC_Zc_OC.lws';
data_center(n,:) = 35:49;
folders{n} = fullfile(script_dir, '../data/capacitor_data/asymmetric/');
files{n} = { ...       
        '91ft_PV_cable_cap_tests_Za_SC_Zb_1pf_Zc_OC.lws', ... 
        '91ft_PV_cable_cap_tests_Za_SC_Zb_2pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_5pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_8pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_15pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_20pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_25pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_47pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_100pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_330pf_Zc_OC.lws', ...
        '91ft_PV_cable_cap_tests_Za_SC_Zb_470pf_Zc_OC.lws', ...
         };
RE{n} = zeros(1,11);
CE{n} = [ 1 2 5 8 15 20 25 47 100 330 470]*1e-12;
LE{n} = zeros(1,11);
% ===========================     

n = n + 1;
% ===========================
% PARALLEL CAPACITORS (SIMULATED INTERFACE IMPEDANCES) 
CI{n} = 1e-9:-1e-12:1e-12;              % Parallel capacitance
RI{n} = 5e-3*ones(1,length(CI{n}));     % Parallel resistance
RLI{n} = 1e4*ones(1,length(CI{n}));     % Parallel resistance
LI{n} = 1e-9*ones(1,length(CI{n}));     % Parallel inductance (parasitic inductance)
Zf{n} = 1;                              % Symmetric = 1, Asymmetric = 1/2

% PARALLEL CAPACITORS (EXPERIMENTAL DATA) 
names{n} = ' Par. Cap.';
types{n} = 'parallel';
baselines{n} = 'parallel_c_Za_OC_Zc_OC.lws';
data_center(n,:) = 35:49;
folders{n} = fullfile(script_dir, '../data/Parallel measurement/parallel_c_measurments/');
files{n} = {...  
        'parallel_c_Za_1pf_Zc_OC.lws', ...
        'parallel_c_Za_2pf_Zc_OC.lws', ...
        'parallel_c_Za_5pf_Zc_OC.lws', ...
        'parallel_c_Za_8pf_Zc_OC.lws', ...
        'parallel_c_Za_15pf_Zc_OC.lws', ...
        'parallel_c_Za_20pf_Zc_OC.lws', ...
        'parallel_c_Za_25pf_Zc_OC.lws', ...
        'parallel_c_Za_47pf_Zc_OC.lws', ...
        'parallel_c_Za_100pf_Zc_OC.lws', ...
        'parallel_c_Za_330pf_Zc_OC.lws', ...
        'parallel_c_Za_470pf_Zc_OC.lws', ...    
        };
RE{n} = zeros(1,11);
CE{n} = [ 1 2 5 8 15 20 25 47 100 330 470]*1e-12;
LE{n} = zeros(1,11);     
% ===========================


%%
% ******************************************************************* %
% NO NEED TO CHANGE ANYTHING BEYOND THIS POINT
% ******************************************************************* %
     
% -------------------------------------------------
% DEFINE DERIVED PARAMETERS
% -------------------------------------------------
Q = PNL*Fsampl/(Fchip);                                 % Length of simulated signal
n = 1:Q;                                                % Sample axis
t = n/Fsampl;                                           % Time axis
r = floor(Q/2)+1;                                       % Center sample in frequency
f = ifftshift((n-r)/Q)*Fsampl;                          % Frequency axis
w = 2*pi*f;                                             % Angular frequency
N = length(files);                                      % Number of file groups

% -------------------------------------------------
% SPECIFY TRANSMISSION LINE IMPEDANCE AND PROPAGATION CONSTANT
% -------------------------------------------------
Z0 = sqrt(1j*w*L0 + R0)./sqrt(1j*w*C0 + G0);            % Characteristic impedance
gm = sqrt(1j*w*L0 + R0).*sqrt(1j*w*C0 + G0);            % Propagation factor

% -------------------------------------------------
% GENERATE SSTDR EXCITATION
% -------------------------------------------------
sqmod = [1 -1].';                                       % Square wave 
sq = kron(sqmod, ones(Fsampl/Fmodu/2,1));               % Square wave (expanded out based on modulation rate)
S  = fft(sq,Q).*conj(fft(sq,Q))/(Fsampl/(Fchip));       % Correlated square wave
S  = S./abs(S(sindx));                                  % Normalize at frequency of interest
s  = real(ifft(S));                                     % Correlated square wave

% INITIALIZE OUTPUTS
sim_maxamp          = cell(N,1);        % Simulated magnitude
sim_maxang          = cell(N,1);        % Simulated angle
exp_sub_maxamp      = cell(N,1);        % Experimental mangitude
exp_sub_maxang      = cell(N,1);        % Experimental angle
time_sig            = cell(N,1);        % Time-domain signal


%%
% -------------------------------------------------
% COMPUTE REFLECTION COEFFICIENTS FROM SIMULATIONS
% ------------------------------------------------
    
for nn = 1:N            % loop over simulated scenarios
    
    % LENGTH OF RESISTOR / CAPACITOR VALUES
    L = max([length(RI{nn}) length(LI{nn}) length(CI{nn})]);
    
     for ii = 1:L        % loop over resistor / capactior values
        
        % INTERFACE IMPEDANCE
        Zi = RI{nn}(ii) + 1j*w.*LI{nn}(ii) + 1./((1j*w.*CI{nn}(ii))) + 1/(RLI{nn}(ii));

        % SPECIFY REFLECTION COEFFICIENT
        if strcmpi(types{nn}, 'series')
            [X,x] = seriesreflsim(S,Z0,Z0,Zi*Zf{nn});   % X,x are frequency and time domain reflected signal respectively using our approach
        else
            [X,x] = parallelreflsim(S,Z0,Z0,Zi);        % X,x are frequency and time domain reflected signal respectively for the parallel case   
        end
        
        % SPECIFY REFLECTION COEFFICIENT [Our approach]
        sim_maxamp{nn}(ii) = abs(X(sindx));          % Reflection Coefficient
        sim_maxang{nn}(ii) = angle(X(sindx));        % Phase 
        
     end
     
 end 

%%
% -------------------------------------------------
% GET REFLECTION COEFFICIENTS FROM EXPERIMENTS
% -------------------------------------------------
for nn = 1:N        % loop over experiments

    % NUMBER OF FILES IN RUN
    L = length(files{nn});
    
    % GET BASELINE DATA (ONLY CONTAINS AN OPEN CIRCUIT)
    btmp = get_data([folders{nn} baselines{nn}]);                % Get data and then zero-pad
    b    = [btmp(data_center(nn,:)).'; zeros(Q-size(data_center,2),1)]; 
    b    = fshift(b, -(range(data_center(nn,:))/2));             % Shift to be centered at t=0
    
    for ii = 1:L      % loop over resistor / capacitor values
        
        % GET EXPERIMENTA DATA
        time_sig{nn}{ii} = get_data([folders{nn} files{nn}{ii}]);% Get data
        y    = [time_sig{nn}{ii}(data_center(nn,:)).'; zeros(Q-size(data_center,2),1)];
        y    = fshift(y, -(range(data_center(nn,:))/2));         % Shift to be center at t=0       
        z    = y-b;                                              % Baseline subtraction
        Z    = fft(z,Q);                                         % fourier transform
                
        % ESTIMATE REFLECTION COEFFICIENT
        [exp_sub_maxamp{nn}(ii)      ] = abs(Z(sindx));          % Reflection Coefficient
        [exp_sub_maxang{nn}(ii)      ] = angle(Z(sindx));        % Phase
        
    end
end



%%
% -------------------------------------------------
% PLOT RESULTS
% -------------------------------------------------

% INITIALIZE NORMALIZATION FACTOR THAT ALIGNS SIMULATION AND EXPERIMENTAL OPEN CIRCUIT 
nrm_fctr = zeros(N,1);

% PLOT MAGNIUTUDES
for nn = 1:N 

    % DETERMINE IF X-AXIS IS RESISTANCE, CAPACITANCE, OR INDUCTANCE
    [~, indx] = max([std(RI{nn}(:)) std(CI{nn}(:)) std(LI{nn}(:))]);
    if indx == 1, ZI = RI{nn}; ZE = RE{nn}; end
    if indx == 2, ZI = CI{nn}; ZE = CE{nn}; end
    if indx == 3, ZI = LI{nn}; ZE = LE{nn}; end
    
    % GENERATE NORMALIZATIONS 
    sim_maxamp_1 = interp1(RI{1}, sim_maxamp{1}, RE{1});        % Find simulation values at experimental impedances for first config
    sim_maxamp_2 = interp1(RI{3}, sim_maxamp{3}, RE{3});        % Find simulation values at experimental impedances for second config    
    sim_maxamp_i = interp1(ZI   , sim_maxamp{nn}, ZE);          % Find simulation values at experimental impedances for dataset nn
    if strcmpi(types{nn}, 'series')
         nrm_fctr(nn) =  max(sim_maxamp_1)./max(exp_sub_maxamp{1});  % Reduce all data sets by the max of the first data set (keep all data sets proportional)
    else
         nrm_fctr(nn) =  max(sim_maxamp_2)./max(exp_sub_maxamp{3});  % Reduce all data sets by the max of the first data set (keep all data sets proportional)
    end

    figure(1)
    set(gcf, 'units', 'inches', 'Position', [0 0 8 5])    
    subplot(2,N/2,nn)
    
    if Zf{nn} ~= 3
        plot(ZI, sim_maxamp{nn}, 'linewidth', 2);
        axis([min(ZI) max(ZI) -1 1])
        hold on;
        scatter(ZE, exp_sub_maxamp{nn}.*nrm_fctr(nn),'s','linewidth', 2);
        set(gca, 'Xscale', 'log')
        set(gca, 'fontsize', fontsize)
        hold off;
    else

        plot(ZI, sim_maxamp{nn}, 'linewidth', 2);
        axis([min(ZI) max(ZI) -1 1])
        hold on;
        scatter(ZE, exp_sub_maxamp{nn}.*nrm_fctr(nn),'s','linewidth', 2)
        set(gca, 'Xscale', 'log')
        set(gca, 'fontsize', fontsize)
        hold off;
    end
    axis([min(ZI) max(ZI) -1 1])
    if indx == 1, xlabel('Series Resistance [\Omega]'); end
    if indx == 2, xlabel('Series Capacitance [F]'); end
    if indx == 3, xlabel('Series Inductance [H]'); end
    if Zf{nn} == 3, xlabel('Parallel Capacitance [F]'); end
    if Zf{nn} == 4, xlabel('Parallel Resistors [\Omega]'); end
    
    ylabel('Ref. Coefficient')
    if indx == 1, legend('Our Theory', 'Experiment', 'Location', 'Southeast'); end
    if indx == 2, legend('Our Theory', 'Experiment',  'Location', 'Southeast'); end
    title(names{nn})
    
    % COMPUTE R^2 VALUES
    R2 = corr(sim_maxamp_i(:), exp_sub_maxamp{nn}(:).*nrm_fctr(nn)); 
    fprintf('%s: Magnitude R^2 value: %f \n', names{nn}, R2);
    
end
saveas(gcf, fullfile(script_dir, '../results/magnitudes.png'))

% PLOT PHASES
for nn = 1:N

    % DETERMINE IF X-AXIS IS RESISTANCE, CAPACITANCE, OR INDUCTANCE
    [~, indx] = max([std(RI{nn}(:)) std(CI{nn}(:)) std(LI{nn}(:))]);
    if indx == 1, ZI = RI{nn}; ZE = RE{nn}; end
    if indx == 2, ZI = CI{nn}; ZE = CE{nn}; end
    if indx == 3, ZI = LI{nn}; ZE = LE{nn}; end
    
    % PLOT RESULTS 
    figure(2)
    set(gcf, 'units', 'inches', 'Position', [0 0 8 5])    
    subplot(2,N/2,nn)

    %set(gcf, 'Units', 'Inches', 'Position', [0 4 4 1.7], 'color', [1 1 1]) %[0 2 4 2]
    plot(ZI, sim_maxang{nn}, 'linewidth', 2);
    axis([min(ZI) max(ZI) -pi pi])
    set(gca, 'fontsize', fontsize)
    hold on; 
    scatter(ZE, exp_sub_maxang{nn}, 's','linewidth', 2)   
    set(gca, 'Xscale', 'log')
    set(gca, 'fontsize', fontsize)
    hold off;
    axis([min(ZI) max(ZI) -pi pi])

    if indx == 1, xlabel('Series Resistance [\Omega]'); end
    if indx == 2, xlabel('Series Capacitance [F]'); end
    if indx == 3, xlabel('Series Inductance [H]'); end
    if Zf{nn} == 3, xlabel('Parallel Capacitance [F]'); end
    if Zf{nn} == 4, xlabel('Parallel Resistors [\Omega]'), legend('Our Theory', 'Experiment',  'Location', 'Southeast'); end
    
    ylabel('Ref. Phase Angle [rad]')
    title(names{nn})
    legend('Our Theory', 'Experiment')
    
    % COMPUTE R^2 VALUES
    sim_maxang_1 = interp1(ZI(:), sim_maxang{nn}(:), ZE(:));        % Find simulation values at experimental impedances for dataset 1
    R2 = corr(sim_maxang_1(:), exp_sub_maxang{nn}(:));
    fprintf('%s: Phase R^2 value: %f \n', names{nn}, R2);
end
saveas(gcf, fullfile(script_dir, '../results/phases.png'))

% PLOT TIME DOMAIN SIGNALS
for nn = 1:N
    % DETERMINE IF X-AXIS IS RESISTANCE, CAPACITANCE, OR INDUCTANCE
    [~, indx] = max([std(RI{nn}(:)) std(CI{nn}(:)) std(LI{nn}(:))]);
    if indx == 1, ZI = RI{nn}; ZE = RE{nn}; end
    if indx == 2, ZI = CI{nn}; ZE = CE{nn}; end
    if indx == 3, ZI = LI{nn}; ZE = LE{nn}; end
    
    % PLOT TIME-DOMAIN RESULTS
    L = length(files{nn});
    figure(3)
    set(gcf, 'units', 'inches', 'Position', [0 0 8 5])    
    subplot(2,N/2,nn)
    plottype = {'-','-','--','--',':b',':','-.r','-.','-.k','-.k*','--m','--ko'};
    for ii = 1:2:L-1
        [~, zero_samp] = max(time_sig{nn}{ii}(1:data_center(nn,1)));
        tdist = 0.5*((-zero_samp/Fsampl):(1/Fsampl):((length(time_sig{nn}{ii})-zero_samp-1)/Fsampl))*vop*2.99792*1e8;
        plot(tdist, time_sig{nn}{ii}*nrm_fctr(nn), plottype{ii}); 
        hold on;
        set(gca, 'fontsize', fontsize)
    end
    hold off;

    if indx   == 1, legend("10\Omega","30\Omega","130\Omega","1000\Omega", 'Location', 'Northeast'); end
    if indx   == 2, legend("1pF","5pF","15pF","25pF","100pF"); end
    if indx   == 3, xlabel('Series Inductance [H]'); end
    if Zf{nn} == 3, legend("1pF","5pF","15pF","25pF","100pF"); end
    if Zf{nn} == 4, legend("10\Omega","30\Omega","130\Omega","1000\Omega", 'Location', 'Northeast'); end
    
    axis([tdist(data_center(nn,1)) tdist(data_center(nn,end)) -0.5 1.5])
    xlabel('Distance (meters)','fontsize', fontsize)
    ylabel('Ref. Coeff.')
    title(names{nn})

end
saveas(gcf, fullfile(script_dir, '../results/time.png'))

