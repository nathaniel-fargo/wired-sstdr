%% New S/SSTDR simulation using predetermined incident signals,
% reflections from multiple genterators,
% interference singals,
% and continuous FFT
% NOTE: This code uses signal files created by Mouad and stored in "CSV_Sig_path"

close all;
CSV_Sig_path = 'H:\My Drive\U of U\Post PhD Research\PN Code Interference Tests';
addpath(CSV_Sig_path)

nextFig = 1; % Initialize figure numbers
saveFigs = 1; % flag to save figures
exportDataCSV = 0; % flag for exporting simulation data

%% Load organize the CSV files into a matrices
% Get directory information from CSV data folder
Sigs = dir(fullfile(CSV_Sig_path, '*.csv'));

% Struct search fields for organizing data
mseq_interf_stringSearch = 'Interf_mseq*';
gold_interf_stringSearch = 'Interf_Gold*';
z_interf_stringSearch = 'Interf_Z*';
mseq_stringSearch = 'mseq*';
gold_stringSearch = 'Gold*';
z_stringSearch = 'Z*';

% Extracting data into matrices
for n = 1:length(Sigs)
    if contains(Sigs(n).name, mseq_interf_stringSearch(1:end-1))
        mseq_interf_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        mseq_inferf_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    elseif contains(Sigs(n).name, gold_interf_stringSearch(1:end-1))
        gold_interf_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        gold_inferf_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    elseif contains(Sigs(n).name, z_interf_stringSearch(1:end-1))
        z_interf_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        z_inferf_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    elseif contains(Sigs(n).name, mseq_stringSearch(1:end-1))
        mseq_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        mseq_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    elseif contains(Sigs(n).name, gold_stringSearch(1:end-1))
        gold_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        gold_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    elseif contains(Sigs(n).name, z_stringSearch(1:end-1))
        z_sig = readmatrix(strcat(Sigs(n).folder,'\',Sigs(n).name)); % Pulls a two column matrix. Col 1 is time, Col 2 is amplitude
        z_name = replace(Sigs(n).name,'.csv',''); % Pulls the name of the signal for use later
    end
end

%% Extract signal information from signal matrices and setup frequency 
delta_T = (mseq_sig(2,1)-mseq_sig(1,1))/2; % All signals have same delta_T,divide by 2 because reflection
sampFreq = (2*delta_T)^(-1); % All signals have same sample rate

% Setup conventional frequency axis
% mseq and gold codes
nfreq_MG = length(mseq_sig(:,1)); % Number of frequency domain points
r_s_MG = floor(nfreq_MG/2);
f_neg_pos_MG = (((((1:nfreq_MG)-r_s_MG)/(nfreq_MG))).*sampFreq)'; % Frequency Axis

% Z codes
nfreq_Z = length(z_sig(:,1)); % Number of frequency domain points
r_s_Z = floor(nfreq_Z/2);
f_neg_pos_Z = (((((1:nfreq_Z)-r_s_Z)/(nfreq_Z))).*sampFreq)'; % Frequency Axis

%% Reflection simulation information
% Get zero distance reflection for each signal
mseq_Idealresponse = get_incident_TD_response(mseq_sig,0,0); % TD Response = funtion(signal, distance, propagation constant)
gold_Idealresponse = get_incident_TD_response(gold_sig,0,0);
z_Idealresponse = get_incident_TD_response(z_sig,0,0);

% System reflection coefficients
    % The system consists of two generators and one oscilloscope. Gen1 is
    % connected directly to the oscilloscope via a T-connector (T1). After 
    % T1, there is 100ft of RG58 until a second T-connector (T2). From T2, 
    % 50ft of RG58 connect to Gen2 (O,S,M simulatable). Also from T2, 100ft
    % of RG58 is terminated with an open circuit (OC). See diagram below:

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gen1-T1-----A=100ft or 30.48m-----T2-----B=102.5ft or 31.24m-----(O,S,M)
    %      |                            |
    %    Scope                          L---C=45ft or 13.716m---Gen2 (O,S,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % NOTE: we are using the reflection from T1 as the zero distance point 
    % for the system

    % Expected reflections @ 0ft, 100ft(30.48m), 145ft(44.2m), [190ft 200ft 202.5ft] ([57.91m  60.96m 61.75]) 

% Cable length (in feet) values for the diagram above
cableA = 100; 
cableB = 102.5;
cableC = 45;

% Cable terminations for reflection coefficients at the appropriate locations
cableA_gamma = -1/3; % Z0 = (Z0||Z0 - Z0)/(Z0||Z0 + Z0)
cableB_gamma = 1; % 1 = OC, -1 = SC, 0=Matched
cableC_gamma = 1; % 1 = OC, -1 = SC, 0=Matched

% The following bounce information was calculated using the image below
    % figure(nextFig)
    % imshow('BounceDiagram_Gen1_ConfigA_PNG.png')
    % title('Bounce Diagram from Generator 1 in Configuration A')
    % nextFig = nextFig+1;

% Gamma from bounces for signal from Gen1
GT1=cableA_gamma; 
GT2=GT1;
GB=cableB_gamma;
GC=cableC_gamma;
V1p=1;
V1m=V1p*GT2;
V2p = V1m*GT1;
V2m=V2p*GT2;
V3p=V2m*GT1;
V4p=V1p*(1+GT2);
V4m=V4p*GB;
V5p=V4m*GT1;
V6m=V4m*(1+GT1);
V6p=V6m*GT1;
V6pp = V6m+GT1;
V7p=V4p;
V7m=V7p*GC;
V8m=V7m*(1+GT2);
V9p=V8m*GT1;
V10p = V7m*GT2;
V10m = V10p+GC;
V11m = V10m*(1+GT2);
V11p = V11m*GT1;

gamma_T1 = 1;
gamma_T2 = GT2; %V1m+V2p;
gamma_AplusC = V8m*(1+GT1);
gamma_V2m_V3p = V2m+V3p;
gamma_V6m_V6p = V6m+V6p;
gamma_V11m_V11p = V11m+V11p;
gamma_AplusB = V2m+V6m+V11m+V3p+V6p+V11p;


% The following bounce information was calculated using the image below
    % figure(nextFig)
    % imshow('BounceDiagram_Gen2_ConfigA_PNG.png')
    % title('Bounce Diagram from Generator 2 in Configuration A')
    % nextFig = nextFig+1;

% Gamma from bounces for interference signals from Gen2 - GC is always zero
% due to Gen2 being matched
iV1p = 1;
iV2p = iV1p*(1+GT2);
iV2m = iV2p*GB;
iV3p = iV2m*GT2;
iV3m = iV3p*GB;
iV4p = iV2p;
iV4m = iV4p*GT1;
iV5p = iV4m*GT2;
iV5m = iV5p*GT1;
iV6p = iV4m*GT2;
iV6m = iV6p*GB;
iV7p = iV3p;
iV7m = iV7p*GT1;

igamma_T2 = GT2;
igamma_T1 = GT1;
igamma_iV4p_iV4m = iV4p + iV4m;
igamma_iV5p_iV5m = iV5p + iV5m;
igamma_iV7p_iV7m = iV7p + iV7m;

% % Propagation speeds - uses ideal lossless characteristics
c = 3e8; % Speed of light m/s
prop_speed_m_s = 2/3*c; % Propagation speed m/s
prop_speed_ft_s = prop_speed_m_s*3.28084; % Propagation speed ft/s

% Propagation speed calculation - uses RG58 chatacteristics
[Z0_MG, VOP_MG, prop_constant_MG] = get_RG58TLCharParam(f_neg_pos_MG); % [Z0, VOP, propagation constant] = func(signal_frequencies)
[Z0_Z, VOP_Z, prop_constant_Z] = get_RG58TLCharParam(f_neg_pos_Z); % [Z0, VOP, propagation constant] = func(signal_frequencies)

% % Time/distance delays to build complete response
%     % We will use time delay from Gen1 to build each reflection at expected
%     % point in time/distance, then sum all signals together.
% T_delay_T2 = cableA/prop_speed_ft_s;
% T_delay_AplusC = (cableA + cableC)/prop_speed_ft_s;
% T_delay_AplusB = (cableA + cableB)/prop_speed_ft_s;

% % Number of samples to delay each response. NOTE: introduces rounding error
% delta_T_samples_T2 = floor(T_delay_T2/delta_T);
% delta_T_samples_AplusC = floor(T_delay_AplusC/delta_T);
% delta_T_samples_AplusB = floor(T_delay_AplusB/delta_T);

%% Get delayed and attenuated incident signals using distance and prop_constant - no interference
% MSEQ Codes
mseq_delayAtt_T1 = gamma_T1*get_incident_TD_response(mseq_sig, 0, prop_constant_MG); % TD Response = funtion(signal, distance, propagation constant)
mseq_delayAtt_T2 = gamma_T2*get_incident_TD_response(mseq_sig, cableA, prop_constant_MG);
mseq_delayAtt_AplusC = gamma_AplusC*get_incident_TD_response(mseq_sig, cableA+cableC, prop_constant_MG);
mseq_delayAtt_V2m = gamma_V2m_V3p*get_incident_TD_response(mseq_sig, 2*cableA, prop_constant_MG);
mseq_delayAtt_V6m = gamma_V6m_V6p*get_incident_TD_response(mseq_sig, cableA+cableB, prop_constant_MG);
mseq_delayAtt_V11m = gamma_V11m_V11p*get_incident_TD_response(mseq_sig, cableA+2*cableC, prop_constant_MG);

mseq_delayAtt_final = (mseq_delayAtt_T1 + mseq_delayAtt_T2 + mseq_delayAtt_V2m + mseq_delayAtt_V6m + mseq_delayAtt_V11m + mseq_delayAtt_AplusC)/max(mseq_delayAtt_T1);

% GOLD Codes
gold_delayAtt_T1 = gamma_T1*get_incident_TD_response(gold_sig, 0, prop_constant_MG); % TD Response = funtion(signal, distance, propagation constant)
gold_delayAtt_T2 = gamma_T2*get_incident_TD_response(gold_sig, cableA, prop_constant_MG);
gold_delayAtt_AplusC = gamma_AplusC*get_incident_TD_response(gold_sig, cableA+cableC, prop_constant_MG);
gold_delayAtt_V2m = gamma_V2m_V3p*get_incident_TD_response(gold_sig, 2*cableA, prop_constant_MG);
gold_delayAtt_V6m = gamma_V6m_V6p*get_incident_TD_response(gold_sig, cableA+cableB, prop_constant_MG);
gold_delayAtt_V11m = gamma_V11m_V11p*get_incident_TD_response(gold_sig, cableA+2*cableC, prop_constant_MG);

gold_delayAtt_final = (gold_delayAtt_T1 + gold_delayAtt_T2 + gold_delayAtt_V2m + gold_delayAtt_V6m + gold_delayAtt_V11m + gold_delayAtt_AplusC)/max(gold_delayAtt_T1);

% ZtoZ codes
z_delayAtt_T1 = gamma_T1*get_incident_TD_response(z_sig, 0, prop_constant_Z); % TD Response = funtion(signal, distance, propagation constant)
z_delayAtt_T2 = gamma_T2*get_incident_TD_response(z_sig, cableA, prop_constant_Z);
z_delayAtt_AplusC = gamma_AplusC*get_incident_TD_response(z_sig, cableA+cableC, prop_constant_Z);
z_delayAtt_V2m = gamma_V2m_V3p*get_incident_TD_response(z_sig, 2*cableA, prop_constant_Z);
z_delayAtt_V6m = gamma_V6m_V6p*get_incident_TD_response(z_sig, cableA+cableB, prop_constant_Z);
z_delayAtt_V11m = gamma_V11m_V11p*get_incident_TD_response(z_sig, cableA+2*cableC, prop_constant_Z);

z_delayAtt_final = (z_delayAtt_T1 + z_delayAtt_T2 + z_delayAtt_V2m + z_delayAtt_V6m + z_delayAtt_V11m + z_delayAtt_AplusC)/max(z_delayAtt_T1);

%% Get delayed and attenuated incident signals using distance and prop_constant - with interference
% MSEQ Codes
interf_mseq_delayAtt_iV4p_iV4m = igamma_iV4p_iV4m*get_interfered_TD_response(mseq_sig, mseq_interf_sig, cableC + cableA, prop_constant_MG); % interf_TD = get_interfered_TD_response(incSig, interfSig, z_ft, prop)
interf_mseq_delayAtt_iV5p_iV5m = igamma_iV5p_iV5m*get_interfered_TD_response(mseq_sig, mseq_interf_sig, cableC + 3*cableA, prop_constant_MG);
interf_mseq_delayAtt_iV7p_iV7m = igamma_iV7p_iV7m*get_interfered_TD_response(mseq_sig, mseq_interf_sig, cableC + cableA + 2*cableB, prop_constant_MG);

interf_mseq_delayAtt_temp = (interf_mseq_delayAtt_iV4p_iV4m + interf_mseq_delayAtt_iV5p_iV5m + interf_mseq_delayAtt_iV7p_iV7m)/max(mseq_delayAtt_T1);
interf_mseq_delayAtt_final = mseq_delayAtt_final + interf_mseq_delayAtt_temp;

% GOLD Codes
interf_gold_delayAtt_iV4p_iV4m = igamma_iV4p_iV4m*get_interfered_TD_response(gold_sig, gold_interf_sig, cableC + cableA, prop_constant_MG); % interf_TD = get_interfered_TD_response(incSig, interfSig, z_ft, prop)
interf_gold_delayAtt_iV5p_iV5m = igamma_iV5p_iV5m*get_interfered_TD_response(gold_sig, gold_interf_sig, cableC + 3*cableA, prop_constant_MG);
interf_gold_delayAtt_iV7p_iV7m = igamma_iV7p_iV7m*get_interfered_TD_response(gold_sig, gold_interf_sig, cableC + cableA + 2*cableB, prop_constant_MG);

interf_gold_delayAtt_temp = (interf_gold_delayAtt_iV4p_iV4m + interf_gold_delayAtt_iV5p_iV5m + interf_gold_delayAtt_iV7p_iV7m)/max(gold_delayAtt_T1);
interf_gold_delayAtt_final = gold_delayAtt_final + interf_gold_delayAtt_temp;

% ZtoZ Codes
interf_z_delayAtt_iV4p_iV4m = igamma_iV4p_iV4m*get_interfered_TD_response(z_sig, z_interf_sig, cableC + cableA, prop_constant_Z); % interf_TD = get_interfered_TD_response(incSig, interfSig, z_ft, prop)
interf_z_delayAtt_iV5p_iV5m = igamma_iV5p_iV5m*get_interfered_TD_response(z_sig, z_interf_sig, cableC + 3*cableA, prop_constant_Z);
interf_z_delayAtt_iV7p_iV7m = igamma_iV7p_iV7m*get_interfered_TD_response(z_sig, z_interf_sig, cableC + cableA + 2*cableB, prop_constant_Z);


interf_z_delayAtt_temp = (interf_z_delayAtt_iV4p_iV4m + interf_z_delayAtt_iV5p_iV5m + interf_z_delayAtt_iV7p_iV7m)/max(z_delayAtt_T1);
interf_z_delayAtt_final = z_delayAtt_final + interf_z_delayAtt_temp;

% Build Ideal response vectors - circshift applies time delay for propagation
% mseq_respVect_T1 = gamma_T1*mseq_Idealresponse;
% mseq_respVect_T2 = gamma_T2*circshift(mseq_Idealresponse,delta_T_samples_T2);
% mseq_respVect_AplusC = gamma_AplusC*circshift(mseq_Idealresponse,delta_T_samples_AplusC);
% 
% mseq_respVect_AplusB = gamma_AplusB*circshift(mseq_Idealresponse,delta_T_samples_AplusB);
% 
% % Add all response vectors and normalize to be plotted
% mseq_response_final = (mseq_respVect_T1 + mseq_respVect_T2 + mseq_respVect_AplusC + mseq_respVect_AplusB)/max(mseq_respVect_T1);

% Determine the refelction from T1 to set as zero distance and determine
% x-axis MSEQ and GOLD
[~,peakInd_MG] = max(mseq_delayAtt_final); % Find zero index
neg_time_MG = -flip(delta_T*linspace(1,peakInd_MG-1,peakInd_MG-1)'); % Set negative time values
pos_time_MG = delta_T*linspace(0,length(mseq_delayAtt_final)+1 - peakInd_MG,length(mseq_delayAtt_final)+1 - peakInd_MG)'; % Set positive time values
x_axis_MG = prop_speed_m_s*[neg_time_MG; pos_time_MG]; % Construct distance x-axis

% x-axis ZtoZ
[~,peakInd_Z] = max(z_delayAtt_final); % Find zero index
neg_time_Z = -flip(delta_T*linspace(1,peakInd_Z-1,peakInd_Z-1)'); % Set negative time values
pos_time_Z = delta_T*linspace(0,length(z_delayAtt_final)+1 - peakInd_Z,length(z_delayAtt_final)+1 - peakInd_Z)'; % Set positive time values
x_axis_Z = prop_speed_m_s*[neg_time_Z; pos_time_Z]; % Construct distance x-axis

%% Plotting in Time Domain
% figure(nextFig);
% plot(x_axis,mseq_response_final)
% xlim([-20 100])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title("Ideal - No attenuation and Time delayed")
% nextFig=nextFig+1;
       
% figure(nextFig)
% plot(x_axis_MG,real(mseq_delayAtt_final))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal: %s',mseq_name))
% nextFig = nextFig+1;
% 
% figure(nextFig)
% plot(x_axis_MG,real(mseq_delayAtt_final + interf_mseq_delayAtt_temp))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal w/interf: %s',mseq_name))
% nextFig = nextFig+1;
% 
% figure(nextFig)
% plot(x_axis_MG,real(gold_delayAtt_final))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal: %s',gold_name))
% nextFig = nextFig+1;
% 
% figure(nextFig)
% plot(x_axis_MG,real(gold_delayAtt_final + interf_gold_delayAtt_temp))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal w/interf: %s',gold_name))
% nextFig = nextFig+1;
% 
% figure(nextFig)
% plot(x_axis_Z,real(z_delayAtt_final))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal: %s',z_name))
% nextFig = nextFig+1;
% 
% figure(nextFig)
% plot(x_axis_Z,real(z_delayAtt_final + interf_z_delayAtt_temp))
% xlim([-20 75])
% xlabel("Distance [m]")
% ylabel("Correlation Magnitude")
% title(sprintf('Test Signal w/interf: %s',z_name))
% nextFig = nextFig+1;

mseq_fig = figure(nextFig);
plot(x_axis_MG,real(mseq_delayAtt_final), 'linewidth', 2)
hold on;
plot(x_axis_MG,real(interf_mseq_delayAtt_final), 'linewidth', 2)
xlim([-20 75])
xlabel("Distance [m]")
ylabel("Correlation Magnitude")
title(sprintf('Test Signal: %s',mseq_name))
legend("No Interference", "With Interference")
nextFig = nextFig+1;

gold_fig = figure(nextFig);
plot(x_axis_MG,real(gold_delayAtt_final), 'linewidth', 2)
hold on;
plot(x_axis_MG,real(interf_gold_delayAtt_final), 'linewidth', 2)
xlim([-20 75])
xlabel("Distance [m]")
ylabel("Correlation Magnitude")
title(sprintf('Test Signal: %s',gold_name))
legend("No Interference", "With Interference")
nextFig = nextFig+1;

z_fig = figure(nextFig);
plot(x_axis_Z,real(z_delayAtt_final), 'linewidth', 2)
hold on;
plot(x_axis_Z,real(interf_z_delayAtt_final), 'linewidth', 2)
xlim([-20 75])
xlabel("Distance [m]")
ylabel("Correlation Magnitude")
title(sprintf('Test Signal: %s',z_name))
legend("No Interference", "With Interference")
nextFig = nextFig+1;

% Exporting simulation data as csv files
if exportDataCSV == 1
    CSV_Path = 'H:\My Drive\U of U\Post PhD Research\CSV Files for Mouad';
    mseq_sim_CSV = [x_axis_MG,mseq_delayAtt_final,interf_mseq_delayAtt_final];
    gold_sim_CSV = [x_axis_MG,gold_delayAtt_final,interf_gold_delayAtt_final];
    z_sim_CSV = [x_axis_Z,z_delayAtt_final,interf_z_delayAtt_final];

    mseq_filePath = fullfile(CSV_Path, sprintf('%s_sim_CSV.csv',mseq_name));
    gold_filePath = fullfile(CSV_Path, sprintf('%s_sim_CSV.csv',gold_name));
    z_filePath = fullfile(CSV_Path, sprintf('%s_sim_CSV.csv',z_name));
    
    writematrix(mseq_sim_CSV, mseq_filePath);
    writematrix(gold_sim_CSV, gold_filePath);
    writematrix(z_sim_CSV, z_filePath);
end

% Saving figures as PNG
if saveFigs == 1
    fig_Path = 'H:\My Drive\U of U\MATLAB\Post PhD Research\Interference Work with Mouad\Output Figs';
    
    mseq_figPath = fullfile(fig_Path, sprintf('TestSignal_%s.png',mseq_name));
    gold_figPath = fullfile(fig_Path, sprintf('TestSignal_%s.png',gold_name));
    z_figPath = fullfile(fig_Path, sprintf('TestSignal_%s.png',z_name));

    exportgraphics(mseq_fig, mseq_figPath);
    exportgraphics(gold_fig, gold_figPath);
    exportgraphics(z_fig, z_figPath);
end
%% Support functions
% Periodic correlation to get time domain response - includes attenuation
% and time delay based on distance in feet
function response = get_incident_TD_response(incSig, z_ft, prop)
    z_meters = z_ft/3.28084;
    inc_fft = fft(incSig(:,2));
    ref_fft = fft(incSig(:,2));
    ref_delay_att = exp(-2*z_meters*prop).*ref_fft;
    inc_ref_corr = inc_fft.*conj(ref_delay_att);
    % inc_ref_corr = inc_fft.*conj(ref_delay_att);
    reflected_TD_temp = ifft(inc_ref_corr);
    response = fftshift(reflected_TD_temp);
end

% This function takes the incident signal, interf signal, distance of
% propgation, propagation constant and determines the TD response
function interf_TD = get_interfered_TD_response(incSig, interfSig, z_ft, prop)
    z_meters = z_ft/3.28084;

    inc_fft = fft(incSig(:,2));
    interf_fft = fft(interfSig(:,2));

    interf_ref_delay_att = exp(-z_meters*prop).*interf_fft;

    inc_interf_corr = inc_fft.*conj(interf_ref_delay_att);
    interf_TD_temp = ifft(inc_interf_corr);
    interf_TD = fftshift(interf_TD_temp);
end








