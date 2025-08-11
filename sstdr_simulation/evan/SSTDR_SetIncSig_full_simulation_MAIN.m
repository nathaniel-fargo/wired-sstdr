%% Simulation code for evaluating S/SSTDR for impedance extraction under various scenarios.
% By Evan Benoit

% Needed Functions:
    % 

close all

fm = 24e6;                  % Modulation frequency
SRtester = 1;                % For scaling the Sampling Rate
fs_ss = SRtester*4*fm;       % Sampling frequency for SSTDR
fs_s = SRtester*2*fm;        % Sampling frequency for STDR

% PN code information
PnL = 15;                               % PN code bit length
mod_method = 1;                         % Modulation method 1=square, 2=sine, 3=sawtooth
nPN_onOFF = 1;                          % Toggle for N-port SS-VNA simulation 
nPN = 4;                                % Number of PN codes to create in simulation (code is set to work with 2^N number of PN codes)
PnData = [PnL, mod_method, nPN_onOFF, nPN];
% grab_nPN_mean_mat = 1;      % Flag for grabbing nPN mean error data for plotting elsewhere

% Attenuation information
att = "off";         % Toggle to include/not include TL attenuation
lengthTL = 0;     % Length in METERS of simulated TL for attenuation
TL_type = "RG58";           % Set the type of TL you are using. Options: "RG58", "PV Cable"

% Frequency domain points
nfreq_ss = 50e3;%fs_ss/fm*500;
nfreq_s = nfreq_ss/2;%fs_s/fm*500;
% nfreq_ss = fs_ss/fm*500;
% nfreq_s = fs_s/fm*500;

% Figure number tracking and figure variables
nextFig = 1;
omega = "\Omega";
gamma = "\Gamma";
% nextFig = lastFig+1;

% Compute frequency axis
r_ss = floor(nfreq_ss/2);                        % Center of frequency axis
r_s = floor(nfreq_s/2);
f_neg_pos_ss = (((((1:nfreq_ss)-r_ss)/(nfreq_ss))).*fs_ss)';  % Conventional freqeuncy axis
f_neg_pos_s = (((((1:nfreq_s)-r_s)/(nfreq_s))).*fs_s)';  % Conventional freqeuncy axis
 
fpos_start_ss = find(f_neg_pos_ss==0)+1;
fpos_start_s = find(f_neg_pos_s==0)+1;
fpos_end_ss = find(f_neg_pos_ss==2*fm);
fpos_end_s = find(f_neg_pos_s==fm);
f_pos_ss = f_neg_pos_ss(fpos_start_ss:fpos_end_ss);
f_pos_s = f_neg_pos_s(fpos_start_s:fpos_end_s);
X_FDaxisEnd = ceil(2*fm/10e6)*10;
% Compute time axis
ts_ss = 1/fs_ss;
ts_s = 1/fs_s;
t_axis_ss = ts_ss*(1:nfreq_ss);
t_axis_s = ts_s*(1:nfreq_s);

% TL parameters
% [Z0_ext,~,~] = get_TLCharParam(TL_type,fshift(f_neg_pos,length(f_neg_pos)/2));
[Z0_ext_ss,~,~] = get_TLCharParam(TL_type,f_neg_pos_ss);
[Z0_ext_s,~,~] = get_TLCharParam(TL_type,f_neg_pos_s);
% Creates the incident S/SSTDR signals
if PnData(3) == 0 % Forces the number of PN codes to be 1 when nPN_onOFF is set to 0
    nPN = 1;
    PnData = [PnL, mod_method, nPN_onOFF, nPN];
end
[~,inc_SSTDR] = get_IncidentSigs(fm,fs_ss,PnData);
[inc_STDR,~] = get_IncidentSigs(fm,fs_s,PnData);

% Insert new signals here. 
% Need new sample rate, 

%% Combines all PN modulated codes into a single signal
if nPN_onOFF == 1 && nPN > 1
    inc_SSTDR_copies = inc_SSTDR;
    inc_SSTDR = sum(inc_SSTDR(:,:)')';
    inc_STDR_copies = inc_STDR;
    inc_STDR = sum(inc_STDR(:,:)')';
elseif nPN_onOFF == 1 && nPN == 1
    inc_SSTDR_copies = inc_SSTDR;
    inc_STDR_copies = inc_STDR;
end

% TD windowing - Uses sample rate - TD_Pulse_Width = 2*(fs/fm)+1
sample_rate_ss = fs_ss/fm; % Assuming standard fs = 4*fm
sample_rate_s = fs_s/fm;
% Windowing assumes the correlated pulse peak located @ length(corr(sig))/2
window_left_ss = length(inc_SSTDR) - sample_rate_ss - 5*SRtester; % sample rate -1 is the normal window: increase the constant to include more of the time domain response
window_right_ss = length(inc_SSTDR) + sample_rate_ss + 5*SRtester;
window_left_s = length(inc_STDR) - sample_rate_s - 5*SRtester; % sample rate -1 is the normal window: increase the constant to include more of the time domain response
window_right_s = length(inc_STDR) + sample_rate_s + 5*SRtester;
stdr_window_comp = 0;

% Compute frequency axis
nfreq2_ss = length(inc_SSTDR);
nfreq2_s = length(inc_STDR);
r2_ss = floor(nfreq2_ss/2);                        % Center of frequency axis
r2_s = floor(nfreq2_s/2);                        % Center of frequency axis
f_neg_pos2_ss = (((((1:nfreq2_ss)-r2_ss)/(nfreq2_ss))).*fs_ss)';  % Conventional freqeuncy axis
f_neg_pos2_s = (((((1:nfreq2_s)-r2_s)/(nfreq2_s))).*fs_s)';  % Conventional freqeuncy axis

fpos2_start_ss = find(f_neg_pos2_ss==0)+1;
fpos2_start_s = find(f_neg_pos2_s==0)+1;
fpos2_end_ss = find(f_neg_pos2_ss==2*fm);
fpos2_end_s = find(f_neg_pos2_s==fm);
f_pos2_ss = f_neg_pos2_ss(fpos2_start_ss:fpos2_end_ss); 
f_pos2_s = f_neg_pos2_s(fpos2_start_s:fpos2_end_s); 

%% Get Load Data
res_Load = [75     ];
ind_Load = [];
cap_Load = [100e-12  ];


% [cap_Load,~] = get_FreqDepLoad(fm,res_Load);

OCSC = 1;   % Flag for adding OC and SC to the list of loads
[load_table,load_mat,load_list] = get_loadList(res_Load,ind_Load,cap_Load,OCSC);

%% Get reflected signals using Load Data

% TL parameters
% [Z0_ref,~,~] = get_TLCharParam(TL_type,fshift(f_neg_pos2,length(f_neg_pos2)/2));
[Z0_ref_ss,~,prop_constant_ss] = get_TLCharParam(TL_type,f_neg_pos2_ss);
[Z0_ref_s,~,prop_constant_s] = get_TLCharParam(TL_type,f_neg_pos2_s);

[~,ref_SSTDR,ZL_mat_ss,~] = get_ReflectedSigs(inc_SSTDR,inc_SSTDR,load_mat,Z0_ref_ss,f_neg_pos2_ss,nfreq2_ss);
[ref_STDR,~,ZL_mat_s,~] = get_ReflectedSigs(inc_STDR,inc_STDR,load_mat,Z0_ref_s,f_neg_pos2_s,nfreq2_s);

% Recalculate ZL Matrix to use with NFREQ FFT vise NFREQ2 FFT
[~,~,ZL_mat_short_ss,Gamma_mat_ss] = get_ReflectedSigs(inc_SSTDR,inc_SSTDR,load_mat,Z0_ext_ss,f_neg_pos_ss,nfreq_ss);
[~,~,ZL_mat_short_s,Gamma_mat_s] = get_ReflectedSigs(inc_STDR,inc_STDR,load_mat,Z0_ext_s,f_neg_pos_s,nfreq_s);

%% Get characteristic TD and FD data
% disp("Getting correlated TD and FD signals...")
dB_vect = [0]; % dB = 0 does not introduce noise - 50 dB is essentially not adding noise

if nPN_onOFF == 0
    [ext_SSTDR_FD,~,ref_SSTDR_TD,ref_SSTDR_FD,~,~]...
        = Analysis_AddSomeNoise(inc_SSTDR,ref_SSTDR,inc_SSTDR,ref_SSTDR,dB_vect,window_left_ss,window_right_ss,nfreq_ss,Z0_ext_ss);
    [~,ext_STDR_FD,~,~,ref_STDR_TD,ref_STDR_FD]...
        = Analysis_AddSomeNoise(inc_STDR,ref_STDR,inc_STDR,ref_STDR,dB_vect,window_left_s-stdr_window_comp,window_right_s+stdr_window_comp,nfreq_s,Z0_ext_s);
    [ext_SSTDR_FD_2,~,ref_SSTDR_TD_2,ref_SSTDR_FD_2,~,~]...
        = Analysis_AddSomeNoise(inc_SSTDR,ref_SSTDR,inc_SSTDR,ref_SSTDR,dB_vect,window_left_ss,window_right_ss,nfreq2_ss,Z0_ref_ss);
    [~,ext_STDR_FD_2,~,~,ref_STDR_TD_2,ref_STDR_FD_2]...
        = Analysis_AddSomeNoise(inc_STDR,ref_STDR,inc_STDR,ref_STDR,dB_vect,window_left_s-stdr_window_comp,window_right_s+stdr_window_comp,nfreq2_s,Z0_ref_s);
elseif nPN_onOFF == 1
    ext_SSTDR_FD_nPN_mat = [];
    ext_STDR_FD_nPN_mat = [];
    ref_SSTDR_TD_nPN_mat = [];
    ref_STDR_TD_nPN_mat = [];
    ref_SSTDR_FD_nPN_mat = [];
    ref_STDR_FD_nPN_mat = [];
    for pn = 1:nPN % Parses through each copy of the S/SSTDR incident signals and performs load extraction using each one
        % Output format of Matrices - MATRIX(nPN,:,nLoad)
        currentINC_SSTDR = inc_SSTDR_copies(:,pn);
        currentINC_STDR = inc_STDR_copies(:,pn);
        
        [ext_SSTDR_FD,~,ref_SSTDR_TD,ref_SSTDR_FD,~,~]...
            = Analysis_AddSomeNoise(currentINC_SSTDR,ref_SSTDR,currentINC_SSTDR,ref_SSTDR,0,window_left_ss,window_right_ss,nfreq_ss,Z0_ext_ss);
        [~,ext_STDR_FD,~,~,ref_STDR_TD,ref_STDR_FD]...
            = Analysis_AddSomeNoise(currentINC_STDR,ref_STDR,currentINC_STDR,ref_STDR,0,window_left_s-stdr_window_comp,window_right_s+stdr_window_comp,nfreq_s,Z0_ext_s);

        ext_SSTDR_FD_nPN_mat(pn,:,:) = squeeze(ext_SSTDR_FD);
        ext_STDR_FD_nPN_mat(pn,:,:) = squeeze(ext_STDR_FD);
        ref_SSTDR_TD_nPN_mat(pn,:,:) = squeeze(ref_SSTDR_TD);
        ref_STDR_TD_nPN_mat(pn,:,:) = squeeze(ref_STDR_TD);
        ref_SSTDR_FD_nPN_mat(pn,:,:) = squeeze(ref_SSTDR_FD);
        ref_STDR_FD_nPN_mat(pn,:,:) = squeeze(ref_STDR_FD);
        
        if pn == nPN
            [~,ref_SSTDR,ZL_mat_ss,~] = get_ReflectedSigs(currentINC_SSTDR,currentINC_SSTDR,load_mat,Z0_ref_ss,f_neg_pos2_ss,nfreq2_ss);
            [ref_STDR,~,ZL_mat_s,~] = get_ReflectedSigs(currentINC_STDR,currentINC_STDR,load_mat,Z0_ref_s,f_neg_pos2_s,nfreq2_s);
            
            [ext_SSTDR_FD,~,ref_SSTDR_TD,ref_SSTDR_FD,~,~]...
                = Analysis_AddSomeNoise(currentINC_SSTDR,ref_SSTDR,currentINC_SSTDR,ref_SSTDR,0,window_left_ss,window_right_ss,nfreq_ss,Z0_ext_ss);
            [~,ext_STDR_FD,~,~,ref_STDR_TD,ref_STDR_FD]...
                = Analysis_AddSomeNoise(currentINC_STDR,ref_STDR,currentINC_STDR,ref_STDR,0,window_left_s-stdr_window_comp,window_right_s+stdr_window_comp,nfreq_s,Z0_ext_s);
        end
    end
end
%% Grab the nonAtt and Att signals for plotting in shortCutForPlottingThings.m
grab_ATT = 1;
if lengthTL == 0 && att == "on" && grab_ATT == 1
    TD_ATT_OC_SSTDR = [];
    TD_ATT_OC_STDR = [];
    TD_ATT_REF_SSTDR = [];
    TD_ATT_REF_STDR = [];
    leg_ATT_OC = [];
    leg_ATT_REF = [];
    Mag_STDR_ATT_ERR = [];
    Phase_STDR_ATT_ERR = [];
    Mag_SSTDR_ATT_ERR = [];
    Phase_SSTDR_ATT_ERR = [];
    leg_ATT_MagPhase_ERR_SSTDR = [];
    leg_ATT_MagPhase_ERR_STDR = [];
end

%% Get attenuated signals
if att == "on"
    Qindex = 2;
    inc_STDR_att = add_SystemEffects(Qindex,inc_STDR,lengthTL,prop_constant_s);
    inc_SSTDR_att = add_SystemEffects(Qindex,inc_SSTDR,lengthTL,prop_constant_ss);


    ref_STDR_att = add_SystemEffects(Qindex,ref_STDR,lengthTL,prop_constant_s);
    ref_SSTDR_att = add_SystemEffects(Qindex,ref_SSTDR,lengthTL,prop_constant_ss);
    
    
    [~,ext_Gamma_STDR_att] = get_LoadFromReflectedSigs(inc_STDR_att,ref_STDR_att,Z0_ext_s,nfreq_s,window_left_s,window_right_s);
    [~,ext_Gamma_SSTDR_att] = get_LoadFromReflectedSigs(inc_SSTDR_att,ref_SSTDR_att,Z0_ext_ss,nfreq_ss,window_left_ss,window_right_ss);

    [ext_SSTDR_FD_att,~,ref_SSTDR_TD_att,ref_SSTDR_FD_att,~,~]...
        = Analysis_AddSomeNoise(inc_SSTDR_att,ref_SSTDR_att,inc_SSTDR_att,ref_SSTDR_att,0,window_left_ss,window_right_ss,nfreq_ss,Z0_ext_ss);
    [~,ext_STDR_FD_att,~,~,ref_STDR_TD_att,ref_STDR_FD_att]...
        = Analysis_AddSomeNoise(inc_STDR_att,ref_STDR_att,inc_STDR_att,ref_STDR_att,0,window_left_s-stdr_window_comp,window_right_s+stdr_window_comp,nfreq_s,Z0_ext_s);
end


%% Grabbing simulation data to plot together in another function/script
% temp_ext_STDR_FD(1,:) = ext_STDR_FD;
% temp_ext_SSTDR_FD(1,:) = ext_SSTDR_FD;

% temp_ext_STDR_FD(2,:) = ext_STDR_FD;
% temp_ext_SSTDR_FD(2,:) = ext_SSTDR_FD;

% temp_ext_STDR_FD(3,:) = ext_STDR_FD;
% temp_ext_SSTDR_FD(3,:) = ext_SSTDR_FD;

% percentBound

%% Plotting flags
plotTitle = 0;          % Fot printing the plot title and load values: plotData = 1, space for Legend = 2
plotLoad_EXT = 1;       % For plotting the extracted FD data 
plotLoad_EXT_ERR = 0;   % For plotting the extracted C error
plotGamma_MagPhase = 0; % For plotting the Mag and Phase of the reflection coefficient
plotGamma_MagPhase_ATT = 0; % For plotting the Mag and Phase of the reflection coefficient
plotZL_EXT = 0;         % For plotting the extracted ZL
plotZL_EXT_ERR = 0;     % For plotting the extracted ZL error
plotZL_MagPhase = 0;    % For plotting the Magnitude and Phase
plotZL_MagPhase_ERR = 0;% For plotting the error in Magnitude and Phase
plotVNA_single_f = 0;   % For checking the VNA response at a specific frequency
plotZL_MMA_ERR = 0;     % For plotting the MMA extracted ZL error
plotZL_MagPhase_MMA = 0;
plotZL_MagPhase_MMA_ERR = 0;         % For plotting the MMA extracted ZL Mag and phase error
plotZL_MagPhase_ATT_ERR = 0;         % For plotting the attenuated ZL Mag and Phase error
plotZL_MagPhase_nPN_ERR = 1;         % For plotting the Mag and Phase of nPN simulation - determine if more PN codes producers more error
plot_TDFD_nPN = 0;      % For plotting all the reflected waveforms using multiple PN codes
plot_TDFD_att = 0;      % For plotting the TDFD data for attenuatied and unattenuated simulations
plot_TDFD = 0;          % For plotting the TD and FD data

% Test of longer extracted ZL figures - Does it actually make a difference?
plotZL_long_EXT = 0;
plotZL_EXT_ERR_long = 0;
plotZL_MagPhase_ERR_long = 0;    % For plotting the error in Magnitude and Phase
plot_TDFD_Long = 0;


if OCSC == 0        % For skipping the OS and SC load extractions but not the TDFD plots
    nLoadStart = 1;     
else
    nLoadStart = 3;
end
twoSubplot = 0;     % For using the default 4 subplot, or combining into 2 subplots
% saveFigs = 1;       % For saving figures
% saveFigPath = 'D:\Google Drive - UofU\U of U\PhD Parts and Papers\Papers\VNA_SSTDR - Simulation\AutoSavedMatlabFigs';

%% MMA Calculation
if plotZL_MMA_ERR == 1 || plotZL_MagPhase_MMA_ERR == 1
    MMAnoise = 5;      % SNR value in dB

    MMA = "on";         % Multi-meas-ave default toggle
    nMMA = 50;          % MMA to average over
    multiple_meas = [{MMA},{nMMA}]; % A cell contianing the instructions for MMA - [T/F, #]


    disp('Performing multiple measurement averaging...');
    for i = nLoadStart:length(load_table.Load)
        Qindex = 3;
        % Average the signals after correlation instead of before
        ref_STDR_mma_load = add_SystemEffects(Qindex,ref_STDR(:,i),multiple_meas{2},MMAnoise);
        ref_SSTDR_mma_load = add_SystemEffects(Qindex,ref_SSTDR(:,i),multiple_meas{2},MMAnoise);

        [ext_ZL_STDR_mma(:,i), ext_Gamma_STDR_mma(:,i)] = get_LoadFromMMA(inc_STDR(:,1),ref_STDR_mma_load,Z0_ext_s,nfreq_s,window_left_s,window_right_s);
        [ext_ZL_SSTDR_mma(:,i), ext_Gamma_SSTDR_mma(:,i)] = get_LoadFromMMA(inc_SSTDR(:,1),ref_SSTDR_mma_load,Z0_ext_ss,nfreq_ss,window_left_ss,window_right_ss);

    end
end

%% PLotting
if plotLoad_EXT == 1    % Needs to use ZL_mat_short to EXT actual R/C
    for nLoad = nLoadStart:length(load_table.Load)
%         legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        for dB = 1:length(dB_vect)

            figure(nextFig);
            subplot(4,1,1)
            plot(f_pos_ss/1e6,abs(real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))),'linewidth',2)
            hold on;

            subplot(4,1,2)
            plot(f_pos_s/1e6,abs(real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))),'linewidth',2)
            hold on;

            c_ext_SSTDR = -1./(imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))'.*f_pos_ss*2*pi);
            c_ext_STDR = -1./(imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))'.*f_pos_s*2*pi);

            subplot(4,1,3)
            semilogy(f_pos_ss/1e6,abs(c_ext_SSTDR),'linewidth',2)
            hold on;

            subplot(4,1,4)
            semilogy(f_pos_s/1e6,abs(c_ext_STDR),'linewidth',2)
            hold on;

            if sum(dB_vect) == 0
                legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                legend_vect_stdr = [legend_vect_stdr, "STDR"];
            else
                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];                    
            end
        end

        figure(nextFig);
        subplot(4,1,1)
        yline(load_table.R(nLoad),':k','linewidth',2)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
%         ylim([load_table.R(nLoad)-.1*(load_table.R(nLoad)+1e-6) load_table.R(nLoad)+.1*(load_table.R(nLoad)+1e-6)])
        ylabel("Resistance [\Omega]")
        xlabel("Frequency [MHz]")
        title(sprintf('%d Ohms',load_table.R(nLoad)))
        legend([legend_vect_sstdr, "Z_{L_{SIM}} = R"])

        subplot(4,1,2)
        yline(load_table.R(nLoad),':k','linewidth',2)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
%         ylim([load_table.R(nLoad)-.1*(load_table.R(nLoad)+1e-6) load_table.R(nLoad)+.1*(load_table.R(nLoad)+1e-6)])
        ylabel("Resistance [\Omega]")
        xlabel("Frequency [MHz]")
        title(sprintf('%d Ohms',load_table.R(nLoad)))
        legend([legend_vect_stdr, "Z_{L_{SIM}} = R"])

        subplot(4,1,3)
        semilogy(f_pos_ss/1e6,load_table.C(nLoad)*ones(length(f_pos_ss),1),':k','linewidth',2)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        ylabel("Capacitance [F]")
        xlabel("Frequency [MHz]")
        legend([legend_vect_sstdr, "{Load_{SIM}} = C"])
        title(sprintf('%d Farads',load_table.C(nLoad)))

        subplot(4,1,4)
        semilogy(f_pos_ss/1e6,load_table.C(nLoad)*ones(length(f_pos_ss),1),':k','linewidth',2)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        ylabel("Capacitance [F]")
        xlabel("Frequency [MHz]")
        legend([legend_vect_stdr, "{Load_{SIM}} = C"])
        title(sprintf('%d Farads',load_table.C(nLoad)))

        sgtitle(sprintf('%sLoad extraction',load_table.Load{nLoad}))

        nextFig = nextFig+1;
    end  
end

if plotLoad_EXT_ERR == 1 % Needs to use ZL_mat_short to get actual ERR in R/C
    for nLoad = nLoadStart:length(load_table.Load)
        if twoSubplot == 1
            legend_vect = [];
            for dB = 1:length(dB_vect)
                c_ext_SSTDR = -1./(imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))'.*f_pos_ss*2*pi);
                c_ext_STDR = -1./(imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))'.*f_pos_s*2*pi);

                r_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad));
                r_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad));
                
                c_sim_SSTDR = -1./(imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)).*f_pos_ss*2*pi);
                c_sim_STDR = -1./(imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)).*f_pos_s*2*pi);
                
                r_sim_SSTDR = real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))';
                r_sim_STDR = real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))';

                SSTDR_err_IM_c_ext = abs((c_ext_SSTDR-c_sim_SSTDR)./c_sim_SSTDR)';
                STDR_err_IM_c_ext = abs((c_ext_STDR-c_sim_STDR)./c_sim_STDR)';

                SSTDR_err_RE = 100*abs((r_ext_SSTDR - r_sim_SSTDR)./r_sim_SSTDR);
                STDR_err_RE = 100*abs((r_ext_STDR - r_sim_STDR)./r_sim_STDR);

                % Grabbing simulation data to plot together in another function/script
    %             temp_ext_STDR_FD_re(1,:) = STDR_err_RE;
    %             temp_ext_STDR_FD_im(1,:) = STDR_err_IM_c_ext;
    %             temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_RE;
    %             temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_IM_c_ext;
    % 
    %             temp_ext_STDR_FD_re(2,:) = STDR_err_RE;
    %             temp_ext_STDR_FD_im(2,:) = STDR_err_IM_c_ext;
    %             temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_RE;
    %             temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_IM_c_ext;

                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_RE,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_RE,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos_ss/1e6,SSTDR_err_IM_c_ext,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_IM_c_ext,'--','linewidth',2)
                hold on;
            
                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR","STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];                   
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Resistance Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Capacitance Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            sgtitle(sprintf('%sLoad Extraction Error',load_table.Load{nLoad}))

            nextFig = nextFig+1;
        else
            legend_vect_sstdr = [];
            legend_vect_stdr = [];
            for dB = 1:length(dB_vect)
                c_ext_SSTDR = -1./(imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))'.*f_pos_ss*2*pi);
                c_ext_STDR = -1./(imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))'.*f_pos_s*2*pi);

                r_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad));
                r_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad));
                
                c_sim_SSTDR = -1./(imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)).*f_pos_ss*2*pi);
                c_sim_STDR = -1./(imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)).*f_pos_s*2*pi);
                
                r_sim_SSTDR = real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))';
                r_sim_STDR = real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))';

                SSTDR_err_IM_c_ext = abs((c_ext_SSTDR-c_sim_SSTDR)./c_sim_SSTDR)';
                STDR_err_IM_c_ext = abs((c_ext_STDR-c_sim_STDR)./c_sim_STDR)';

                SSTDR_err_RE = 100*abs((r_ext_SSTDR - r_sim_SSTDR)./r_sim_SSTDR);
                STDR_err_RE = 100*abs((r_ext_STDR - r_sim_STDR)./r_sim_STDR);
                
                % Grabbing simulation data to plot together in another function/script
    %             temp_ext_STDR_FD_re(1,:) = STDR_err_RE;
    %             temp_ext_STDR_FD_im(1,:) = STDR_err_IM_c_ext;
    %             temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_RE;
    %             temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_IM_c_ext;
    % 
    %             temp_ext_STDR_FD_re(2,:) = STDR_err_RE;
    %             temp_ext_STDR_FD_im(2,:) = STDR_err_IM_c_ext;
    %             temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_RE;
    %             temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_IM_c_ext;

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_RE,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos_s/1e6,STDR_err_RE,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos_ss/1e6,SSTDR_err_IM_c_ext,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos_s/1e6,STDR_err_IM_c_ext,'linewidth',2)
                hold on;

    %             legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];                    
                end
            end

            figure(nextFig);
            subplot(4,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Resistance Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_sstdr)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Resistance Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_stdr)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Capacitance %Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            subplot(4,1,4)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Capacitance %Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            sgtitle(sprintf('%sLoad Extraction Error',load_table.Load{nLoad}))

            nextFig = nextFig+1;
        end
    end
end

if plotZL_EXT == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos_ss/1e6,real_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,real_ZL_ext_STDR,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                semilogy(f_pos_ss/1e6,imag_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                semilogy(f_pos_s/1e6,imag_ZL_ext_STDR,'--','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR","STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect 'Z_{L_{SIM}}'],'Location','best')

            subplot(2,1,2)
            plot(f_pos_ss/1e6,-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(-imag(ZL_mat_short_ss(3*(nfreq_ss/4)-10:3*(nfreq_ss/4)+10,nLoad)));
        %         ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect 'Z_{L_{SIM}}'],'Location','best')
            
            if plotTitle == 1
                subplot(2,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)));
                
                subplot(2,1,2)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('Z_L Extraction for %s',load_table.Load{nLoad}))
            end

            nextFig = nextFig+1;
         
        else
            
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';
        %             imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start:fpos2_end,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos_ss/1e6,real_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos_s/1e6,real_ZL_ext_STDR,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos_ss/1e6,imag_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos_s/1e6,imag_ZL_ext_STDR,'linewidth',2)
                hold on;

%                 legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];                    
                end
            end

            figure(nextFig);
            subplot(4,1,1)
            plot(f_pos_ss/1e6,real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("SSTDR Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect_sstdr 'Z_{L_{SIM}}'])
            
            subplot(4,1,2)
            plot(f_pos_s/1e6,real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("STDR Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect_stdr 'Z_{L_{SIM}}'])
            
            subplot(4,1,3)
            plot(f_pos_ss/1e6,-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(-imag(ZL_mat_short_ss(3*(nfreq_ss/4)-10:3*(nfreq_ss/4)+10,nLoad)));
        %         ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("SSTDR Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect_sstdr 'Z_{L_{SIM}}'])
            
            subplot(4,1,4)
            plot(f_pos_s/1e6,-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        %         ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("STDR Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect_stdr 'Z_{L_{SIM}}'])
            
            if plotTitle == 1
                subplot(4,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))
                
                subplot(4,1,2)
                title(sprintf('%d Ohms',load_table.R(nLoad)))
                
                subplot(4,1,3)
                title(sprintf('%d Farads',load_table.C(nLoad)))
                
                subplot(4,1,4)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction',load_table.Load{nLoad}))
            end

            nextFig = nextFig+1;
        end
    end 
end

if plotGamma_MagPhase == 1
    legend_vect = [];
    for nLoad = nLoadStart:length(load_table.Load)
        
        figure(nextFig)
        subplot(2,1,1)
        plot(f_pos_ss/1e6, abs(Gamma_mat_ss(fpos_start_ss:fpos_end_ss,nLoad)),'linewidth',2)
        hold on;
        
        subplot(2,1,2)
        plot(f_pos_ss/1e6, angle(Gamma_mat_ss(fpos_start_ss:fpos_end_ss,nLoad)),'linewidth',2)
        hold on;
        
        legend_vect = [legend_vect, sprintf("%sLoad = %d%s %d pF",load_list{nLoad,1},load_list{nLoad,2},omega,load_list{nLoad,4}*1e12)];
    end
    
    figure(nextFig)
    subplot(2,1,1)
    xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(sprintf("%s Magnitude",gamma))
    xlabel("Frequency [MHz]")

    subplot(2,1,2)
    xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(sprintf("%s Phase",gamma))
    xlabel("Frequency [MHz]")
    legend(legend_vect)
    
    nextFig = nextFig + 1;
end

if plotGamma_MagPhase_ATT == 1
    legend_vect = [];
    for nLoad = nLoadStart:length(load_table.Load)
        
        figure(nextFig)
        subplot(2,1,1)
        p = plot(f_pos_ss/1e6, abs(Gamma_mat_ss(fpos_start_ss:fpos_end_ss,nLoad)),'--','linewidth',2);
        hold on;
        col = get(p,'Color');
        plot(f_pos_ss/1e6, abs(ext_Gamma_SSTDR_att(fpos_start_ss:fpos_end_ss,nLoad)),'linewidth',2,'Color',col)
        hold on;
        
        subplot(2,1,2)
        plot(f_pos_ss/1e6, angle(Gamma_mat_ss(fpos_start_ss:fpos_end_ss,nLoad)),'--','linewidth',2)
        hold on;
        plot(f_pos_ss/1e6, -angle(ext_Gamma_SSTDR_att(fpos_start_ss:fpos_end_ss,nLoad)),'linewidth',2,'Color',col)
        hold on;
        
        legend_vect = [legend_vect, sprintf("%sLoad = %d%s %d pF",load_list{nLoad,1},load_list{nLoad,2},omega,load_list{nLoad,4}*1e12),sprintf("%sLoad_{ATT} = %d%s %d pF",load_list{nLoad,1},load_list{nLoad,2},omega,load_list{nLoad,4}*1e12)];
    end
    
    figure(nextFig)
    subplot(2,1,1)
    xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(sprintf("%s Magnitude",gamma))
    xlabel("Frequency [MHz]")
    ylim([0 1.25])
    xlim([0 X_FDaxisEnd])

    subplot(2,1,2)
    xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(sprintf("%s Phase",gamma))
    xlabel("Frequency [MHz]")
    ylim([-pi/2 pi/2])
    yticks([-pi/2 0 pi/2])
    yticklabels({ '-\pi/2', '0', '\pi/2'})
    xlim([0 X_FDaxisEnd])
    
    legend(legend_vect)
    
    nextFig = nextFig + 1;
end

if plotZL_EXT_ERR == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                % Grabbing simulation data to plot together in another function/script
%                 if fm == 24e6
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_IM_ZL_ext;
%                     temp_ref_STDR_FD(1,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(1,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,1) = f_pos_ss;
%                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fm == 12e6
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_IM_ZL_ext;
%                     temp_ref_STDR_FD(2,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(2,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,2) = f_pos_ss;
%                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fm == 6e6
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_IM_ZL_ext;
%                     temp_ref_STDR_FD(3,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(3,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,3) = f_pos_ss;
%                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%                 
%% Sample Rate comparison
%                 if fs_ss == 4*fm
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_IM_ZL_ext;
%                     temp_f_pos_ss(:,1) = f_pos_ss;
%                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fs_ss == 16*fm
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_IM_ZL_ext;
%                     temp_f_pos_ss(:,2) = f_pos_ss;
%                     temp_f_pos_s(:,2) = f_pos_s;
% %                 elseif fs_ss == 16*fm
% %                     temp_ext_STDR_FD_re(3,:) = STDR_err_RE_ZL_ext;
% %                     temp_ext_STDR_FD_im(3,:) = STDR_err_IM_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_RE_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_IM_ZL_ext;
% %                     temp_f_pos_ss(:,3) = f_pos_ss;
% %                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%% PN Codde LEngth
%                 if PnL == 5
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_IM_ZL_ext;
%                 elseif PnL == 15
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_RE_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_IM_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_RE_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_IM_ZL_ext;
%                 end

                % Grabbing non-MMA Error to plot with MMA noise
%                 temp_SSTDR_nonMMA_re(nLoad,:) = SSTDR_err_RE_ZL_ext;
%                 temp_STDR_nonMMA_re(nLoad,:) = STDR_err_RE_ZL_ext;
% 
%                 temp_SSTDR_nonMMA_im(nLoad,:) = SSTDR_err_IM_ZL_ext;
%                 temp_STDR_nonMMA_im(nLoad,:) = STDR_err_IM_ZL_ext;
                
                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_RE_ZL_ext,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_IM_ZL_ext,'--','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR", "STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect,'Location','best')

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect,'Location','best')

            if plotTitle == 1
                subplot(2,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,1,2)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('Z_L Extraction Error for %s',load_table.Load{nLoad}))
            end

            nextFig = nextFig+1;
        else        
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos_s/1e6,STDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos_s/1e6,STDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];                    
                end
            end

            figure(nextFig);
            subplot(4,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_sstdr)

            subplot(4,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_stdr)

            subplot(4,1,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)

            subplot(4,1,4)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)

            if plotTitle == 1
                subplot(4,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,2)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,3)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                subplot(4,1,4)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error',load_table.Load{nLoad}))
            end

            nextFig = nextFig+1;
        end
    end
end

if plotZL_MagPhase == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos_ss/1e6,mag_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,mag_ZL_ext_STDR,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos_ss/1e6,phase_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,phase_ZL_ext_STDR,'--','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR","STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
%             ylimBound = mean(abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
%             ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylim([-1 1e3])
            ylabel("Magnitude")
            xlabel("Frequency [MHz]")  
            legend([legend_vect 'Z_{L_{SIM}}'],'Location','best')

            subplot(2,1,2)
            plot(f_pos_ss/1e6,angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
%             ylimBound = mean(-imag(ZL_mat_short_ss(3*(nfreq_ss/4)-10:3*(nfreq_ss/4)+10,nLoad)));
        %         ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylim([-pi pi])
            yticks([-pi -pi/2 0 pi/2 pi])
            yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
            ylabel("Phase")
            xlabel("Frequency [MHz]")
            legend([legend_vect 'Z_{L_{SIM}}'],'Location','best')
            
            if plotTitle == 1
                sgtitle(sprintf('Z_L Mag&Phase: %sLoad = %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
         
        else
            
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';
        %             phase_ZL_ext_STDR = angle(ext_STDR_FD_2(dB,fpos2_start:fpos2_end,nLoad))';

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                figure(nextFig);
                subplot(2,2,1)
                plot(f_pos_ss/1e6,mag_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6,mag_ZL_ext_STDR,'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(f_pos_ss/1e6,phase_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_s/1e6,phase_ZL_ext_STDR,'linewidth',2)
                hold on;

%                 legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];                    
                end
            end

            figure(nextFig);
            subplot(2,2,1)
            plot(f_pos_ss/1e6,abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 1e3])
            xlim([0 X_FDaxisEnd])
            ylabel("SSTDR Magnitude")
            xlabel("Frequency [MHz]") 
            title("(a)")
            
            subplot(2,2,2)
            plot(f_pos_s/1e6,abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',1.5,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 1e3]);
            xlim([0 X_FDaxisEnd/2])
            ylabel("STDR Magnitude")
            xlabel("Frequency [MHz]")
            title("(b)")
            
            subplot(2,2,3)
            plot(f_pos_ss/1e6,angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','right')
            ylim([-pi pi])
            yticks([-pi -pi/2 0 pi/2 pi])
            yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
            xlim([0 X_FDaxisEnd])
            ylabel("SSTDR Phase")
            xlabel("Frequency [MHz]")
            title("(c)")
            legend([legend_vect_sstdr 'Z_{L_{SIM}}'],'orientation','horizontal')
            
            subplot(2,2,4)
            plot(f_pos_s/1e6,angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([-pi pi])
            yticks([-pi -pi/2 0 pi/2 pi])
            yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
            xlim([0 X_FDaxisEnd/2])
            ylabel("STDR Phase")
            xlabel("Frequency [MHz]")
            title("(d)")
            legend([legend_vect_stdr 'Z_{L_{SIM}}'],'orientation','horizontal')
            
            if plotTitle == 1
                sgtitle(sprintf('Z_L Extraction: %sLoad = %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            elseif plotTitle == 2
                sgtitle("")
            end

            nextFig = nextFig+1;
        end
    end 
end

if plotZL_MagPhase_ERR == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';
                
                phase_ZL_sim_SSTDR = angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad));

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

%                 SSTDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_SSTDR+angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%                 STDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_STDR+angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
                SSTDR_err_phase_ZL_ext = abs(phase_ZL_ext_SSTDR-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_phase_ZL_ext = abs(phase_ZL_ext_STDR-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
                

                SSTDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_SSTDR-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_STDR-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                % Grabbing simulation data to plot together in another function/script
%                 if fm == 24e6
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(1,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(1,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,1) = f_pos_ss;
%                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fm == 12e6
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(2,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(2,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,2) = f_pos_ss;
%                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fm == 6e6
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(3,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(3,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,3) = f_pos_ss;
%                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%                 
%% Sample Rate comparison
%                 if fs_ss == 4*fm
% %                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(1,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(1,:) = ref_STDR_FD;
%                     temp_f_pos_ss1 = f_neg_pos_ss;
%                     temp_f_pos_s1 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,1) = f_pos_ss;
% %                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fs_ss == 8*fm
% %                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(2,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(2,:) = ref_STDR_FD;
%                     temp_f_pos_ss2 = f_neg_pos_ss;
%                     temp_f_pos_s2 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,2) = f_pos_ss;
% %                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fs_ss == 16*fm
% %                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(3,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(3,:) = ref_STDR_FD;
%                     temp_f_pos_ss3 = f_neg_pos_ss;
%                     temp_f_pos_s3 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,3) = f_pos_ss;
% %                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%% PN Codde LEngth
%                 if PnL == 5
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 15
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                 end                
                
                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_mag_ZL_ext,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos_s/1e6,STDR_err_phase_ZL_ext,'--','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR", "STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            yline(1,'--','linewidth',2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel([{"Impedance"},{"Magnitude Error"}])
            xlabel("Frequency [MHz]")  
            legend(legend_vect)

            subplot(2,1,2)
            yline(1,'--','linewidth',2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 pi])
            yticks([0 pi/4 pi/2 pi*3/4 pi])
            yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
            ylabel([{"Impedance"},{"Phase Error"}])
            xlabel("Frequency [MHz]")
            legend(legend_vect)

            if plotTitle == 1
                sgtitle(sprintf('Z_L Mag&Phase Error: %sLoad = %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
        else        
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad))';

%                 SSTDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_SSTDR+angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%                 STDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_STDR+angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
                SSTDR_err_phase_ZL_ext = abs(phase_ZL_ext_SSTDR-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
                STDR_err_phase_ZL_ext = abs(phase_ZL_ext_STDR-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));                
                
                SSTDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_SSTDR-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
                STDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_STDR-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));

                % Grabbing non-MMA Error to plot with MMA noise
                temp_SSTDR_nonMMA_re(nLoad,:) = SSTDR_err_mag_ZL_ext;
                temp_STDR_nonMMA_re(nLoad,:) = STDR_err_mag_ZL_ext;

                temp_SSTDR_nonMMA_im(nLoad,:) = SSTDR_err_phase_ZL_ext;
                temp_STDR_nonMMA_im(nLoad,:) = STDR_err_phase_ZL_ext;
                
                figure(nextFig);
                subplot(2,2,1)
                plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext,'linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6,STDR_err_mag_ZL_ext,'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext,'linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_s/1e6,STDR_err_phase_ZL_ext,'linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                    titleSpace = "";
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR_{SNR = %d dB}",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR_{SNR = %d dB}",dB_vect(dB))]; 
                    titleSpace = ["";""];
                end
            end

            figure(nextFig);
            subplot(2,2,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]")  
            title("(a)")

            subplot(2,2,2)
            xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]")  
            title("(b)")

            subplot(2,2,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 pi])
            yticks([0 pi/4 pi/2 pi*3/4 pi])
            yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
%             ylim([0 10])
%             yticks([0 5 10])
%             yticklabels({'0%', '5%', '10%'})
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")  
            title("(c)")
            legend(legend_vect_sstdr)

            subplot(2,2,4)
            xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([0 pi])
            yticks([0 pi/4 pi/2 pi*3/4 pi])
            yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
%             ylim([0 10])
%             yticks([0 5 10])
%             yticklabels({'0%', '5%', '10%'})
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")  
            title("(d)")
            legend(legend_vect_stdr)

            if plotTitle == 1
                subplot(2,2,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,2,2)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,2,3)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                subplot(2,2,4)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Mag&Phase Extraction Error',load_table.Load{nLoad}))
            elseif plotTitle == 2
                sgtitle(titleSpace)
            end

            nextFig = nextFig+1;
        end
    end
end

if plotVNA_single_f == 1
    fm_loc_ss = find(f_neg_pos_ss==fm);     %Finds index of fm
    Z0_vna_fms = Z0_ext_ss(fm_loc_ss);    %Grabs Z0 at fm
    for nLoad = nLoadStart:length(load_table.Load)
        ZL_vna_fm = load_table.R(nLoad) + 1/(1i*2*pi*fm*load_table.C(nLoad));
        Gamma_vna_ss = (ZL_vna_fm-Z0_vna_fms)/(ZL_vna_fm+Z0_vna_fms);
        inc_vna_ss = sin(2*pi*t_axis_ss*fm);
        ref_vna_ss = Gamma_vna_ss*inc_vna_ss;
        GammaEXT_vna_ss = ref_vna_ss/inc_vna_ss; % How do VNA's report a single number when the time vector is long
        ZL_vna_ext = (-Z0_vna_fms.*(GammaEXT_vna_ss+1))./(GammaEXT_vna_ss-1);
        mag_diff_vna = abs(ZL_vna) - abs(ZL_vna_ext);
        phaase_diff_vna = angle(ZL_vna) - angle(ZL_vna_ext);
    end
end

if plotZL_MMA_ERR == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            imag_ZL_ext_SSTDR = imag(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            imag_ZL_ext_STDR = imag(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

            real_ZL_ext_SSTDR = real(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            real_ZL_ext_STDR = real(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

            SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
            STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

            SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
            STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

            % Grabbing MMA Error to plot with non-MMA noise
            temp_SSTDR_MMA_re(nLoad,:) = SSTDR_err_RE_ZL_ext;
            temp_STDR_MMA_re(nLoad,:) = STDR_err_RE_ZL_ext;
            
            temp_SSTDR_MMA_im(nLoad,:) = SSTDR_err_IM_ZL_ext;
            temp_STDR_MMA_im(nLoad,:) = STDR_err_IM_ZL_ext;
            
            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_RE_ZL_ext,'--','linewidth',2)
            hold on;

            subplot(2,1,2)
            plot(f_pos_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_IM_ZL_ext,'--','linewidth',2)
            hold on;

            legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",MMAnoise), sprintf("STDR %d dB Noise",MMAnoise)];

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect,'Location','best')
            set(gca,'yticklabel',["0%", "5%", "10%"])

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect,'Location','best')
            set(gca,'yticklabel',["0%", "5%", "10%"])

            if plotTitle == 1
                subplot(2,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,1,2)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error - %d MMA',load_table.Load{nLoad},nMMA))
            end

            nextFig = nextFig+1;
        else        
%             for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad));
                imag_ZL_ext_STDR = imag(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad));

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD(dB,fpos_start_ss:fpos_end_ss,nLoad));
                real_ZL_ext_STDR = real(ext_STDR_FD(dB,fpos_start_s:fpos_end_s,nLoad));

                SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR'+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR'+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR'-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
                STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR'-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos_s/1e6,STDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos_s/1e6,STDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",MMAnoise)];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",MMAnoise)];
%             end

            figure(nextFig);
            subplot(4,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_sstdr)
            set(gca,'yticklabel',["0%", "5%", "10%"])

            subplot(4,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_stdr)
            set(gca,'yticklabel',["0%", "5%", "10%"])

            subplot(4,1,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)
            set(gca,'yticklabel',["0%", "5%", "10%"])

            subplot(4,1,4)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)
            set(gca,'yticklabel',["0%", "5%", "10%"])

            if plotTitle == 1
                subplot(4,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,2)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,3)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                subplot(4,1,4)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error - %d MMA',load_table.Load{nLoad},nMMA))
            end

            nextFig = nextFig+1;
        end
    end
end

if plotZL_MagPhase_MMA == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            phase_ZL_ext_SSTDR = angle(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            phase_ZL_ext_STDR = angle(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

            mag_ZL_ext_SSTDR = abs(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            mag_ZL_ext_STDR = abs(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

%             SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%             STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
% 
%             SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%             STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

            % Grabbing MMA Error to plot with non-MMA noise
%             temp_SSTDR_MMA_re(nLoad,:) = SSTDR_err_RE_ZL_ext;
%             temp_STDR_MMA_re(nLoad,:) = STDR_err_RE_ZL_ext;
%             
%             temp_SSTDR_MMA_im(nLoad,:) = SSTDR_err_IM_ZL_ext;
%             temp_STDR_MMA_im(nLoad,:) = STDR_err_IM_ZL_ext;
            
            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_RE_ZL_ext,'--','linewidth',2)
            hold on;

            subplot(2,1,2)
            plot(f_pos_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_IM_ZL_ext,'--','linewidth',2)
            hold on;

            legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",MMAnoise), sprintf("STDR %d dB Noise",MMAnoise)];

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect,'Location','best')
            set(gca,'yticklabel',["0%", "5%", "10%"])

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect,'Location','best')
            set(gca,'yticklabel',["0%", "5%", "10%"])

            if plotTitle == 1
                subplot(2,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,1,2)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error - %d MMA',load_table.Load{nLoad},nMMA))
            end

            nextFig = nextFig+1;
        else        
%             for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR_MMA = -angle(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
                phase_ZL_ext_STDR_MMA = -angle(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

                mag_ZL_ext_SSTDR_MMA = abs(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
                mag_ZL_ext_STDR_MMA = abs(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

%                 SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-imag(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%                 STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-imag(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
% 
%                 SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./real(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
%                 STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./real(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

                figure(nextFig);
                subplot(2,2,1)
                plot(f_pos_ss/1e6,abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':')
                hold on;
                plot(f_pos_ss/1e6,mag_ZL_ext_SSTDR_MMA,'linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6,abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':')
                hold on;
                plot(f_pos_s/1e6,mag_ZL_ext_STDR_MMA,'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(f_pos_ss/1e6,angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),':')
                hold on;
                plot(f_pos_ss/1e6,phase_ZL_ext_SSTDR_MMA,'linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_s/1e6,angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)),':')
                hold on;
                plot(f_pos_s/1e6,phase_ZL_ext_STDR_MMA,'linewidth',2)
                hold on;

                legend_vect_sstdr = [legend_vect_sstdr, "ZL_{SIM}", sprintf("SSTDR %d dB Noise",MMAnoise)];                    
                legend_vect_stdr = [legend_vect_stdr, "ZL_{SIM}", sprintf("STDR %d dB Noise",MMAnoise)];
%             end

            figure(nextFig);
            subplot(2,2,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([-1 1e3])
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]")  
%             legend(legend_vect_sstdr)

            subplot(2,2,2)
            xline(fm/1e6,'-',{'f_{pn} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([-1 1e3])
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]")  
%             legend(legend_vect_stdr)

            subplot(2,2,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([-pi pi])
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)

            subplot(2,2,4)
            xline(fm/1e6,'-',{'f_{pn} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([-pi pi])
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)

            if plotTitle == 1
                subplot(4,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,2)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(4,1,3)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                subplot(4,1,4)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error - %d MMA',load_table.Load{nLoad},nMMA))
            end

            nextFig = nextFig+1;
        end
    end

end

if plotZL_MagPhase_MMA_ERR == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        
            phase_ZL_ext_SSTDR = -angle(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            phase_ZL_ext_STDR = -angle(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

            mag_ZL_ext_SSTDR = abs(ext_ZL_SSTDR_mma(fpos_start_ss:fpos_end_ss,nLoad));
            mag_ZL_ext_STDR = abs(ext_ZL_STDR_mma(fpos_start_s:fpos_end_s,nLoad));

            SSTDR_err_phase_ZL_ext = abs((phase_ZL_ext_SSTDR-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))))';
            STDR_err_phase_ZL_ext = abs((phase_ZL_ext_STDR-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))))';
            
%             SSTDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_SSTDR+angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))';
%             STDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_STDR+angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';
%             
            SSTDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_SSTDR-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
            STDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_STDR-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

            % Grabbing MMA Error to plot with non-MMA noise
            temp_SSTDR_MMA_re(nLoad,:) = SSTDR_err_mag_ZL_ext;
            temp_STDR_MMA_re(nLoad,:) = STDR_err_mag_ZL_ext;
            
            temp_SSTDR_MMA_im(nLoad,:) = SSTDR_err_phase_ZL_ext;
            temp_STDR_MMA_im(nLoad,:) = STDR_err_phase_ZL_ext;
            
            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_mag_ZL_ext,'--','linewidth',2)
            hold on;

            subplot(2,1,2)
%             plot(f_pos_ss/1e6,angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)),'--')
%             hold on;
            plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_phase_ZL_ext,'--','linewidth',2)
            hold on;

            legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",MMAnoise), sprintf("STDR %d dB Noise",MMAnoise)];

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 50])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel([{"Impedance"},{"Magnitude Error"}])
            xlabel("Frequency [MHz]")  
            legend(legend_vect)

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 50])
%             ylim([0 10])
%             yticks([0 5 10])
%             yticklabels({'0%', '5%', '10%'})
            ylim([0 pi])
            yticks([0 pi/4 pi/2 pi*3/4 pi])
            yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
            ylabel([{"Impedance"},{"Phase Error"}])
            xlabel("Frequency [MHz]")
            legend(legend_vect)

            if plotTitle == 1
                subplot(2,1,1)
                title(sprintf('%d Ohms',load_table.R(nLoad)))

                subplot(2,1,2)
                title(sprintf('%d Farads',load_table.C(nLoad)))

                sgtitle(sprintf('%sZ_L Extraction Error - %d MMA',load_table.Load{nLoad},nMMA))
            end

            nextFig = nextFig+1;
    end
end

if att == "on" && plotZL_MagPhase_ATT_ERR == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
%         legend_vect_sstdr = [];
%         legend_vect_stdr = [];
        if twoSubplot == 1
            phase_ZL_ext_SSTDR_att = angle(ext_SSTDR_FD_att(1,fpos_start_ss:fpos_end_ss,nLoad))';
            phase_ZL_ext_STDR_att = angle(ext_STDR_FD_att(1,fpos_start_s:fpos_end_s,nLoad))';

            mag_ZL_ext_SSTDR_att = abs(ext_SSTDR_FD_att(1,fpos_start_ss:fpos_end_ss,nLoad))';
            mag_ZL_ext_STDR_att = abs(ext_STDR_FD_att(1,fpos_start_s:fpos_end_s,nLoad))';

            SSTDR_err_phase_ZL_ext_att = abs((phase_ZL_ext_SSTDR_att-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))))';
            STDR_err_phase_ZL_ext_att = abs((phase_ZL_ext_STDR_att-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))))';

            SSTDR_err_mag_ZL_ext_att = 100*abs((mag_ZL_ext_SSTDR_att-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))';
            STDR_err_mag_ZL_ext_att = 100*abs((mag_ZL_ext_STDR_att-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))';

            % Grabbing the Attenuated Mag and Phase for S/SSTDR
            if grab_ATT == 1
                Mag_STDR_ATT_ERR = [Mag_STDR_ATT_ERR STDR_err_mag_ZL_ext_att'];
                Phase_STDR_ATT_ERR = [Phase_STDR_ATT_ERR STDR_err_phase_ZL_ext_att'];
                Mag_SSTDR_ATT_ERR = [Mag_SSTDR_ATT_ERR SSTDR_err_mag_ZL_ext_att'];
                Phase_SSTDR_ATT_ERR = [Phase_SSTDR_ATT_ERR SSTDR_err_phase_ZL_ext_att'];
                leg_ATT_MagPhase_ERR_SSTDR = [leg_ATT_MagPhase_ERR_SSTDR sprintf("SSTDR %s_{TL = %dm}", load_table.Load{nLoad},lengthTL)];
                leg_ATT_MagPhase_ERR_STDR = [leg_ATT_MagPhase_ERR_STDR sprintf("STDR %s_{TL = %dm}", load_table.Load{nLoad},lengthTL)];
            end
            
            % Grabbing simulation data to plot together in another function/script
%                 if fm == 24e6
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(1,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(1,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,1) = f_pos_ss;
%                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fm == 12e6
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(2,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(2,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,2) = f_pos_ss;
%                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fm == 6e6
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(3,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(3,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,3) = f_pos_ss;
%                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%                 
%% Sample Rate comparison
%                 if fs_ss == 4*fm
% %                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(1,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(1,:) = ref_STDR_FD;
%                     temp_f_pos_ss1 = f_neg_pos_ss;
%                     temp_f_pos_s1 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,1) = f_pos_ss;
% %                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fs_ss == 8*fm
% %                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(2,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(2,:) = ref_STDR_FD;
%                     temp_f_pos_ss2 = f_neg_pos_ss;
%                     temp_f_pos_s2 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,2) = f_pos_ss;
% %                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fs_ss == 16*fm
% %                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(3,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(3,:) = ref_STDR_FD;
%                     temp_f_pos_ss3 = f_neg_pos_ss;
%                     temp_f_pos_s3 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,3) = f_pos_ss;
% %                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%% PN Codde LEngth
%                 if PnL == 5
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 8
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 11
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 15
%                     temp_ext_STDR_FD_re(4,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(4,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(4,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(4,:) = SSTDR_err_phase_ZL_ext;
%                 end

            % Grabbing non-MMA Error to plot with MMA noise
%                 temp_SSTDR_nonMMA_re(nLoad,:) = SSTDR_err_mag_ZL_ext;
%                 temp_STDR_nonMMA_re(nLoad,:) = STDR_err_mag_ZL_ext;
% 
%                 temp_SSTDR_nonMMA_im(nLoad,:) = SSTDR_err_phase_ZL_ext;
%                 temp_STDR_nonMMA_im(nLoad,:) = STDR_err_phase_ZL_ext;

            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext_att,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_mag_ZL_ext_att,'--','linewidth',2)
            hold on;

            subplot(2,1,2)
            plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext_att,'linewidth',2)
            hold on;
            plot(f_pos_s/1e6,STDR_err_phase_ZL_ext_att,'--','linewidth',2)
            hold on;

            legend_vect = [legend_vect, "ATT SSTDR", "ATT STDR"];

        end
        
        figure(nextFig);
        subplot(2,1,1)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 X_FDaxisEnd])
        ylim([0 10])
        yticks([0 5 10])
        yticklabels({'0%', '5%', '10%'})
        ylabel([{"Impedance"},{"Magnitude Error"}])
        xlabel("Frequency [MHz]")  
        legend(legend_vect)
%             set(gca,'yticklabel',["0%", "5%", "10%"])

        subplot(2,1,2)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 X_FDaxisEnd])
        ylim([0 pi])
        yticks([0 pi/4 pi/2 pi*3/4 pi])
        yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
        ylabel([{"Impedance"},{"Phase Error"}])
        xlabel("Frequency [MHz]")
        legend(legend_vect)
%             set(gca,'yticklabel',["0%", "5%", "10%"])

        if plotTitle == 1
            subplot(2,1,1)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(2,1,2)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            sgtitle(sprintf('Z_L Mag&Phase Extraction Error w/ATT for %s',load_table.Load{nLoad}))
        end

        nextFig = nextFig+1;
    end
end

nPN_Gray = [175 175 175]/255;
nPN_Green = [0.4660 0.6740 0.1880];%[78 193 37]/255;
nPN_Orange = [216   82   25]/255;
nPN_color = nPN_Gray;
if nPN_onOFF == 1 && plotZL_MagPhase_nPN_ERR == 1
%     S_nPN_mean_err_Mag_mat = [];
%     S_nPN_mean_err_Phase_mat = [];
%     SS_nPN_mean_err_Mag_mat = [];
%     SS_nPN_mean_err_Phase_mat = [];
    for nLoad = nLoadStart:length(load_table.Load)
%         legend_vect = [];
%         legend_vect_sstdr = [];
%         legend_vect_stdr = [];
        S_mag_temp_mat = [];
        SS_mag_temp_mat = [];
        S_phase_temp_mat = [];
        SS_phase_temp_mat = [];
        
        % Get 1PN code data to compare to nPN data
        phase_ZL_ext_SSTDR_1PN = angle(ext_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,nLoad));
        phase_ZL_ext_STDR_1PN = angle(ext_STDR_FD(1,fpos_start_s:fpos_end_s,nLoad));
        
        mag_ZL_ext_SSTDR_1PN = abs(ext_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,nLoad));
        mag_ZL_ext_STDR_1PN = abs(ext_STDR_FD(1,fpos_start_s:fpos_end_s,nLoad));

%         SSTDR_err_phase_ZL_ext_1PN = 100*abs((phase_ZL_ext_SSTDR_1PN' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
%         STDR_err_phase_ZL_ext_1PN = 100*abs((phase_ZL_ext_STDR_1PN' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));
        
        SSTDR_err_phase_ZL_ext_1PN = abs((phase_ZL_ext_SSTDR_1PN' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))));
        STDR_err_phase_ZL_ext_1PN = abs((phase_ZL_ext_STDR_1PN' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))));

        SSTDR_err_mag_ZL_ext_1PN = 100*abs((mag_ZL_ext_SSTDR_1PN'-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
        STDR_err_mag_ZL_ext_1PN = 100*abs((mag_ZL_ext_STDR_1PN'-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));
    
        % Get mean error from ALL PN codes
        phase_ZL_ext_SSTDR_mean = mean(angle(ext_SSTDR_FD_nPN_mat(:,fpos_start_ss:fpos_end_ss,nLoad)));
        phase_ZL_ext_STDR_mean = mean(angle(ext_STDR_FD_nPN_mat(:,fpos_start_s:fpos_end_s,nLoad)));

        mag_ZL_ext_SSTDR_mean = mean(abs(ext_SSTDR_FD_nPN_mat(:,fpos_start_ss:fpos_end_ss,nLoad)));
        mag_ZL_ext_STDR_mean = mean(abs(ext_STDR_FD_nPN_mat(:,fpos_start_s:fpos_end_s,nLoad)));

%         SSTDR_err_phase_ZL_ext_mean = 100*abs((phase_ZL_ext_SSTDR_mean' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
%         STDR_err_phase_ZL_ext_mean = 100*abs((phase_ZL_ext_STDR_mean' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));

        SSTDR_err_phase_ZL_ext_mean = abs((phase_ZL_ext_SSTDR_mean' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))));
        STDR_err_phase_ZL_ext_mean = abs((phase_ZL_ext_STDR_mean' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))));
        
        SSTDR_err_mag_ZL_ext_mean = 100*abs((mag_ZL_ext_SSTDR_mean'-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
        STDR_err_mag_ZL_ext_mean = 100*abs((mag_ZL_ext_STDR_mean'-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));        
        
        xSSTDR = find(SSTDR_err_mag_ZL_ext_1PN <= 2);             % Finds the 1PN data below 2% error
        xSTDR = find(STDR_err_mag_ZL_ext_1PN <= 2);               % Finds the 1PN data below 2% error
        percentBound = [min(xSTDR) max(xSTDR); min(xSSTDR) max(xSSTDR)];   % Percent bounds based on 1PN data [S,S;SS,SS]
        
        % Parse through the PN codes and plot the loads together
        for pn = 1:nPN
            phase_ZL_ext_SSTDR_nPN = angle(ext_SSTDR_FD_nPN_mat(pn,fpos_start_ss:fpos_end_ss,nLoad));
            phase_ZL_ext_STDR_nPN = angle(ext_STDR_FD_nPN_mat(pn,fpos_start_s:fpos_end_s,nLoad));

            mag_ZL_ext_SSTDR_nPN = abs(ext_SSTDR_FD_nPN_mat(pn,fpos_start_ss:fpos_end_ss,nLoad));
            mag_ZL_ext_STDR_nPN = abs(ext_STDR_FD_nPN_mat(pn,fpos_start_s:fpos_end_s,nLoad));

%             SSTDR_err_phase_ZL_ext_nPN = 100*abs((phase_ZL_ext_SSTDR_nPN' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./-angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
%             STDR_err_phase_ZL_ext_nPN = 100*abs((phase_ZL_ext_STDR_nPN' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./-angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));
            
            SSTDR_err_phase_ZL_ext_nPN = abs((phase_ZL_ext_SSTDR_nPN' - angle(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad))));
            STDR_err_phase_ZL_ext_nPN = abs((phase_ZL_ext_STDR_nPN' - angle(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad))));
                        
            SSTDR_err_mag_ZL_ext_nPN = 100*abs((mag_ZL_ext_SSTDR_nPN'-abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)))./abs(ZL_mat_short_ss(fpos_start_ss:fpos_end_ss,nLoad)));
            STDR_err_mag_ZL_ext_nPN = 100*abs((mag_ZL_ext_STDR_nPN'-abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)))./abs(ZL_mat_short_s(fpos_start_s:fpos_end_s,nLoad)));
    
            S_mag_temp_mat(:,pn) = STDR_err_mag_ZL_ext_nPN;
            SS_mag_temp_mat(:,pn) = SSTDR_err_mag_ZL_ext_nPN;
            S_phase_temp_mat(:,pn) = STDR_err_phase_ZL_ext_nPN;
            SS_phase_temp_mat(:,pn) = SSTDR_err_phase_ZL_ext_nPN;
            
            % Grab each PN code in MAT to find MEAN
            
            figure(nextFig)
            subplot(2,2,1)
            plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext_nPN,'linewidth',2,'Color', (nPN_color));
            hold on;

            subplot(2,2,2)
            plot(f_pos_s/1e6,STDR_err_mag_ZL_ext_nPN,'linewidth',2,'Color', (nPN_color));
            hold on;

            subplot(2,2,3)
            h_nPN_ss = plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext_nPN,'linewidth',2,'Color', (nPN_color));
            hold on;

            subplot(2,2,4)
            h_nPN_s = plot(f_pos_s/1e6,STDR_err_phase_ZL_ext_nPN,'linewidth',2,'Color', (nPN_color));
            hold on;
            
%             legend_vect_sstdr = [legend_vect_sstdr sprintf("SSTDR nPN = %d",pn)];
%             legend_vect_stdr = [legend_vect_stdr sprintf("STDR nPN = %d",pn)];   
        end
        % Collect mean vs 1PN error values
%         percentBound = get_bound(SSTDR_err_mag_ZL_ext_1PN,STDR_err_mag_ZL_ext_1PN); % 1 Percent bounds based on 1PN data [S,S;SS,SS]
%         percentBound = [6736,19759;7730,34769]; % nfreq = 2e3 use these bounds [100,400;150,702]; 
        mat_index = log(nPN)/log(2)+1;
        if mat_index == 1
            S_nPN_mean_err_Mag_mat(:,mat_index) = STDR_err_mag_ZL_ext_1PN;
            S_nPN_mean_err_Phase_mat(:,mat_index) = STDR_err_phase_ZL_ext_1PN;
            SS_nPN_mean_err_Mag_mat(:,mat_index) = SSTDR_err_mag_ZL_ext_1PN;
            SS_nPN_mean_err_Phase_mat(:,mat_index) = SSTDR_err_phase_ZL_ext_1PN;
        else        
            S_nPN_mean_err_Mag_mat(:,mat_index) = STDR_err_mag_ZL_ext_mean;
            S_nPN_mean_err_Phase_mat(:,mat_index) = STDR_err_phase_ZL_ext_mean;
            SS_nPN_mean_err_Mag_mat(:,mat_index) = SSTDR_err_mag_ZL_ext_mean;
            SS_nPN_mean_err_Phase_mat(:,mat_index) = SSTDR_err_phase_ZL_ext_mean;
        end
        % Grab the STD from S/SS MagPhase 
        S_mag_std(:,mat_index) = std(S_mag_temp_mat,0,2);
        SS_mag_std(:,mat_index) = std(SS_mag_temp_mat,0,2);
        S_phase_std(:,mat_index) = std(S_phase_temp_mat,0,2);
        SS_phase_std(:,mat_index) = std(SS_phase_temp_mat,0,2);
        
%         max(S_nPN_mean_err_Mag_mat(percentBound(1,1):percentBound(1,2),:))
%         max(S_nPN_mean_err_Mag_mat(percentBound(1,1)-1:percentBound(1,2),:))
        
        legend_vect_sstdr = ["SSTDR Single PN code" sprintf("%d Parallel PN codes",nPN) sprintf("SSTDR Mean %d PN codes",nPN)];
        legend_vect_stdr = ["STDR Single PN code" sprintf("%d Parallel PN codes",nPN) sprintf("STDR Mean %d PN codes",nPN)];
        
        subplot(2,2,1)
        plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext_mean,':','linewidth',1.5,'color', 'k');%,'Color',([0 0.4470 0.7410]));
        hold on;
        plot(f_pos_ss/1e6,SSTDR_err_mag_ZL_ext_1PN,'--','linewidth',2,'Color',([0 0.4470 0.7410]));
        hold on;
%         yline(1,'m:','linewidth',1.5)
        yline(2,'m:','linewidth',1.5)
        xline(f_pos_ss(percentBound(2,1))/1e6,'k:','linewidth',1)
        xline(f_pos_ss(percentBound(2,2))/1e6,'k:','linewidth',1)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 50])
        ylim([0 5])
        yticks([0 2 5])
        yticklabels({'0%', '2%', '5%'})
        ylabel("Magnitude Error")
        xlabel("Frequency [MHz]")
        title("(a)")
        
        subplot(2,2,2)
        plot(f_pos_s/1e6,STDR_err_mag_ZL_ext_mean,':','linewidth',1.5,'color', 'k');%,'Color',([0 0.4470 0.7410]));
        hold on;
        plot(f_pos_s/1e6,STDR_err_mag_ZL_ext_1PN,'--','linewidth',2,'Color',([0 0.4470 0.7410]));
        hold on;
%         yline(1,'m:','linewidth',1.5)
        yline(2,'m:','linewidth',1.5)
        xline(f_pos_s(percentBound(1,1))/1e6,'k:','linewidth',1)
        xline(f_pos_s(percentBound(1,2))/1e6,'k:','linewidth',1)
        xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 25])
        ylim([0 5])
        yticks([0 2 5])
        yticklabels({'0%', '2%', '5%'})
        ylabel("Magnitude Error")
        xlabel("Frequency [MHz]")
        title("(b)")
        
        subplot(2,2,3)
        h_mean_ss = plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext_mean,':','linewidth',1.5,'color', 'k');%,'Color',([0 0.4470 0.7410]));
        hold on;
        h_1PN_ss = plot(f_pos_ss/1e6,SSTDR_err_phase_ZL_ext_1PN,'--','linewidth',2,'Color',([0 0.4470 0.7410]));
        hold on;
%         yline(1,'m:','linewidth',1.5)
%         yline(2,'m:','linewidth',1.5)
        xline(f_pos_ss(percentBound(2,1))/1e6,'k:','linewidth',1)
        xline(f_pos_ss(percentBound(2,2))/1e6,'k:','linewidth',1)
        xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 50])
        ylim([0 pi])
        yticks([0 pi/4 pi/2 pi*3/4 pi])
        yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
        ylabel("Phase Error")
        xlabel("Frequency [MHz]")
        title("(c)")
        legend([h_1PN_ss h_nPN_ss h_mean_ss],legend_vect_sstdr)
        
        subplot(2,2,4)
        h_mean_s = plot(f_pos_s/1e6,STDR_err_phase_ZL_ext_mean,':','linewidth',1.5,'color', 'k');%,'Color',([0 0.4470 0.7410]));
        hold on;
        h_1PN_s = plot(f_pos_s/1e6,STDR_err_phase_ZL_ext_1PN,'--','linewidth',2,'Color',([0 0.4470 0.7410]));
        hold on;
%         yline(1,'m:','linewidth',1.5)
%         yline(2,'m:','linewidth',1.5)
        xline(f_pos_s(percentBound(1,1))/1e6,'k:','linewidth',1)
        xline(f_pos_s(percentBound(1,2))/1e6,'k:','linewidth',1)
        xline(fm/1e6,'-',{'f_{PN} =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        xlim([0 25])
        ylim([0 pi])
        yticks([0 pi/4 pi/2 pi*3/4 pi])
        yticklabels({'0', '1/4\pi', '1/2\pi', '3/4\pi', '\pi'})
        ylabel("Phase Error")
        xlabel("Frequency [MHz]")
        title("(d)")
        legend([h_1PN_s h_nPN_s h_mean_s],legend_vect_stdr)
        
        if plotTitle == 1
            sgtitle(sprintf('Z_L Mag&Phase Extraction Error w/nPN = %d; for %s Load',nPN,load_table.Load{nLoad}));
        elseif plotTitle == 2
            sgtitle(["";""])
        end
    
        nextFig = nextFig +1;
    end
end

if nPN_onOFF == 1 && plot_TDFD_nPN == 1
    leg = [];
    figure(nextFig);
    for nLoad = 1:length(load_table.Load)
        if OCSC == 1
            if nLoad == 1 % Open Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),'--','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),'--','linewidth',2)
                hold on;
            elseif nLoad == 2 % Short Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),':','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),':','linewidth',2)
                hold on;
            else % All the rest
                for pn = 1: nPN
                    subplot(2,2,1)
                    plot(squeeze(ref_STDR_TD_nPN_mat(pn,:,nLoad)),'linewidth',2)
                    hold on;

                    subplot(2,2,2)
                    plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_nPN_mat(pn,:,nLoad))),'linewidth',2)
                    hold on;

                    subplot(2,2,3)
                    plot(squeeze(ref_SSTDR_TD_nPN_mat(pn,:,nLoad)),'linewidth',2)
                    hold on;

                    subplot(2,2,4)
                    plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_nPN_mat(pn,:,nLoad))),'linewidth',2)
                    hold on;


        %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                    leg = [leg sprintf("%s_{PN %d} = 75%s 100pF", load_table.Load{nLoad},pn,omega)];
                end
            end
        else
            for pn = 1:nPN
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD_nPN_mat(pn,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_nPN_mat(pn,:,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD_nPN_mat(pn,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_nPN_mat(pn,:,nLoad))),'linewidth',2)
                hold on;


    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s_{PN %d} = 75%s 100pF", load_table.Load{nLoad},pn,omega)];
            end
        end
    end

    if OCSC == 1
        leg_vect = ["OC", "SC", leg];
    else
        leg_vect = leg;
    end
    
    figure(nextFig);
    subplot(2,2,1)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")

    subplot(2,2,2)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    
    subplot(2,2,3)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")

    subplot(2,2,4)
    legend(leg_vect)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    
    if plotTitle == 1
        sgtitle({'TDFD nPN Plots';sprintf("%dPN codes",pn)})
    end
    nextFig = nextFig+1;
end

% For plotting the time and frequency domain data
if att == "on" && plot_TDFD_att == 1
    leg = [];
    figure(nextFig);
    for nLoad = 1:length(load_table.Load)              
        if OCSC == 1
            if nLoad == 1 % Open Circuit
                
                if grab_ATT == 1
                    TD_ATT_OC_SSTDR = [TD_ATT_OC_SSTDR squeeze(ref_SSTDR_TD_att(1,:,nLoad))'];
                    TD_ATT_OC_STDR = [TD_ATT_OC_STDR squeeze(ref_STDR_TD_att(1,:,nLoad))'];
                    leg_ATT_OC = [leg_ATT_OC sprintf("%s_{TL = %dm}", load_table.Load{nLoad},lengthTL)];
                end
                
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;
                plot(squeeze(ref_STDR_TD_att(1,:,nLoad)),'--','linewidth',2)
                hold on;
                

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),'--','linewidth',2)
                hold on;
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_att(1,:,nLoad))),'--','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;
                plot(squeeze(ref_SSTDR_TD_att(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),'--','linewidth',2)
                hold on;
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_att(1,:,nLoad))),'--','linewidth',2)
                hold on;
            elseif nLoad == 2 % Short Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;
                plot(squeeze(ref_STDR_TD_att(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),':','linewidth',2)
                hold on;
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_att(1,:,nLoad))),':','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;
                plot(squeeze(ref_SSTDR_TD_att(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),':','linewidth',2)
                hold on;
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_att(1,:,nLoad))),':','linewidth',2)
                hold on;
            else % All the rest
                
                if grab_ATT == 1
                    TD_ATT_REF_SSTDR = [TD_ATT_REF_SSTDR squeeze(ref_SSTDR_TD_att(1,:,nLoad))'];
                    TD_ATT_REF_STDR = [TD_ATT_REF_STDR squeeze(ref_STDR_TD_att(1,:,nLoad))'];
%                     leg_ATT_REF = [leg_ATT_REF sprintf("%s_{TL = %dm} = 75%s 100pF", load_table.Load{nLoad},lengthTL,omega)];
                    leg_ATT_REF = [leg_ATT_REF sprintf("%s_{TL = %dm}", load_table.Load{nLoad},lengthTL)];
                end
                
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'linewidth',2)
                hold on;
                plot(squeeze(ref_STDR_TD_att(1,:,nLoad)),'linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_att(1,:,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'linewidth',2)
                hold on;
                plot(squeeze(ref_SSTDR_TD_att(1,:,nLoad)),'linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_att(1,:,nLoad))),'linewidth',2)
                hold on;

                
    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s = 75%s 100pF", load_table.Load{nLoad},omega) sprintf("%s_{ATT} = 75%s 100pF", load_table.Load{nLoad},omega)];
            end
        else
            subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'-','linewidth',2)
                hold on;
                plot(squeeze(ref_STDR_TD_att(1,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;
                plot(f_neg_pos_s/1e6, abs(squeeze(ref_STDR_FD_att(1,:,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'-','linewidth',2)
                hold on;
                plot(squeeze(ref_SSTDR_TD_att(1,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;
                plot(f_neg_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD_att(1,:,nLoad))),'linewidth',2)
                hold on;

               
    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s = 75%s 100pF", load_table.Load{nLoad},omega) sprintf("%s_{ATT} = 75%s 100pF", load_table.Load{nLoad},omega)];
        end
    end

    if OCSC == 1
        leg_vect = ["OC", "OC_{ATT}", "SC", "SC_{ATT}", leg];
    else
        leg_vect = leg;
    end

    figure(nextFig);
    subplot(2,2,1)
    % xlim([0 12])
%     legend(leg_vect)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")
%     axis off

    subplot(2,2,2)
%     legend(leg_vect)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
%     ylim([0 max(abs(squeeze(ref_STDR_FD(1,:,nLoad))))])
%     axis off

    subplot(2,2,3)
    % xlim([0 12])
%     legend(leg_vect)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")
%     axis off

    subplot(2,2,4)
    legend(leg_vect)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
%     axis off

    if plotTitle == 1
        sgtitle({sprintf('w/ATT - %d MHz Modulation Frequency',fm/1e6);sprintf("%d bit PN code",PnL)})
    end
    nextFig = nextFig+1;
end


if plotZL_long_EXT == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';
    %             imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start:fpos2_end,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos2_ss/1e6,real_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,real_ZL_ext_STDR,'linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos2_ss/1e6,imag_ZL_ext_SSTDR,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,imag_ZL_ext_STDR,'linewidth',2)
                hold on;

                legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];
            end

            figure(nextFig);
            subplot(2,1,1)
            plot(f_pos2_ss/1e6,real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("Z_L Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect 'Z_{L_{SIM}}'])

            subplot(2,1,2)
            plot(f_pos2_ss/1e6,-imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(-imag(ZL_mat_ss(3*(nfreq2_ss/4)-10:3*(nfreq2_ss/4)+10,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("Z_L Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect 'Z_{L_{SIM}}'])

            if plotTitle == 1
                sgtitle(sprintf('%sZ_L_{long} Extraction - %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
        else        
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';
    %             imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start:fpos2_end,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos2_ss/1e6,real_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos2_s/1e6,real_ZL_ext_STDR,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos2_ss/1e6,imag_ZL_ext_SSTDR,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos2_s/1e6,imag_ZL_ext_STDR,'linewidth',2)
                hold on;

                legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];
            end

            figure(nextFig);
            subplot(4,1,1)
            plot(f_pos2_ss/1e6,real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("SSTDR Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect_sstdr 'Z_{L_{SIM}}'])
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,2)
            plot(f_pos2_s/1e6,real(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("STDR Real")
            xlabel("Frequency [MHz]")  
            legend([legend_vect_stdr 'Z_{L_{SIM}}'])
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,3)
            plot(f_pos2_ss/1e6,-imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylimBound = mean(-imag(ZL_mat_ss(3*(nfreq2_ss/4)-10:3*(nfreq2_ss/4)+10,nLoad)));
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("SSTDR Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect_sstdr 'Z_{L_{SIM}}'])
            title(sprintf('%d Farads',load_table.C(nLoad)))

            subplot(4,1,4)
            plot(f_pos2_s/1e6,-imag(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)),':k','linewidth',2);
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([ylimBound - .1*(ylimBound + 1) ylimBound + .1*(ylimBound + 1)]);
            ylabel("STDR Imaginary")
            xlabel("Frequency [MHz]")
            legend([legend_vect_stdr 'Z_{L_{SIM}}'])
            title(sprintf('%d Farads',load_table.C(nLoad)))

            sgtitle(sprintf('%sZ_L Extraction - Long Z_L',load_table.Load{nLoad}))

            nextFig = nextFig+1;
        end
    end 
end

if plotZL_EXT_ERR_long == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./-imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./-imag(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./real(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos2_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,STDR_err_RE_ZL_ext,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos2_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,STDR_err_IM_ZL_ext,'--','linewidth',2)
                hold on;

                legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];
            end

            figure(nextFig);
            subplot(2,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            xlim([0 X_FDaxisEnd])
            ylabel("Real Z_L Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect)

            subplot(2,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            xlim([0 X_FDaxisEnd])
            ylabel("Imaginary Z_L Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect)
            
            if plotTitle == 1
                sgtitle(sprintf('%sZ_{L_{Long}} Extraction Error - %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
        else
            for dB = 1:length(dB_vect)
                imag_ZL_ext_SSTDR = imag(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                imag_ZL_ext_STDR = imag(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                real_ZL_ext_SSTDR = real(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                real_ZL_ext_STDR = real(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                SSTDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_SSTDR+imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./-imag(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_IM_ZL_ext = 100*abs((imag_ZL_ext_STDR+imag(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./-imag(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                SSTDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_SSTDR-real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./real(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_RE_ZL_ext = 100*abs((real_ZL_ext_STDR-real(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./real(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                figure(nextFig);
                subplot(4,1,1)
                plot(f_pos2_ss/1e6,SSTDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,2)
                plot(f_pos2_s/1e6,STDR_err_RE_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,3)
                plot(f_pos2_ss/1e6,SSTDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                subplot(4,1,4)
                plot(f_pos2_s/1e6,STDR_err_IM_ZL_ext,'linewidth',2)
                hold on;

                legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR %d dB Noise",dB_vect(dB))];                    
                legend_vect_stdr = [legend_vect_stdr, sprintf("STDR %d dB Noise",dB_vect(dB))];
            end

            figure(nextFig);
            subplot(4,1,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_sstdr)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Real Error")
            xlabel("Frequency [MHz]")  
            legend(legend_vect_stdr)
            title(sprintf('%d Ohms',load_table.R(nLoad)))

            subplot(4,1,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary %Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            subplot(4,1,4)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            ylim([0 10])
            ylabel("Imaginary %Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)
            title(sprintf('%d Farads',load_table.C(nLoad)))

            sgtitle(sprintf('%sZ_L Extraction Error - Long Z_L',load_table.Load{nLoad}))

            nextFig = nextFig+1;
        end
    end
end

if plotZL_MagPhase_ERR_long == 1
    for nLoad = nLoadStart:length(load_table.Load)
        legend_vect = [];
        legend_vect_sstdr = [];
        legend_vect_stdr = [];
        if twoSubplot == 1
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                SSTDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_SSTDR+angle(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./-angle(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_STDR+angle(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./-angle(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                SSTDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_SSTDR-abs(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./abs(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_STDR-abs(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./abs(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                % Grabbing simulation data to plot together in another function/script
%                 if fm == 24e6
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(1,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(1,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,1) = f_pos_ss;
%                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fm == 12e6
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(2,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(2,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,2) = f_pos_ss;
%                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fm == 6e6
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_ref_STDR_FD(3,:) = squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,1));
%                     temp_ref_SSTDR_FD(3,:) = squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,1));
%                     temp_f_pos_ss(:,3) = f_pos_ss;
%                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%                 
%% Sample Rate comparison
%                 if fs_ss == 4*fm
% %                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(1,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(1,:) = ref_STDR_FD;
%                     temp_f_pos_ss1 = f_neg_pos_ss;
%                     temp_f_pos_s1 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,1) = f_pos_ss;
% %                     temp_f_pos_s(:,1) = f_pos_s;
%                 elseif fs_ss == 8*fm
% %                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(2,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(2,:) = ref_STDR_FD;
%                     temp_f_pos_ss2 = f_neg_pos_ss;
%                     temp_f_pos_s2 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,2) = f_pos_ss;
% %                     temp_f_pos_s(:,2) = f_pos_s;
%                 elseif fs_ss == 16*fm
% %                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
% %                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
% %                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
% %                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                     temp_SSTDR_FD(3,:) = ref_SSTDR_FD;
%                     temp_STDR_FD(3,:) = ref_STDR_FD;
%                     temp_f_pos_ss3 = f_neg_pos_ss;
%                     temp_f_pos_s3 = f_neg_pos_s;
% %                     temp_f_pos_ss(:,3) = f_pos_ss;
% %                     temp_f_pos_s(:,3) = f_pos_s;
%                 end
%% PN Codde LEngth
%                 if PnL == 5
%                     temp_ext_STDR_FD_re(1,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(1,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(1,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(1,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 8
%                     temp_ext_STDR_FD_re(2,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(2,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(2,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(2,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 11
%                     temp_ext_STDR_FD_re(3,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(3,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(3,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(3,:) = SSTDR_err_phase_ZL_ext;
%                 elseif PnL == 15
%                     temp_ext_STDR_FD_re(4,:) = STDR_err_mag_ZL_ext;
%                     temp_ext_STDR_FD_im(4,:) = STDR_err_phase_ZL_ext;
%                     temp_ext_SSTDR_FD_re(4,:) = SSTDR_err_mag_ZL_ext;
%                     temp_ext_SSTDR_FD_im(4,:) = SSTDR_err_phase_ZL_ext;
%                 end

                % Grabbing non-MMA Error to plot with MMA noise
%                 temp_SSTDR_nonMMA_re(nLoad,:) = SSTDR_err_mag_ZL_ext;
%                 temp_STDR_nonMMA_re(nLoad,:) = STDR_err_mag_ZL_ext;
% 
%                 temp_SSTDR_nonMMA_im(nLoad,:) = SSTDR_err_phase_ZL_ext;
%                 temp_STDR_nonMMA_im(nLoad,:) = STDR_err_phase_ZL_ext;
                
                figure(nextFig);
                subplot(2,1,1)
                plot(f_pos2_ss/1e6,SSTDR_err_mag_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,STDR_err_mag_ZL_ext,'--','linewidth',2)
                hold on;

                subplot(2,1,2)
                plot(f_pos2_ss/1e6,SSTDR_err_phase_ZL_ext,'linewidth',2)
                hold on;
                plot(f_pos2_s/1e6,STDR_err_phase_ZL_ext,'--','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect = [legend_vect, "SSTDR", "STDR"];
                else
                    legend_vect = [legend_vect, sprintf("SSTDR %d dB Noise",dB_vect(dB)), sprintf("STDR %d dB Noise",dB_vect(dB))];
                end
            end

            figure(nextFig);
            subplot(2,1,1)
            yline(1,'--','linewidth',2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel([{"Impedance"},{"Magnitude Error"}])
            xlabel("Frequency [MHz]")  
            legend(legend_vect,'Location','best')

            subplot(2,1,2)
            yline(1,'--','linewidth',2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel([{"Impedance"},{"Phase Error"}])
            xlabel("Frequency [MHz]")
            legend(legend_vect,'Location','best')

            if plotTitle == 1
                sgtitle(sprintf('%sZ_{L_{Long}} Mag&Phase Extraction Error - %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
        else        
            for dB = 1:length(dB_vect)
                phase_ZL_ext_SSTDR = angle(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                phase_ZL_ext_STDR = angle(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                mag_ZL_ext_SSTDR = abs(ext_SSTDR_FD_2(dB,fpos2_start_ss:fpos2_end_ss,nLoad))';
                mag_ZL_ext_STDR = abs(ext_STDR_FD_2(dB,fpos2_start_s:fpos2_end_s,nLoad))';

                SSTDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_SSTDR+angle(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./angle(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_phase_ZL_ext = 100*abs((phase_ZL_ext_STDR+angle(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))./angle(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                SSTDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_SSTDR-abs(ZL_mat_short_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))./abs(ZL_mat_ss(fpos2_start_ss:fpos2_end_ss,nLoad)))';
                STDR_err_mag_ZL_ext = 100*abs((mag_ZL_ext_STDR-abs(ZL_mat_short_s(fpos2_start_s:fpos2_end_s,nLoad)))./abs(ZL_mat_s(fpos2_start_s:fpos2_end_s,nLoad)))';

                figure(nextFig);
                subplot(2,2,1)
                plot(f_pos2_ss/1e6,SSTDR_err_mag_ZL_ext','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos2_s/1e6,STDR_err_mag_ZL_ext','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(f_pos2_ss/1e6,SSTDR_err_phase_ZL_ext','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos2_s/1e6,STDR_err_phase_ZL_ext','linewidth',2)
                hold on;

                if sum(dB_vect) == 0
                    legend_vect_sstdr = [legend_vect_sstdr, "SSTDR"];                    
                    legend_vect_stdr = [legend_vect_stdr, "STDR"];
                else
                    legend_vect_sstdr = [legend_vect_sstdr, sprintf("SSTDR_{SNR = %d dB}",dB_vect(dB))];                    
                    legend_vect_stdr = [legend_vect_stdr, sprintf("STDR_{SNR = %d dB}",dB_vect(dB))];                    
                end
            end

            figure(nextFig);
            subplot(2,2,1)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]")

            subplot(2,2,2)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Magnitude Error")
            xlabel("Frequency [MHz]") 

            subplot(2,2,3)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_sstdr)

            subplot(2,2,4)
            xline(fm/1e6,'-',{'f_m =', sprintf('%d MHz',fm/1e6)},'linewidth',2,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            xlim([0 X_FDaxisEnd/2])
            ylim([0 10])
            yticks([0 5 10])
            yticklabels({'0%', '5%', '10%'})
            ylabel("Phase Error")
            xlabel("Frequency [MHz]")
            legend(legend_vect_stdr)

            if plotTitle == 1
                sgtitle(sprintf('%sZ_{L_{Long}} Mag&Phase Extraction Error - %d%s %dpF',load_table.Load{nLoad},load_table.R(nLoad),omega,load_table.C(nLoad)*1e12))
            end

            nextFig = nextFig+1;
        end
    end
end

% For plotting the time and frequency domain data
if plot_TDFD == 1
    leg = [];
    figure(nextFig);
    for nLoad = 1:length(load_table.Load)    
        if OCSC == 1
            if nLoad == 1 % Open Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,nLoad))),'--','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,nLoad))),'--','linewidth',2)
                hold on;
            elseif nLoad == 2 % Short Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,nLoad))),':','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,nLoad))),':','linewidth',2)
                hold on;
            else % All the rest
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,fpos_start_s:fpos_end_s,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,fpos_start_ss:fpos_end_ss,nLoad))),'linewidth',2)
                hold on;

                
    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s = 75%s 100pF", load_table.Load{nLoad},omega)];
            end
        else
            subplot(2,2,1)
                plot(squeeze(ref_STDR_TD(1,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_pos_s/1e6, abs(squeeze(ref_STDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD(1,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_pos_ss/1e6, abs(squeeze(ref_SSTDR_FD(1,:,nLoad))),'linewidth',2)
                hold on;

                
    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s = 75%s 100pF", load_table.Load{nLoad},omega)];
        end
    end

    if OCSC == 1
        leg_vect = ["OC", "SC", leg];
    else
        leg_vect = leg;
    end
    
    figure(nextFig);
    subplot(2,2,1)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")
    title("(a)")

    subplot(2,2,2)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    title("(b)")
    
    subplot(2,2,3)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")
    title("(c)")

    subplot(2,2,4)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    title("(d)")
    legend(leg_vect,'orientation','horizontal','location','north')
%     legend(leg_vect,'Position',[0.2 0.2 0.1 0.2],'Box','off')
%     title(["";""])

%     legend(leg_vect,'Location','EastOutside','Orientation','vertical','Box','off')

    if plotTitle == 1
        sgtitle({sprintf('%d MHz Modulation Frequency',fm/1e6);sprintf("%d bit PN code",PnL)})
    elseif plotTitle == 2
        sgtitle("")
    end
    nextFig = nextFig+1;
end

if plot_TDFD_Long == 1
    leg = [];
    figure(nextFig);
    for nLoad = 1:length(load_table.Load)
        if OCSC == 1
            if nLoad == 1 % Open Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD_2(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos2_s/1e6, abs(squeeze(ref_STDR_FD_2(1,:,nLoad))),'--','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD_2(1,:,nLoad)),'--','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos2_ss/1e6, abs(squeeze(ref_SSTDR_FD_2(1,:,nLoad))),'--','linewidth',2)
                hold on;
            elseif nLoad == 2 % Short Circuit
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD_2(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos2_s/1e6, abs(squeeze(ref_STDR_FD_2(1,:,nLoad))),':','linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD_2(1,:,nLoad)),':','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos2_ss/1e6, abs(squeeze(ref_SSTDR_FD_2(1,:,nLoad))),':','linewidth',2)
                hold on;
            else % All the rest
                for dB = 1:length(dB_vect)
                    subplot(2,2,1)
                    plot(squeeze(ref_STDR_TD_2(dB,:,nLoad)),'linewidth',2)
                    hold on;

                    subplot(2,2,2)
                    plot(f_neg_pos2_s/1e6, abs(squeeze(ref_STDR_FD_2(dB,:,nLoad))),'linewidth',2)
                    hold on;

                    subplot(2,2,3)
                    plot(squeeze(ref_SSTDR_TD_2(dB,:,nLoad)),'linewidth',2)
                    hold on;

                    subplot(2,2,4)
                    plot(f_neg_pos2_ss/1e6, abs(squeeze(ref_SSTDR_FD_2(dB,:,nLoad))),'linewidth',2)
                    hold on;


        %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                    leg = [leg sprintf("%s_{%ddB Noise} = %d%s %dpF",load_table.Load{nLoad},dB_vect(dB),load_table.R(nLoad),omega,load_table.C(nLoad)*1e12)];
                end
            end
        else
            for dB = 1:length(dB_vect)
                subplot(2,2,1)
                plot(squeeze(ref_STDR_TD_2(dB,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,2)
                plot(f_neg_pos2_s/1e6, abs(squeeze(ref_STDR_FD_2(dB,:,nLoad))),'linewidth',2)
                hold on;

                subplot(2,2,3)
                plot(squeeze(ref_SSTDR_TD_2(dB,:,nLoad)),'-','linewidth',2)
                hold on;

                subplot(2,2,4)
                plot(f_neg_pos2_ss/1e6, abs(squeeze(ref_SSTDR_FD_2(dB,:,nLoad))),'linewidth',2)
                hold on;


    %             leg = [leg sprintf("%s = [%d %d %d]", load_table.Load{nLoad},load_table.R(nLoad),load_table.L(nLoad),load_table.C(nLoad))];
                leg = [leg sprintf("%s_{%ddB Noise} = %d%s %dpF",load_table.Load{nLoad},dB_vect(dB),load_table.R(nLoad),omega,load_table.C(nLoad)*1e12)];
            end
        end
    end


    if OCSC == 1
        leg_vect = ["OC", "SC", leg];
    else
        leg_vect = leg;
    end
    
    figure(nextFig);
    subplot(2,2,1)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")

    subplot(2,2,2)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    
    subplot(2,2,3)
    xlabel("Sample Points")
    ylabel("Correlation Magnitude")

    subplot(2,2,4)
    legend(leg_vect)
    xlabel("Frequency [MHz]")
    ylabel("Frequency Magnitude")
    
    if plotTitle == 1
        sgtitle({sprintf('%d MHz Modulation Frequency - Long',fm/1e6);sprintf("%d bit PN code",PnL)})
    end
    nextFig = nextFig+1;
end
















