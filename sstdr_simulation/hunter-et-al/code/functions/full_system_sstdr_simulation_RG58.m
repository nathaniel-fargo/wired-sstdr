%==========================================================================
% Full System Simulation Using the SSP method
% The creates a simulated SSTDR response of a system containing an n number
% of loads with any impedance. The loads can be arranged in series of
% parallel. Rg58 cable is the transmission line that is modeled.
%
% inputs:
%   num_z = the number of elements to be modeled (integer)
% 
%   imped = the impedance of the elements (array with num_z number of
%   columns and 40960 rows)
% 
%   len = distance between elements in meters. The first len entry is for the
%   distance from the SSTDR to the closest element, and the second entry is
%   for the distance from the element closest to the SSTDR to the next closest
%   element... (vector of length num_z)
% 
%   s_p = shows if the element is in series or parallel. This can change
%   for each element. (vector of length num_z - 1) write 's' for series and
%   'p' for parallel. the first entry corresponds to the element closest to
%   the SSTDR and the second corresponds to the next closest element...
%   there is no entry for the final element because this element can not be
%   in series or parallel its just at the end of the transmission line
%   (assumed)
%
%   Ns = number of data points to be outputted (integer)
% 
% outputs:
%   ref_sig = the normalized reflection signature (SSTDR response) of the system under test
%   (vector of length Ns (must not exceed 1500)) ((y data))
% 
%   dist_m = the corresponding distance in meters for the normalized reflection
%   signature (vector of length Ns (must not exceed 1500)) ((x data))
%==========================================================================
% Hunter Ellis
% Department of Electrical and Computer Engineering
% University of Utah, Salt Lake City, Utah
% Last Updated June 16, 2020
 
function [ref_sig, dist_m] = full_system_sstdr_simulation_RG58(num_z, imped, len, s_p, Ns) 
 
%==========================================================================
% Step 0: Define SSTDR paramiters
%==========================================================================
 
Fsampl = 2*48e6;        %sampling rate at 48 MHz
Fmodu  = Fsampl/4;      % Modulation rate
Fchip  = Fsampl/4;      % Chip rate
PNL    = 10240;         % Number of chips
                     %   THIS IS FOR PLOTTING PURPOSES ONLY
                     Ns = 1590;            % Number of samples to plot
 
amp = 1;
shft = 10;           % How much to shift signal to left to visualize (so that the signal is not centered around 0)
                     %   THIS IS FOR PLOTTING PURPOSES ONLY
 
 
%==========================================================================
% Step 1: create discretized frequency axis. This is only needed if you are
% not using LCR data for the component impedance.
%==========================================================================
Q = PNL*Fsampl/(Fchip);         % Length of simulated signal
n = 1:Q;                        % Sample axis
t = n/Fsampl;                   % Time axis
r = floor(Q/2)+1;               % Center sample in frequency
f = ifftshift((n-r)/Q)*Fsampl;  % Frequency axis
nshft = circshift(ifftshift((n-r)), shft);      % Shifted sample axis
tshft = nshft/Fsampl;                           % Sifted time axis
sqmod = [1 -1].';                                               % Square wave 
sqmod_long = kron(sqmod, ones(Fsampl/Fmodu/2,1));               % Square wave (expanded out based on modulation rate)
scorr_sig = ifft(fft(sqmod_long,Q).*conj(fft(sqmod_long,Q)))/(Fsampl/(Fchip));   
scorr_sig_shift = circshift(scorr_sig, shft); 
 
 
%==========================================================================
% Step 2: Define transmission line parameters for Rg58 cable.
%==========================================================================
rs = 1/(5.8*10^7);                   % impedance of copper
a=.445*10^-3;                        % inner radius of RG58
b=1.765*10^-3;                       %outer radius of RG58
 
%lumped element piecewise values for RG58
r = (rs/(2*pi))*((1/a)+(1/b));      %resistance of RG58 per m
l = 2.7557e-7;                      %inductance of RG58 per m
g = 4.5602e-4;                      % admittance of RG58 per m
c = 9.1164e-11;                     %capacitance of rg58 per m
vop = .641;                          % vop of rg58
x = r+(1i*2*pi*f(:)*l);               %intermediate calculation
y = g+(1i*2*pi*f(:)*c);             %intermediate calculation
gam = sqrt(x.*y);                   % calculating gamma 
z0 = 54.976-(1i*.0368);             %Z0 for rg58
len = len.*3.28084;                  % lengths of the different transmission lines
 
 
 
%==========================================================================
% Step 5: Define number of loads, and load impedances
%==========================================================================
num_resis = num_z;                      %number of components separated by T-lines
sstdr_imped = 68;
zl = imped; %load impedance for resistor
%==========================================================================
% Step 6: define variables for simulation
%==========================================================================
TX = 1;                              %initialize current transfer function
Z_equiv = 0;                        %initialize equivalent resistance
 
 
%==========================================================================
% Step 7: create simulation using SSP method
%==========================================================================
 beta = 2*pi*f(:)*sqrt(l*c);         % betta of T-line
alpha =  .5*((r/z0) + (g*z0));      %alpha of T-line
for i = 1:length(f)                 %loop over frequency
    Z_equiv = 0; %reset equivalent impedance every time you retest system with new freq
    for j = num_resis:-1:1 %loop over number of elements separated by T-lines
        if j == num_resis %if you are at the beginning of transmission line (furthest from measuring device) do this
            gamma1 = (zl(i,j)-z0)/(zl(i,j)+z0); %find reflection coefficient
        else                %if you are not at the beginning of transmission line do this
            
            if s_p(j) == 'p'
                TX = zl(i,j)/(Z_equiv+zl(i,j)); %calculate current transfer function
                Z_equiv = TX*Z_equiv;       %find equivalent impedance
            end
            
            if s_p(j) == 's'
                TX = Z_equiv/(Z_equiv+zl(i,j)); %calculate voltage transfer function
                Z_equiv = Z_equiv/TX;       %find equivalent impedance
            end
            
        gamma1 = (Z_equiv-z0)/(Z_equiv+z0); %find reflection coefficient using equivalent impedance
        end
        Gamma2(i)=gamma1*exp(2*(alpha+(1i*beta(i)))*(-len(j)*0.3048)); % move reflection coeficient through T-line
        Z_equiv = z0*(Gamma2(i)+1)/(1-Gamma2(i)); %calculate new equivalent impedance
    end
    Gamma3(i) = (Z_equiv-sstdr_imped)/(Z_equiv+sstdr_imped); %reflection coefficient at the input of the SSTDR
end
 
Gamma = transpose(Gamma3); %convert to row vector B/c Dr. Harleys code makes a row vector
scorr_y = real(ifft(fft(scorr_sig,[],1).*Gamma,[],1)).*amp(:).';  % multiply signal by response in freq domain then take ifft
scorr_y_shift = circshift(scorr_y,shft,1);              % FOR PLOTTING PURPOSES ONLY
 
 
 
 
sample_size=.2;  % This variable is used as the sampling size for spline interpolation  %
 
count = 0;
    tshft = resample(tshft,5,1);
 dist = tshft*3.28084*2.99792*1e8*vop*.5; %convert x axis to distance
for k = 1:Ns
    if dist(k) <0
        count = count + 1;
    end
end
 
 
Gamma = transpose(Gamma3); 
scorr_y = real(ifft(fft(scorr_sig,[],1).*Gamma,[],1)).*amp(:).';  % multiply signal by response in freq domain then take ifft
scorr_y_shift = circshift(scorr_y,shft,1);              % FOR PLOTTING PURPOSES ONLY
x = 1:length(scorr_y_shift);
    t_ip=1:sample_size:length(scorr_y_shift);%% This is the new interpolated samples for the x-axis of SSTDR data %%
    ip_raw_data=interp1(x,scorr_y_shift,t_ip,'spline'); %% First Spline Interpolation on the Raw data %%
[x_sim,y_sim,~,~,~,~]=slope_calc(tshft(count:Ns+count),ip_raw_data(count:Ns+count)); % Smooth reflection signature 
 
dist_m = dist(count:Ns+count).*.3048;
ref_sig = y_sim./sstdr_sim_open(sum(len));              % Normalize the reflection coefficient
end
