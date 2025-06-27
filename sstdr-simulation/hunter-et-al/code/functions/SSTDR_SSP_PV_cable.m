%==========================================================================
% Full System Simulation Using the SSP method
% The creates a simulated SSTDR response of a system containing an n number
% of loads with any impedance. The loads can be arranged in series of
% parallel. PV cable is the transmission line that is modeled.
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
%   there is no entry for the final element because this element cannot be
%   in series or parallel it’s just at the end of the transmission line
%   (assumed)
%
%   Ns = number of data points to be outputted (integer)
% 
% outputs:
%   ref_sig = the normalized reflection signature (SSTDR response) of the system under test
%   (vector of length Ns (must not exceed 20480)) ((y data))
% 
%   dist_m = the corresponding distance in meters for the normalized reflection
%   signature (vector of length Ns (must not exceed 20480)) ((x data))
%==========================================================================
% Hunter Ellis
% Department of Electrical and Computer Engineering
% University of Utah, Salt Lake City, Utah
% Last Updated June 16, 2020


function [ref_sig, dist_m] = SSTDR_SSP_PV_cable(num_z, imped, len, s_p, Ns) 

%==========================================================================
% Step 0: Define SSTDR parameters
%==========================================================================
Fsampl = 4*24e6;          % Sampling rate
Fmodu  = Fsampl/4;        % Modulation rate
Fchip  = Fsampl/4;        % Chip rate
PNL    = 10240;           % Number of chips

amp = 1;
ds = 0;
shft = 1;                % How much to shift signal to left to visualize (so that the signal is not centered around 0)
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
nshft = ifftshift((n-r));
tshft = nshft/Fsampl;                                           % Shifted time axis
sqmod = [1 -1].';                                               % Square wave 
sqmod_long = kron(sqmod, ones(Fsampl/Fmodu/2,1));               % Square wave (expanded out based on modulation rate)
scorr_sig = ifft(fft(sqmod_long,Q).*conj(fft(sqmod_long,Q)))/(Fsampl/(Fchip));   
scorr_sig_shift = circshift(scorr_sig, shft); 



%==========================================================================
% Step 2: Define transmission line parameters for PV cable.
%==========================================================================
w=2*pi*f(:);                                         % radial frequency of operation (rad/s)
 
len= len.*3.28084;                                  %length of transmission line cable in foot (some lengths and impedance combinations can result in a distorted signal)
 
% Conductor Parameters
sigma_c = 5.98e7;                                   %conductivity of copper (S/m)
mu_r = .999994;                                     %relative permeability of copper 
mu_c = mu_r*pi*4e-7;                                
mu=pi*4e-7;
dc = 3.15e-3;                                      %Conductor diameter (m)
Dc = 12e-3;                                        %Twin lead wire separation (m)
 
% Insulator Parameters
tanD = 4e-4;                                       % Loss tangent or dissipation factor
eps_r = 2.5;                                       %Relative permittivity
eps_i = eps_r*8.854e-12;
sigma_i = eps_i.*w.*tanD;                          %conductivity of XLPE insulator


%==========================================================================
% Step 3: Calculate RLGC parameters of PV cable (twin lead)
%==========================================================================
for i=1:length(f)
    Rs(i) = sqrt((pi*f(i)*mu_c)/(sigma_c));       %Skin depth
    R(i) = (2*Rs(i))/(pi*dc);                     % Characteristic Resistance
    L(i)= mu/pi *log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    G(i) = pi*sigma_i(i)/log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    C(i) = pi*eps_i/log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    
end


%==========================================================================
% Step 5: Define number of loads, and load impedances
%==========================================================================
num_resis = num_z;                      %number of components separated by T-lines
sstdr_imped = 68;                       %define impedance of SSTDR
zl = imped;                             %impedance of loads


%==========================================================================
% Step 6: define variables for simulation
%==========================================================================
TX = 1;                           %initialize current transfer function
Z_equiv = 0;                      %initialize equivalent resistance
 
%==========================================================================
% Step 7: create simulation using SSP method
%==========================================================================
for i = 1:length(f)                                              %loop over frequency
     Gam(i)=sqrt((R(i)+1i*w(i)*L(i))*(G(i)+1i*w(i)*C(i)));       %complex propagation constant for the TL
    z0(i)= sqrt((R(i)+1i*w(i)*L(i))/(G(i)+1i*w(i)*C(i)));        %characteristic impedance of the transmission line
    z0(1) = .001;                                                %redefine z0 at 0Hz (or else its 0 and it leads to infinities later in the code
    Z_equiv = 0;                                                 %reset equivalent impedance every time you retest system with new freq
    Z_equiv1 = 0;
    for j = num_resis:-1:1                                      %loop over number of elements separated by T-lines
        if j == num_resis                                       %if you are at the beginning of transmission line (furthest from measuring device) do this
            gamma1 = (zl(i,j)-z0(i))/(zl(i,j)+z0(i));           %find reflection coefficient
        else
                                                                %if you are not at the beginning of transmission line do this
        if s_p(j) == 's'                                        % if the system is in series do this
            TX = Z_equiv/(Z_equiv+zl(i,j));                           %calculate voltage transfer function
            Z_equiv = Z_equiv/TX;                                   %find equivalent impedance
        end
        
        if s_p(j) == 'p'                                        % if the system is in parallel do this
            TX = zl(i,j)/(Z_equiv+zl(i,j));                          %calculate current transfer function
            Z_equiv = Z_equiv*TX;                                 %find equivalent impedance
        end
 
        gamma1 = (Z_equiv-z0(i))/(Z_equiv+z0(i));               %find reflection coefficient using equivalent impedance
        end
        Gamma2(i)=gamma1*exp(2*(Gam(i))*(-len(j)*0.3048));      % move reflection coefficient through T-line
        Z_equiv = z0(i)*(Gamma2(i)+1)/(1-Gamma2(i));            %calculate new equivalent impedance
    end
   Gamma3(i) = (Z_equiv-sstdr_imped)/(Z_equiv+sstdr_imped); %reflection coefficient at input of SSTDR (impulse response)
 
end

%==========================================================================
% Step 8: combine input signal data with reflection signature
%==========================================================================
sample_size=.2;                                                 % This variable is used as the sampling size for spline interpolation  %
vop=0.7;                                                        % (Arnold board) approximate VOP for PV cable (Subjected to change as we optimize the algorithm) & corresponding Z0=157ohms%%
Gamma = transpose(Gamma3); 
scorr_y = real(ifft(fft(scorr_sig,[],1).*Gamma,[],1)).*amp(:).';% multiply signal by response in freq domain then take ifft
scorr_y_shift = circshift(scorr_y,shft,1);              
x = 1:length(scorr_y_shift);
    t_ip=1:sample_size:length(scorr_y_shift);                   % This is the new interpolated samples for the x-axis of SSTDR data %%
    ip_raw_data=interp1(x,scorr_y_shift,t_ip,'spline');         % First Spline Interpolation on the Raw data %%
[x_sim,y_sim,~,~,~,~]=slope_calc(tshft(1:Ns),ip_raw_data(1:Ns));% Smooth reflection signature 
dist_m = .5*nshft.*2.99792*1e8*vop/(Fsampl*6);            %resample distance axis
dist_m = dist_m(1:Ns);                                          
ref_sig = y_sim./ocsim_pv(sum(len));                            % Normalize the reflection coefficient
end
