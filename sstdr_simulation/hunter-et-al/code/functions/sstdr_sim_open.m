function [oc_sim] = sstdr_sim_open(len)
% Dr. Harleys sstdr signal simulation simulation

% Fsampl = 192e6;          % Sampling rate 48mhz
Fsampl = 4*24e6;
Fmodu  = Fsampl/4;      % Modulation rate
Fchip  = Fsampl/4;      % Chip rate
PNL    = 10240;         % Number of chips
shft = 10;           % How much to shift signal to left to visualize (so that the signal is not centered around 0)
                     %   THIS IS FOR PLOTTING PURPOSES ONLY

amp = 1;
ds = 0;
shft = 10;           % How much to shift signal to left to visualize (so that the signal is not centered around 0)
                     %   THIS IS FOR PLOTTING PURPOSES ONLY
Ns = 200;            % Number of samples to plot
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




vop = .66; % vop of rg58






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system simulation

% the reason i used rs a and b was to find the resistance in the lumped
% element model. i initally used a value off of the ECE 3300 website but it
% gave me bad values and when i calculated it my self it worked

rs = 1/(5.8*10^7); % impedance of copper
a=.445*10^-3; % inner radius of RG58
b=1.765*10^-3; %outer radius of RG58
%lumped element piecwies values for RG58
r = (rs/(2*pi))*((1/a)+(1/b)); %resistance of RG58 per m
l = 2.7557e-7; %inductance of RG58 per m
g = 4.5602e-4; % admitance of RG58 per m
c = 9.1164e-11; %capacitance of rg58 per m
% vop = .695; % vop of rg58
x = r+(1i*2*pi*f(:)*l);
y = g+(1i*2*pi*f(:)*c);
gam = sqrt(x.*y); % calculating gamma 
z0 = 52; %54.976-(1i*.0368); %Z0 for rg58
% zl = [100,110,110];

TX = 1; %initilize current transfer function
Z_equiv = 0; %initilize equivilent resistance
num_resis = 1;%number of components seperated by T-lines
% sstdr_imped = sstdr_impedance(Q,f);
%I used alpha and beta initially to debug a problem i was having. i could
%switch it back to gamma and there should be no differnece
% myfiles = {'RC_45.8ft.lws'};
 beta = 2*pi*f(:)*sqrt(l*c); % betta of T-line
alpha =  .5*((r/z0) + (g*z0)); %alpha of T-line 6.6274e-04;



% zl = impedance(Q,f); %rlc meter
% zl = impedance_cap(Q,f);

for i = 1:length(f) %loop over frequency
%      w(i) = 2*pi*f(i);
%      zl(i) = ZCell(0.4, 300, w(i), 161.29, 0.02, 1, 0.0248, 0, 10000, 10^19, 10^16, 10^10, 11.7, 95.2217, 410.7919, 5e-8, 5e-5);
%         f(1) = 20;
%        zl = [-1i/(f(i)*2*pi*100e-12),-1i/(f(i)*2*pi*100e-12),-1i/(f(i)*2*pi*100e-12)]; %load impedance for resistor
%          
%        zl_test(i) = zl;
       zl(i)=10000;
        sstdr_imped(i) = 68;
    Z_equiv = 0; %reset equivilent impedance every time you retest system with new freq
    for j = num_resis:-1:1 %loop over number of elements seperated by T-lines
        if j == num_resis %if you are at the begining of transmission line (furthest from measuring device) do this
            gamma1 = (zl(i)-z0)/(zl(i)+z0); %find reflection coeficient
        else                %if you are not at the begining of transmission line do this
        TX = zl(i)/(Z_equiv+zl(i)); %calculate current transfer function
        Z_equiv = TX*Z_equiv;       %find equivilent impedance
        gamma1 = (Z_equiv-z0)/(Z_equiv+z0); %find reflection coeficient using equivilent impedance
        end
        Gamma2(i)=gamma1*exp(2*(alpha+(1i*beta(i)))*(-len(j)*0.3048)); % move reflection coeficient through T-line
        Z_equiv = z0*(Gamma2(i)+1)/(1-Gamma2(i)); %calculate new equivilent impedance
    end
    question(i) = Z_equiv;
    Gamma3(i) = (Z_equiv-sstdr_imped(i))/(Z_equiv+sstdr_imped(i));
    
end





sample_size=.2;  % This variable is used as the sampling size for spline interpolation  %
% if m ==1
count = 0;
%     vop=0.7; %% (Arnold board) approximate VOP for PV cable (Subjected to change as we optimize the algorithm) & corresponding Z0=157ohms%%
    tshft = resample(tshft,5,1);
 dist = tshft*3.28084*2.99792*1e8*vop*.5; %convert x axis to distance
for k = 1:Ns
    if dist(k) <0
        count = count + 1;
    end
end
% end

Gamma = transpose(Gamma3); %convert to row vector B/c Dr. Harleys code makes a row vector
scorr_y = real(ifft(fft(scorr_sig,[],1).*Gamma,[],1)).*amp(:).';  % multiply signal by respone in freq domain then take ifft
scorr_y_shift = circshift(scorr_y,shft,1);              % FOR PLOTTING PURPOSES ONLY
x = 1:length(scorr_y_shift);
    t_ip=1:sample_size:length(scorr_y_shift);%% This is the new interpolated samples for the x-axis of SSTDR data %%
    ip_raw_data=interp1(x,scorr_y_shift,t_ip,'spline'); %% First Spline Interpolation on the Raw data %%
[x_sim,y_sim,~,~,~,~]=slope_calc(tshft(count:Ns+count),ip_raw_data(count:Ns+count));
y_sim = y_sim(1:Ns+1);
dist_final = dist(count:Ns + count);
oc_sim = max(y_sim);
end