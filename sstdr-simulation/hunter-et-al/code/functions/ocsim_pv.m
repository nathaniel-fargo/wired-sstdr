function [oc_sim] = ocsim_pv(len)
% Dr. Harleys sstdr signal simulation simulation

Fsampl = 96e6;          % Sampling rate
Fmodu  = Fsampl/4;      % Modulation rate
Fchip  = Fsampl/4;      % Chip rate
PNL    = 10240;         % Number of chips
shft = 10;           % How much to shift signal to left to visualize (so that the signal is not centered around 0)
                     %   THIS IS FOR PLOTTING PURPOSES ONLY
Ns = 120;            % Number of samples to plot
amp = 1;
ds = 0;
shft = 10;           % How much to shift signal to left to visualize (so that the signal is not centered around 0)
                     %   THIS IS FOR PLOTTING PURPOSES ONLY
Ns = 1234;            % Number of samples to plot
Q = PNL*Fsampl/(Fchip);         % Length of simulated signal
n = 1:Q;                        % Sample axis
t = n/Fsampl;                   % Time axis
r = floor(Q/2)+1;               % Center sample in frequency
f = ifftshift((n-r)/Q)*Fsampl;  % Frequency axis
nshft = circshift(ifftshift((n-r)), shft);      % Shifted sample axis
tshft = nshft/Fsampl;                           % Shifted time axis
sqmod = [1 -1].';                                               % Square wave 
sqmod_long = kron(sqmod, ones(Fsampl/Fmodu/2,1));               % Square wave (expanded out based on modulation rate)
scorr_sig = ifft(fft(sqmod_long,Q).*conj(fft(sqmod_long,Q)))/(Fsampl/(Fchip));   
scorr_sig_shift = circshift(scorr_sig, shft); 




% vop = .695; % vop of rg58





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system simulation





w=2*pi*f(:);  % radial frequency of operation (rad/s)


% Conductor Parameters
sigma_c = 5.98e7; %conductivity of copper (S/m)
mu_r = .999994; %relative permiability of copper ()
mu_c = mu_r*pi*4e-7;
mu=pi*4e-7;
dc = 3.15e-3; %Conductor diameter (m)
Dc = 12e-3; %Twin lead wire seperation (m)

% Insulator Parameters
tanD = 4e-4; % Losstangent or dissapation factor
eps_r = 2.5; %Relative permativity
eps_i = eps_r*8.854e-12;
sigma_i = eps_i.*w.*tanD; %conductivity of XLPE insulator

for i=1:length(f)
    Rs(i) = sqrt((pi*f(i)*mu_c)/(sigma_c)); %Skin depth
    R(i) = (2*Rs(i))/(pi*dc); % Charateristic Resistance
    L(i)= mu/pi *log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    G(i) = pi*sigma_i(i)/log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    C(i) = pi*eps_i/log((Dc/dc)+sqrt((Dc/dc)^2 -1));
    
end
TX = 1; %initilize current transfer function
Z_equiv = 0; %initilize equivilent resistance
num_resis = 1 ;%number of components seperated by T-lines
sstdr_imped = 68;


for i = 1:length(f) %loop over frequency
     zl(i) = 1e17; %8.81e-10
     Gam(i)=sqrt((R(i)+1i*w(i)*L(i))*(G(i)+1i*w(i)*C(i)));       %complex propagation constant for the TL
       %zl = 100; %load impedance for resistor
    z0(i)= sqrt((R(i)+1i*w(i)*L(i))/(G(i)+1i*w(i)*C(i)));  %characteristic impedance of the transmission line
    z0(1) = 1.5310e+02 - 1.1381e+01i;
    Z_equiv = 0; %reset equivilent impedance every time you retest system with new freq
    Z_equiv1 = 0;
    for j = num_resis:-1:1 %loop over number of elements seperated by T-lines
        if j == num_resis %if you are at the begining of transmission line (furthest from measuring device) do this
            gamma1 = (zl(i)-z0(i))/(zl(i)+z0(i)); %find reflection coeficient
        else                %if you are not at the begining of transmission line do this
        TX = Z_equiv/(Z_equiv+zl(i)); %calculate voltage transfer function
        Z_equiv = Z_equiv/TX;       %find equivilent impedance
        

%         TX = zl(i)/(Z_equiv+zl(i)); %calculate current transfer function
%         Z_equiv = Z_equiv*TX;       %find equivilent impedance
%         z_tracker2(i) = Z_equiv;
        z_equiv = 1;
        gamma1 = (Z_equiv-z0(i))/(Z_equiv+z0(i)); %find reflection coeficient using equivilent impedance
        end
        Gamma2(i)=gamma1*exp(2*(Gam(i))*(-len(j)*0.3048)); % move reflection coeficient through T-line
        Z_equiv = z0(i)*(Gamma2(i)+1)/(1-Gamma2(i)); %calculate new equivilent impedance
    end
     Gamma3(i) = (Z_equiv-68)/(Z_equiv+68);
% induct_sstdr = 6.1576e-06;
% Gamma3(i) = (Z_equiv-(sstdr_imped+((1i*induct_sstdr*w(i)))))/(Z_equiv+(sstdr_imped+((1i*induct_sstdr*w(i)))));
% induct_sstdr = 10.51e-10;%2.51e-10
% w(1) = w(2);
% Gamma3(i) = (Z_equiv-(sstdr_imped+(1/(1i*induct_sstdr*w(i)))))/(Z_equiv+(sstdr_imped+(1/(1i*induct_sstdr*w(i)))));
end

% plot(f,abs(z_tracker2))
% hold on

sample_size=.2;  % This variable is used as the sampling size for spline interpolation  %

count = 0;
    vop=0.7; %% (Arnold board) approximate VOP for PV cable (Subjected to change as we optimize the algorithm) & corresponding Z0=157ohms%%
    tshft = resample(tshft,5,1);
%         tshft = tshft * 5;
 dist = tshft*3.28084*2.99792*1e8*vop*.5; %convert x axis to distance
for k = 1:Ns
    if dist(k) <0
        count = count + 1;
    end
end

Gamma = transpose(Gamma3); %convert to row vector B/c Dr. Harleys code makes a row vector
scorr_y = real(ifft(fft(scorr_sig,[],1).*Gamma,[],1)).*amp(:).';  % multiply signal by respone in freq domain then take ifft
scorr_y_shift = circshift(scorr_y,shft,1);              % FOR PLOTTING PURPOSES ONLY
x = 1:length(scorr_y_shift);
    t_ip=1:sample_size:length(scorr_y_shift);%% This is the new interpolated samples for the x-axis of SSTDR data %%
    ip_raw_data=interp1(x,scorr_y_shift,t_ip,'spline'); %% First Spline Interpolation on the Raw data %%
[x_sim,y_sim,~,~,~,~]=slope_calc(tshft(count:Ns+count),ip_raw_data(count:Ns+count));
oc_sim = max(y_sim);



end