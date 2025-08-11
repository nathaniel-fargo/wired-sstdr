function [Z0, VOP, prop_constant] = get_RG58TLCharParam(f)
    w = 2*pi*f;

    % RG 58 - using Suhner Coaxial Cable item: 22510015 (Data sheet on ELN)

        % Conductor Parameters
        sigma_c = 5.98E7;       % Conductivity of copper (S/m)
        mu_r_cu = .999994;      % Relative permeability of copper (mu/mu_0)
        mu_0 = pi*4E-7;         % Permeability of free space
        mu_cu = mu_r_cu*mu_0;   % Permeability of copper
        a = (.9e-3)/2;              % Inner Conductor Radius (m) cable from RG58 datasheet on lab archives
        b = (3.6e-3)/2;             % Outer Conductor Radius (m) cable from RG58 datasheet on lab archives

        % Insulator Parameters - Solid PE (Polyethylene)
        mu_r_i = .999994;       % Magnetic permeability of non-ferrous insulator
        tanD = 4E-4;            % Losstangent or dissapation factor
        eps_r = 2.2;           % Relative permittivity
        eps_i = eps_r*8.854E-12;% Permittivity of insulator
        sigma_i = eps_i.*abs(w)*tanD;% Conductivity of insulator
        mu_i = mu_r_i*mu_0;     % Magnetic permeability of the insulator separating conductors a and b

        % RLGC parameters
        Rs = sqrt((pi*(f)*mu_cu)/(sigma_c));          % Skin Depth
        R = Rs./(2*pi)*(1/a + 1/b);                  % Charateristic Resistance      (ohm/m)
        L = mu_i/(2*pi)*log(b/a);                   % Chatrateristic Inductance     (H/m)
        G = (2*pi*sigma_i)/log(b/a);                % Charateristic conductance     (S/m)
        C = (2*pi*eps_i)/log(b/a);                  % Charateristic Capacitance     (F/m)

    %% Charcteristic components using the RLGC parameters
    Z0 = sqrt((R + 1i.*w*L)./(G + 1i.*w*C));             % Characteristic impedance of TL
    % Z0(1) = mean(abs/(Z0(2:length(Z0)/2-1)));               % Removes the infinite impedance at DC and repaces with Z0   
    % NOTE for potentially broken code
    % The followingline was added because I needed to calculate Z0 for
    % positive and negative freqs (5/6/21) 
    % For broken code, try passing the whole positive and negative frequency vector
    Z0(length(Z0)/2) = (mean(abs(Z0((1:length(Z0)/2-1))))+mean(abs(Z0(length(Z0)/2+1:end))))/2;
    prop_constant_temp = sqrt((R + 1i.*w*L).*(G + 1i.*w*C));  % Propogation constant of TL | REAL - att_const (Np/m) | IMAG - phase_const (rad/m)
    prop_constant = [prop_constant_temp(length(w)/2:end); flip(prop_constant_temp(length(w)/2+1:end-1))];
    % prop_constant(length(prop_constant)/2) = (prop_constant(length(prop_constant)-1) + prop_constant(length(prop_constant)+1))/2; % Gets the average of prop_constant around DC to remove NaN
    VOP = w./imag(prop_constant);                        % Velocity of propagation
end

