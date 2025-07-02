% Frequency sweep
f = logspace(3, 9, 500);  % 1 kHz to 1 GHz
w = 2*pi*f;

% Geometry and material parameters (customize as needed)
mu0 = 4*pi*1e-7;           % Vacuum permeability
eps0 = 8.854e-12;          % Vacuum permittivity
sigma_c = 5.8e7;           % Copper conductivity [S/m]
eps_r = 2.2;               % Relative permittivity of dielectric
tan_delta = 0.002;         % Dielectric loss tangent
d = 0.01;                  % Distance between conductors [m]
r = 1e-3;                  % Conductor radius [m]
l = 1;                     % Line length [m] (optional for scaling)

% Frequency-dependent per-unit-length RLCG
delta = sqrt(2 ./ (w * mu0 * sigma_c));  % Skin depth
R = 1 ./ (2*pi*r*sigma_c*delta);         % Skin-effect resistance
L = mu0 / (2*pi) * log(d/r);             % Approximate inductance

C = 2*pi*eps0*eps_r / log(d/r);          % Capacitance
G = 2*pi*f.*C.*tan_delta;                % Dielectric loss

% Construct Z and Y (per unit length)
ZcPls = R + 1j*w.*L;  % Series impedance
YcPls = G + 1j*w.*C;  % Shunt admittance

% Bundle into struct
LineData.Frequency = f;
LineData.ZcPls = ZcPls;
LineData.YcPls = YcPls;

% Send to workspace for Simulink
assignin('base','LineData',LineData);


figure;
subplot(2,1,1);
semilogx(f, abs(ZcPls), 'LineWidth', 1.5); grid on;
ylabel('|Z| [Ohm/m]');
title('Series Impedance per Unit Length');

subplot(2,1,2);
semilogx(f, abs(YcPls), 'LineWidth', 1.5); grid on;
ylabel('|Y| [S/m]');
xlabel('Frequency [Hz]');
title('Shunt Admittance per Unit Length');