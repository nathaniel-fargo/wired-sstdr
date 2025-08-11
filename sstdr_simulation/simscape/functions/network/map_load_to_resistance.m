function [R_value, element_type] = map_load_to_resistance(load_value, varargin)
%MAP_LOAD_TO_RESISTANCE Convert reflection coefficient to resistance for Simscape faults
% This helper converts a *reflection coefficient* (Γ) in the range [-1,1]
% into an equivalent lumped resistance that can be inserted either *in series*
% with the transmission line (for Γ > 0, "open"-type fault) or *in shunt*
% to ground (for Γ < 0, "short"-type fault).
%
% Reflection coefficient definitions used (assuming purely resistive fault):
%   • Series element  :  Γ =  R_s  / (2·Z₀ + R_s)
%                        R_s = 2·Z₀·Γ / (1-Γ)
%   • Shunt element   :  Γ = − 2·Z₀ / (R_p + 2·Z₀)
%                        R_p = −2·Z₀·(1+Γ)/Γ
%
% where Z₀ is the characteristic impedance (default 50 Ω).
%
% Special cases:
%   Γ ≈ 0    → no extra element (transmission line continues)
%   Γ →  1   → R_s → ∞   (capped at 10 MΩ for numerical stability)
%   Γ → −1   → R_p → 0   (floored at 0.1 Ω for numerical stability)
%
% Usage examples:
%   [R, type] = map_load_to_resistance( 0.5)        % R ≈ 2·Z0  series fault
%   [R, type] = map_load_to_resistance(-0.5)        % R ≈ 2·Z0  shunt  fault
%   [R, ~   ] = map_load_to_resistance( 0  , 'Z0', 75) % Pass-through (no fault)
%
% Input:
%   load_value – Scalar Γ in range [-1, 1]. Positive ⇒ series, negative ⇒ shunt
%   Name-value pairs:
%     'Z0' – Characteristic impedance (default 50 Ω)
%
% Output:
%   R_value      – Equivalent resistance (Ohms). NaN if no element (Γ≈0)
%   element_type – 'none', 'series', or 'shunt'

% Parse optional arguments
p = inputParser;
addParameter(p, 'Z0', 50, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

% Validate input
if load_value < -1 || load_value > 1
    error('map_load_to_resistance:OutOfRange', ...
        'Load value must be in range [-1, 1], got %.3f', load_value);
end

% Determine element type and resistance value
if abs(load_value) < 1e-6
    % Transmission line segment (no extra element)
    R_value = NaN;
    element_type = 'none';
    
elseif load_value > 0
    % Series fault (positive load = more open)
    element_type = 'series';
    Z0 = p.Results.Z0;
    % Avoid division-by-zero when Γ → 1
    if abs(1 - load_value) < eps
        R_value = inf;
    else
        % R_s = 2·Z₀·Γ / (1-Γ)
        R_value = 2 * Z0 * load_value / (1 - load_value);
    end
    
else
    % Shunt fault (negative load = more short)
    element_type = 'shunt';
    Z0 = p.Results.Z0;
    % Avoid division-by-zero when Γ → −1
    if abs(load_value + 1) < eps
        R_value = 0;
    else
        % R_p = −2·Z₀·(1+Γ)/Γ   (note Γ<0)
        R_value = -2 * Z0 * (1 + load_value) / load_value;
    end
end

% Numerical caps for simulator stability
if ~isnan(R_value)
    if R_value < 0.1          % avoid perfect short (R→0)
        R_value = 0.1;
    elseif R_value > 1e7      % cap very large opens to 10 MΩ
        R_value = 1e7;
    end
end

end 