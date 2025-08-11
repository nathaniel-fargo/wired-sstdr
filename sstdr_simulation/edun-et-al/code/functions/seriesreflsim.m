function [X,x] = seriesreflsim(S,Z0,Z1,Zi)
%SERIESREFLSIM  Simulate the response to a series interface impedance
%   [X,x] = SERIESREFLSIM(S,Z0,Z1,Zi) calculates the reflection  
%   coefficient for a series interface impedance Zi in a transmission
%   line with characteristic impedance of Z0 before the interface and a 
%   characteristic impedance of Z1 after it.
%
%   INPUTS:
%       S: A Q-by-1 frequency domain excitation signal at the interface
%      Z0: A Q-by-1 frequency domain characteristic impedance for segment
%          before the interface
%      Z1: A Q-by-1 frequency domain characteristic impedance for segment 
%          after the interface
%      Zi: A Q-by-1 frequency domain interface impedance
%
%   OUTPUTS:
%       X: A Q-by-1 frequency domain reflection signal at the interface
%       x: A Q-by-1 time domain reflection signal at the interface
% 

    % SPECIFY REFLECTION COEFFICIENT
    Gamma = (Z1 + 2*Zi - Z0 )./(Z1 + 2*Zi + Z0);            % Reflection Coefficient
    if isnan(Gamma(1)), Gamma(1) = Gamma(2); end            % Correct if DC gives NaN

    % BUILD SIGNAL
    X = S.*Gamma(:);                                        % Frequency domain signal
    x = real(ifft(X));                                      % Received signal after correlation

end