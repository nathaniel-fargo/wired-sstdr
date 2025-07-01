
%% 1 a.  Choose a primitive polynomial for a 10-bit m-sequence
% x^10 + x^3 + 1  → exponent list [10 3 0]
fs              = 1e6;
Ts              = 1/fs;
interpFactor    = 8;

seq     = comm.PNSequence( ...
        'Polynomial',[10 3 0], ...              % ✔ legal primitive
        'SamplesPerFrame',1023, ...
        'InitialConditions',[zeros(1,9) 1]);    % any non-zero seed

mag     = 1;
chips   = 2*mag*seq() - mag;                    % ±1 PN chips (one period)

code_interp     = repelem(chips, interpFactor);

assignin('base','code_seq',chips');

t       = (0:numel(chips)-1)'/fs;           % time vector
codeTS  = timeseries(chips, t);             % TimeSeries for From Workspace

assignin('base','code_output',codeTS);
