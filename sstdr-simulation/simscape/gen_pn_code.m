%% 1 a.  Choose a primitive polynomial for a 10-bit m-sequence
% x^10 + x^3 + 1  → exponent list [10 3 0]
fs  = 25e6;                       % 25 MHz sample rate  (Ts = 40 ns)

seq = comm.PNSequence( ...
        'Polynomial',[10 3 0], ...       % ✔ legal primitive
        'SamplesPerFrame',1023, ...
        'InitialConditions',[zeros(1,9) 1]);    % any non-zero seed

chips   = 2*seq() - 1;            % ±1 PN chips (one period)

assignin('base','chips',chips)

%% 1 b.  Pad with zeros so the drive idles
tailLen = 2000;                   % 2 000 zeros  (≈80 µs at 25 MHz)
paddedChips   = [chips ; zeros(tailLen,1)];

t       = (0:numel(paddedChips)-1)'/fs;           % time vector
codeTS  = timeseries(paddedChips, t);             % TimeSeries for From Workspace

assignin('base','code',codeTS);
