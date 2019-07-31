function rTrace = dff(trace, t1, t2, Fs, nPlanes)
% trace: 1D vector of F values
% t1: smoothing window (seconds)
% t2: min(trace) window (seconds)
% Fs: image sampling rate
% nPlanes: number of planes used during imaging

planeRate = Fs/nPlanes;
fr1 = ceil(t1 * planeRate);
fr2 = ceil(t2 * planeRate);

% 1. smooth the raw trace
smoothTrace = movmean(trace, [fr1 fr1]);
minS = min(smoothTrace);
if minS < 0
    smoothTrace = smoothTrace + 2*abs(minS);
end

% 2. find min(trace) for a specified window before each timepoint
minTrace = movmin(smoothTrace,[fr2 0]);

% 3. compute dF/F
fTrace = (smoothTrace - minTrace)./ minTrace;

% 4. exponential moving average (alpha = 0.5)
rTrace = filter(.5, [1 .5-1], fTrace);
end
