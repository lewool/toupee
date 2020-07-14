function gw = fsGausswin(stdev, a, Fs)
% function gw = myGaussWin(stdev, Fs)
% A gaussian window with specified stdev in units relative to the sampling
% frequency, and normalized in amplitude
stdevSamps = round(stdev*Fs);
gw = wheel.gausswin(stdevSamps, a);
gw = gw./sum(gw);

