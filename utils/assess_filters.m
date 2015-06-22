function [filters] = assess_filters( signal, sample_rate, lo_bounds, hi_bounds )

% Assumes -- 
% addpath(genpath('/Users/andrew/Projects/HMMPAC/toolbox'));
% rmpath('/Users/andrew/Projects/HMMPAC/toolbox/archive/');
% addpath('/Users/andrew/Software/Matlab/eeglab13_3_2b/functions/sigprocfunc');
% addpath(genpath('/Users/andrew/Software/Matlab/HMM-MAR2/HMM-MAR'));

[nsamples,nsignals] = size(signal);


%% Assess spectra
opt.Fs = sample_rate;
opt.fpass = [1 150];

fit = hmmspectramt(signal,nsamples,opt);


%% Design filters
filters = [];

% Low frequency filter
[lo_signal, filters.low_wts] = eegfilt_silent(signal',sample_rate,lo_bounds(1),lo_bounds(2),0,[],0,'fir1');
% High frequency filter
[hi_signal, filters.high_wts] = eegfilt_silent(signal',sample_rate,hi_bounds(1),hi_bounds(2),0,[],0,'fir1');

lo_fit = hmmspectramt(lo_signal',nsamples,opt)
hi_fit = hmmspectramt(hi_signal',nsamples,opt)

figure;
subplot(311)
plot(fit.state.f,fit.state.psd); hold on;
line([hi_bounds(1) hi_bounds(1)],get(gca,'YLim'))
line([hi_bounds(2) hi_bounds(2)],get(gca,'YLim'))
line([lo_bounds(1) lo_bounds(1)],get(gca,'YLim'))
line([lo_bounds(2) lo_bounds(2)],get(gca,'YLim'))
subplot(312)
plot(lo_fit.state.f,lo_fit.state.psd); hold on;
line([lo_bounds(1) lo_bounds(1)],get(gca,'YLim'))
line([lo_bounds(2) lo_bounds(2)],get(gca,'YLim'))
subplot(313)
plot(hi_fit.state.f,hi_fit.state.psd); hold on;
line([hi_bounds(1) hi_bounds(1)],get(gca,'YLim'))
line([hi_bounds(2) hi_bounds(2)],get(gca,'YLim'))
xlabel('Frequency (Hz)');
