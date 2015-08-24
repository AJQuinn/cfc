function assess_filters( signal, sample_rate, lo_bounds, hi_bounds, lo_trans, hi_trans )

% Assumes --
% addpath(genpath('/Users/andrew/Projects/HMMPAC/toolbox'));
% rmpath('/Users/andrew/Projects/HMMPAC/toolbox/archive/');
% addpath('/Users/andrew/Software/Matlab/eeglab13_3_2b/functions/sigprocfunc');
% addpath(genpath('/Users/andrew/Software/Matlab/HMM-MAR2/HMM-MAR'));


if nargin < 6 || isempty(lo_trans)
    % Use a constant transition band
    lo_trans = [lo_bounds(1)-5 lo_bounds(2)+5];
end

if nargin < 5 || isempty(hi_trans)
    hi_trans = [hi_bounds(1)-7.5 hi_bounds(2)+7.5];
end

[nchannels,nsamples] = size(signal);

%% Assess spectra
opt.Fs = sample_rate;
opt.fpass = [1 200];

fit = hmmspectramt(signal',nsamples,opt);


%% Design filters

theta_cfg.order = 512;
theta_cfg.sample_rate = sample_rate;
theta_cfg.centre_freq = (lo_bounds(1)+lo_bounds(2))/2;
theta_cfg.pass_width = lo_bounds(2)-lo_bounds(1);
theta_cfg.trans_width = lo_trans(2) - lo_trans(1);
theta_cfg.method = 'twopass';

theta = cfc_filt_fir(signal(1,:),theta_cfg);

gamma_cfg.order = 512;
gamma_cfg.sample_rate = sample_rate;
gamma_cfg.centre_freq = (hi_bounds(1)+hi_bounds(2))/2;
gamma_cfg.pass_width = hi_bounds(2)-hi_bounds(1);
gamma_cfg.trans_width = hi_trans(2) - hi_trans(1);
gamma_cfg.method = 'twopass';

if nchannels == 1
    gamma = cfc_filt_fir(signal(1,:),gamma_cfg);
elseif nchannels == 2
    gamma = cfc_filt_fir(signal(2,:),gamma_cfg);
end

lo_fit = hmmspectramt(theta',nsamples,opt);
hi_fit = hmmspectramt(gamma',nsamples,opt);

figure;
subplot(311)
plot(fit.state.f,log(fit.state.psd)); hold on;
line([hi_bounds(1) hi_bounds(1)],get(gca,'YLim'),'Color',[1 0 0])
line([hi_bounds(2) hi_bounds(2)],get(gca,'YLim'),'Color',[1 0 0])
line([hi_trans(1) hi_trans(1)],get(gca,'YLim'),'Color',[.5 0 0])
line([hi_trans(2) hi_trans(2)],get(gca,'YLim'),'Color',[.5 0 0])

line([lo_bounds(1) lo_bounds(1)],get(gca,'YLim'),'Color',[0 1 0])
line([lo_bounds(2) lo_bounds(2)],get(gca,'YLim'),'Color',[0 1 0])
line([lo_trans(1) lo_trans(1)],get(gca,'YLim'),'Color',[0 .5 0])
line([lo_trans(2) lo_trans(2)],get(gca,'YLim'),'Color',[0 .5 0])
title('Original data')
subplot(312)
plot(lo_fit.state.f,log(lo_fit.state.psd)); hold on;
line([lo_bounds(1) lo_bounds(1)],get(gca,'YLim'),'Color',[0 1 0])
line([lo_bounds(2) lo_bounds(2)],get(gca,'YLim'),'Color',[0 1 0])
line([lo_trans(1) lo_trans(1)],get(gca,'YLim'),'Color',[0 .5 0])
line([lo_trans(2) lo_trans(2)],get(gca,'YLim'),'Color',[0 .5 0])
title('Low Filters')
subplot(313)
plot(hi_fit.state.f,log(hi_fit.state.psd)); hold on;
line([hi_bounds(1) hi_bounds(1)],get(gca,'YLim'),'Color',[1 0 0])
line([hi_bounds(2) hi_bounds(2)],get(gca,'YLim'),'Color',[1 0 0])
line([hi_trans(1) hi_trans(1)],get(gca,'YLim'),'Color',[.5 0 0])
line([hi_trans(2) hi_trans(2)],get(gca,'YLim'),'Color',[.5 0 0])
title('High Filters')
xlabel('Frequency (Hz)');
