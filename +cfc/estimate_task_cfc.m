function cfc_results = estimate_task_cfc( signal, cfg )
%% Estimate cross frequency coupling from a signal
%
% cfc_results = estimate_cfc(signal, cfg)
%
% signal is an array [nchannels x nsamples x nepochs] containing a time-series. If
% one channel is passed in both the modulating and modulated time series will
% be extracted from it. If two channels are passed in the modulating
% time-series will be extracted from the first channels and the modulated
% time-series will be extracted from the second.
%
% cfg is a struct with:
%     sr: sampling rate
%     cfg.lo_freq: frequency for modulating signal (eg 6)
%     cfg.lo_bandwidth: filter bandwidth for modulating signal (eg 2)
%     cfg.hi_freq: frequency for modulated signal (eg 60)
%     cfg.hi_bandwidth: bandwidth for modulated signal
%     cfg.metric: cell array with metrics to compute, defaults to Canolty. eg {'MI','PLV'};
%            options are:
%                MI       - Modulation Index: Canolty et al 2008
%                MI_NORM  - Normalised Modulation Index: Canolty et al 2008
%                ESC      - Envelope Signal Correlation
%                NESC     - Normalised Envelope Signal Correlation
%                GLM      - General Linear Model Method
%                PLV      - Phase Locking Value
%                AEC      - Amplitude Envelope Correlation
%     cfg.window_size: length in seconds for sliding window
%
% options for statistical assessment (optional)
%     cfg.nperms: integer denoting the number of surrogate time series to test
%
% and optionally:
%     true_timecourse: 1d signal indicating where pac exists
%     zero_pad: the number of samples to pad the time series when filtering
%     ret_cfc_signals: return the phase and amplitude signals used in cfc estimation. true|false


%% Housekeeping

if ndims(signal) == 2
    error('Please pass in 3d data!, [nchannels,nsamples,nepochs]');
else
    [nchannels,nsamples,nepochs] = size(signal);
end

if ~isfield(cfg,'ret_cfc_signals')
    cfg.ret_cfc_signals = false;
end

if nchannels > 2
    error('%s channels were passed in, two is the maximum', nchannels)
end
cfg
if ~isfield(cfg,'pad')
    cfg.pad = 100;
end

if ~isfield(cfg,'pretrig')
    cfg.pretrig = 0;
end

if ~isfield(cfg,'true_timecourse')
    cfg.true_timecourse = zeros(size(signal,1),size(signal,2));
end

% Generate time vector
if isfield(cfg,'signal')
    [nsamples,~] = size(signal);
end
time_vect = (0:1/cfg.sr:(nsamples-1) * (1/cfg.sr)) - cfg.pretrig;

% Generate frequency bounds
lo_bounds = [cfg.lo_freq - cfg.lo_bandwidth/2, cfg.lo_freq + cfg.lo_bandwidth/2];
hi_bounds = [cfg.hi_freq - cfg.hi_bandwidth/2, cfg.hi_freq + cfg.hi_bandwidth/2];

% Set sliding window values if requested
if isfield(cfg,'window_size')
    % Convert window and step to samples
    window_size = floor( cfg.window_size / ( 1000 / cfg.sr ) );

    % Find number of windows
    nwindows = floor ( (nsamples - window_size) / window_size);
else
    error('Please define sliding window parameters, window_size')
end

% Permutation defaults
if ~isfield( cfg, 'nperms')
    nperms = 0;
else
    nperms = cfg.nperms;
end

% Create results output
cfc_results = [];
cfc_results.cfg = cfg;
cfc_results.signal = signal;

%% Main body

% Get the observed signals
cfc_signals = cfc.util.basesignals( signal, cfg.sr, ...
    hi_bounds, ...
    lo_bounds, ...
    time_vect, ...
    cfg.true_timecourse );

% Add to output if requested
if cfg.ret_cfc_signals == true
    cfc_results.cfc_signals = cfc_signals;
end

% Stack data for analysis
cfc_signals = cfc.util.stacktrials(cfc_signals,cfg.window_size,'stack');

for win_idx = 1:nwindows

    % Estimate observed CFC
    for met_idx = 1:length(cfg.metrics)
        if strcmp(cfg.metrics{met_idx},'ESC')
            esc = cfc.est.corr(cfc_signals.theta,cfc_signals.gamma_amp);
            cfc_results.esc = esc;
        elseif strcmp(cfg.metrics{met_idx},'NESC')
            nesc = ncfc.est.corr(cos(cfc_signals.theta_phase),cfc_signals.gamma_amp);
            cfc_results.nesc = nesc;
        elseif strcmp(cfg.metrics{met_idx},'AEC')
            aec = cfc.est.corr(cfc_signals.theta_amp,cfc_signals.gamma_amp);
            cfc_results.aec = aec;
        elseif strcmp(cfg.metrics{met_idx},'PLV')
            plv = cfc.est.plv(cfc_signals.theta_phase,cfc_signals.gamma_amp_phase);
            cfc_results.plv = plv;
        elseif strcmp(cfg.metrics{met_idx},'GLM')
            glm = cfc.est.glm(cfc_signals.theta_phase,cfc_signals.gamma_amp);
            cfc_results.glm = glm;
        elseif strcmp(cfg.metrics{met_idx},'MI')
            mi = cfc.est.mi(cfc_signals.theta_phase,cfc_signals.gamma_amp);
            cfc_results.mi = mi;
        elseif strcmp(cfg.metrics{met_idx},'MI_NORM')
            mi_norm = cfc.est.minorm(cfc_signals.theta_phase,cfc_signals.gamma_amp);
            cfc_results.mi_norm = mi_norm;
        else
            fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
        end
    end

end
