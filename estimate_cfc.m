function cfc_results = estimate_cfc( signal, cfg )
%% Estimate cross frequency coupling from a signal
%
% cfc_results = estimate_cfc(signal, cfg)
%
% signal is a 1d or 2d array [nchannels x nsamples] containing a time-series. If
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
%                MI   - Modulation Index: Canolty et al 2008
%                ESC  - Envelope Signal Correlation
%                NESC - Normalised Envelope Signal Correlation
%                GLM  - General Linear Model Method
%                PLV  - Phase Locking Value
%                AEC  - Amplitude Envelope Correlation
%
% options for sliding windows (optional)
%     cfg.window_size: scalar indicating sliding window length in ms
%     cfg.window_step: scalar indicating sliding window step size in ms
%
% options for statistical assessment (optional)
%     cfg.nperms: integer denoting the number of surrogate time series to test
%
% and optionally:
%     true_timecourse: 1d signal indicating where pac exists
%     zero_pad: the number of samples to pad the time series when filtering
%     window_size: length in seconds for sliding window
%     step: step size between windows in seconds


%% Housekeeping
[nchannels,nsamples] = size(signal);
if nchannels > 2
    error('%s channels were passed in, two is the maximum', nchannels)
end

if ~isfield(cfg,'pad')
    cfg.pad = 100;
end

if ~isfield(cfg,'true_timecourse')
    cfg.true_timecourse = zeros(size(signal,1),size(signal,2));
end

% Generate time vector
if isfield(cfg,'signal')
    [nsamples,~] = size(signal);
end
time_vect = (0:1/cfg.sr:(nsamples-1) * (1/cfg.sr));

% Generate frequency bounds
lo_bounds = [cfg.lo_freq - cfg.lo_bandwidth, cfg.lo_freq + cfg.lo_bandwidth];
hi_bounds = [cfg.hi_freq - cfg.hi_bandwidth, cfg.hi_freq + cfg.hi_bandwidth];

% Set sliding window values if requested
if isfield(cfg,'window_size') && isfield(cfg,'window_step')
    % Convert window and step to samples
    window_size = floor( cfg.window_size / ( 1000 / cfg.sr ) );
    window_step = floor( cfg.window_step / ( 1000 / cfg.sr ) );

    % Find number of windows
    nwindows = floor ( (nsamples - window_size) / window_step);
else
    nwindows = 1;
end

% Permutation defaults
if ~isfield( cfg, 'nperms')
    nperms = 0;
else
    nperms = cfg.nperms;
end

%% Create results output
cfc_results = [];
cfc_results.cfg = cfg;
cfc_results.signal = signal;

%% Main body

% Get the observed signals
cfc_signals = make_pac_signals( signal, cfg.sr, ...
    hi_bounds, ...
    lo_bounds, ...
    time_vect, ...
    cfg.true_timecourse );

if nwindows > 1
    cfc_signals  = make_sw_data(cfc_signals,window_size,window_step);
end

cfc_results.time_vect = cfc_signals.time_vect;

% Estimate observed CFC
for met_idx = 1:length(cfg.metrics)
    if strcmp(cfg.metrics{met_idx},'ESC')
        esc = esc_estimator(cfc_signals.theta,cfc_signals.gamma_amp);
        cfc_results.esc = esc;
    elseif strcmp(cfg.metrics{met_idx},'NESC')
        nesc = nesc_estimator(cfc_signals.theta_phase,cfc_signals.gamma_amp);
        cfc_results.nesc = nesc;
    elseif strcmp(cfg.metrics{met_idx},'AEC')
        aec = aec_estimator(cfc_signals.theta_amp,cfc_signals.gamma_amp);
        cfc_results.aec = aec;
    elseif strcmp(cfg.metrics{met_idx},'PLV')
        plv = plv_estimator(cfc_signals.theta_phase,cfc_signals.gamma_amp_phase);
        cfc_results.plv = plv;
    elseif strcmp(cfg.metrics{met_idx},'GLM')
        glm = glm_estimator(cfc_signals.theta_phase,cfc_signals.gamma_amp);
        cfc_results.glm = glm;
    elseif strcmp(cfg.metrics{met_idx},'MI')
        mi = mi_estimator(cfc_signals.theta_phase,cfc_signals.gamma_amp);
        cfc_results.mi = mi;
    else
        fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
    end
end



% Get the surrogate signals if requested
if nperms > 0
    % We only need the surrogates for the modulating time-course
    surrogates = generate_phase_surrogates(signal(1,:),nperms);

    cfc_results.esc_null = nan(nperms,1);
    ncfc_results.esc_null = nan(nperms,1);
    cfc_results.aec_null = nan(nperms,1);
    cfc_results.mi_null = nan(nperms,1);
    cfc_results.glm_null = nan(nperms,1);
    cfc_results.plv_null = nan(nperms,1);

    % I need to compute the maximum statistic across sliding windows for each
    % permutatation. It will therefore be simpler to compute every window for each
    % permutation and add the max across windows to the null

    % Estimate surrogate CFC
    msg = ''; tic; % Setup for update message
    for idx = 1:nperms

        % print update message
        fprintf(repmat(char(8),1,length(msg)));
        msg = sprintf('Computing permutation: %d/%d, elapsed time: %d seconds', ...
            idx,nperms,fix(toc));
        fprintf(msg);

        if nchannels == 2
            surr_signal = cat(1,surrogates(idx,:),signal.signal(2,:));
        else
            surr_signal = surrogates(idx,:);
        end

        surr_signals = make_pac_signals( surr_signal, cfg.sr, ...
            hi_bounds, ...
            lo_bounds, ...
            time_vect, ...
            cfg.true_timecourse );

        if nwindows > 1
            surr_signals  = make_sw_data(surr_signals,window_size,window_step);
        end

        % Estimate surrogate CFC
        for met_idx = 1:length(cfg.metrics)
            if strcmp(cfg.metrics{met_idx},'ESC')
                tmp = esc_estimator(surr_signals.theta,surr_signals.gamma_amp);
                cfc_results.esc_null(idx) = max(tmp);
            elseif strcmp(cfg.metrics{met_idx},'NESC')
                tmp = nesc_estimator(surr_signals.theta_phase,surr_signals.gamma_amp);
                cfc_results.nesc_null(idx) = max(tmp);
            elseif strcmp(cfg.metrics{met_idx},'AEC')
                tmp = aec_estimator(surr_signals.theta_amp,surr_signals.gamma_amp);
                cfc_results.aec_null(idx) = max(tmp);
            elseif strcmp(cfg.metrics{met_idx},'PLV')
                tmp = plv_estimator(surr_signals.theta_phase,surr_signals.gamma_amp_phase);
                cfc_results.plv_null(idx) = max(tmp);
            elseif strcmp(cfg.metrics{met_idx},'GLM')
                tmp = glm_estimator(surr_signals.theta_phase,surr_signals.gamma_amp);
                cfc_results.glm_null(idx) = max(tmp);
            elseif strcmp(cfg.metrics{met_idx},'MI')
                tmp = mi_estimator(surr_signals.theta_phase,surr_signals.gamma_amp);
                cfc_results.mi_null(idx) = max(tmp);
            else
                fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
            end
        end
    end

    % Get critical values
    for met_idx = 1:length(cfg.metrics)
        if strcmp(cfg.metrics{met_idx},'ESC')
            cfc_results.esc_z = ( cfc_results.esc - mean(cfc_results.esc_null) ) / std(cfc_results.esc_null);
            cfc_results.esc_thresh(1) = prctile(cfc_results.esc_null,95);
            cfc_results.esc_thresh(2) = prctile(cfc_results.esc_null,99);
            cfc_results.esc_thresh(3) = prctile(cfc_results.esc_null,99.9);
        elseif strcmp(cfg.metrics{met_idx},'NESC')
            cfc_results.nesc_z = ( cfc_results.nesc - mean(cfc_results.nesc_null) ) / std(cfc_results.nesc_null);
            cfc_results.nesc_thresh(1) = prctile(ncfc_results.esc_null,95);
            cfc_results.nesc_thresh(2) = prctile(ncfc_results.esc_null,99);
            cfc_results.nesc_thresh(3) = prctile(ncfc_results.esc_null,99.9);
        elseif strcmp(cfg.metrics{met_idx},'AEC')
            cfc_results.aec_z = ( cfc_results.aec - mean(cfc_results.aec_null) ) / std(cfc_results.aec_null);
            cfc_results.aec_thresh(1) = prctile(cfc_results.aec_null,95);
            cfc_results.aec_thresh(2) = prctile(cfc_results.aec_null,99);
            cfc_results.aec_thresh(3) = prctile(cfc_results.aec_null,99.9);
        elseif strcmp(cfg.metrics{met_idx},'PLV')
            cfc_results.plv_z = ( cfc_results.plv - mean(cfc_results.plv_null) ) / std(cfc_results.plv_null);
            cfc_results.plv_thresh(1) = prctile(cfc_results.plv_null,95);
            cfc_results.plv_thresh(2) = prctile(cfc_results.plv_null,99);
            cfc_results.plv_thresh(3) = prctile(cfc_results.plv_null,99.9);
        elseif strcmp(cfg.metrics{met_idx},'GLM')
            cfc_results.glm_z = ( cfc_results.glm - mean(cfc_results.glm_null) ) / std(cfc_results.glm_null);
            cfc_results.glm_thresh(1) = prctile(cfc_results.glm_null,95);
            cfc_results.glm_thresh(2) = prctile(cfc_results.glm_null,99);
            cfc_results.glm_thresh(3) = prctile(cfc_results.glm_null,99.9);
        elseif strcmp(cfg.metrics{met_idx},'MI')
            cfc_results.mi_z = ( cfc_results.mi - mean(cfc_results.mi_null) ) / std(cfc_results.mi_null);
            cfc_results.mi_thresh(1) = prctile(cfc_results.mi_null,95);
            cfc_results.mi_thresh(2) = prctile(cfc_results.mi_null,99);
            cfc_results.mi_thresh(3) = prctile(cfc_results.mi_null,99.9);
        else
            fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
        end
    end
end
fprintf('\n');

end % function

