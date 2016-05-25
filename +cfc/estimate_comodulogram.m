function cmg = estimate_comodulogram( signal, cfg )
%% Estimate the co-modulogram from a signal
%
% This function can be called in two ways.
%
% cmg = estimate_comodulogram(signal, cfg)
%
% signal is an array [nchannels x nsamples x nrealisations] containing a time-series. If
% one channel is passed in both the modulating and modulated time series will
% be extracted from it. If two channels are passed in the modulating
% time-series will be extracted from the first channels and the modulated
% time-series will be extracted from the second.
%
% cfg is a struct with:
%     sr: sampling rate
%     cfg.lo_bounds: frequency range for modulating signal (eg [5 7])
%     cfg.lo_step: frequency steps for modulating signal (eg 1)
%     cfg.lo_bandwidth: filter bandwidth for modulating signal (eg 2)
%     cfg.hi_bounds: frequency range for modulated signal (eg [50 60])
%     cfg.hi_step: frequency steps for modulated signal (eg 10)
%     cfg.hi_bandwidth: bandwidth for modulated signal, either int or 'adaptive'
%     cfg.metric: cell array with metrics to compute, defaults to Canolty. eg {'MI','PLV'};
%            options are:
%                MI       - Modulation Index: Canolty et al 2008
%                MI_NORM  - Normalised Modulation Index: Canolty et al 2008
%                ESC      - Envelope Signal Correlation
%                NESC     - Normalised Envelope Signal Correlation
%                GLM      - General Linear Model Method
%                PLV      - Phase Locking Value
%                AEC      - Amplitude Envelope Correlation
%
% and optionally:
%     cfg.true_timecourse: 1d signal indicating where pac exists
%     cfg.zero_pad: the number of samples to pad the time series when filtering
%     cfg.theta_interp: bool to indicate that cfc_util_thetawaveform should be used
%     cfg.window_size: length in seconds for sliding window
%     cfg.step: step size between windows in seconds
%     cfg.filter_method: method for applying filter (onepass|twopass|eeglab)

%% Housekeeping
if ndims(signal) == 2
    [nchannels,nsamples] = size(signal);
    ntrials = 1;
else
    [nchannels,nsamples,ntrials] = size(signal);
end

if nchannels > 2
    error('%s channels were passed in, two is the maximum', nchannels)
end

if nargin < 2
    verbose = false;
end

if ~isfield(cfg,'filter_method');
    cfg.filter_method = 'twopass';
end

if ~isfield(cfg,'theta_interp')
    cfg.theta_interp = false;
end

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

if isfield(cfg,'window_size') && isfield(cfg,'step')
    % Convert window and step to samples
    window_size = floor( cfg.window_size / ( 1 / cfg.sr ) );
    step = floor( cfg.step / ( 1 / cfg.sr ) );

    % Find number of windows
    nwindows = floor ( (nsamples - window_size) / step);
else
    nwindows = 1;
end


%% Create frequency ranges
%
% Create a [2 x nsteps] matrix of low frequency bounds to scan. lo_freqs(1,:)
% and lo_freqs(2,:) contain the high and low bounds to filter on based on the
% center frequencies in lo_bounds, lo_steps and lo_bandwidth. High frequency
% bandlimits are created in a [2, n_hi_steps, n_lo_steps] matrix. This is to
% allow for an adaptive bandwidth based on the lowest low frequency (if
% requested).

n_lo_steps = (cfg.lo_bounds(2)-cfg.lo_bounds(1))/cfg.lo_step + 1;
lo_freqs = ones(2,n_lo_steps) .* repmat(cfg.lo_bounds(1):cfg.lo_step:cfg.lo_bounds(2),2,1);
lo_freqs(1,:) = lo_freqs(1,:) - cfg.lo_bandwidth/2;
lo_freqs(2,:) = lo_freqs(2,:) + cfg.lo_bandwidth/2;

n_hi_steps = floor((cfg.hi_bounds(2)-cfg.hi_bounds(1))/cfg.hi_step + 1);
hi_freqs = ones(2,n_hi_steps,n_lo_steps) .* repmat(cfg.hi_bounds(1):cfg.hi_step:cfg.hi_bounds(2),[2,1,n_lo_steps]);

if strcmp(cfg.hi_bandwidth,'adaptive')
    for idx = 1:n_lo_steps
        hi_bandwidth = lo_freqs(1,idx)+2;
        hi_freqs(1,:,idx) = hi_freqs(1,:,idx) - ones(1,n_hi_steps,1)*hi_bandwidth;
        hi_freqs(2,:,idx) = hi_freqs(2,:,idx) + ones(1,n_hi_steps,1)*hi_bandwidth;
    end
else
    size(hi_freqs)
    size(ones(n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2)
     hi_freqs(1,:,:) = hi_freqs(1,:,:) - ones(1,n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2;
     hi_freqs(2,:,:) = hi_freqs(2,:,:) + ones(1,n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2;

end

% Preallocate metrics

esc = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
nesc = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
plv = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
mi = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
mi_norm = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
mi_tort = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
glm = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
glmdv = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);
aec = nan(nwindows,n_lo_steps,n_hi_steps,ntrials);


%%%%%%%%%%%%%%%%%%
%% Begin Main Loop
%%%%%%%%%%%%%%%%%%
msg = '';

for lo_idx = 1:n_lo_steps
    for hi_idx = 1:n_hi_steps

        fprintf(repmat(char(8),1,length(msg)));
        msg = sprintf('Computing low - %.2f:%.2f high - %.2f:%.2f', ...
            lo_freqs(1,lo_idx),lo_freqs(2,lo_idx),hi_freqs(1,hi_idx,lo_idx),hi_freqs(2,hi_idx,lo_idx));
        fprintf(msg);

        if hi_freqs(1,hi_idx,lo_idx) <= lo_freqs(2,lo_idx)
            % Skip estimation if the hi and lo bands overlap
            continue
        end

        %% Create PAC signal -

        if cfg.theta_interp == true
            % Extract phase through interpolation
            lo = mean([lo_freqs(1,lo_idx), lo_freqs(2,lo_idx)]);
        else
            lo = [lo_freqs(1,lo_idx), lo_freqs(2,lo_idx)];
        end

        signals = cfc_util_basesignals(signal, ...
                                   cfg.sr, ...
                                   [hi_freqs(1,hi_idx,lo_idx) hi_freqs(2,hi_idx,lo_idx)], ...
                                   lo,...
                                   time_vect, ...
                                   cfg.true_timecourse,[],[],[],cfg.filter_method);


        if isfield(cfg,'window_size')
            [ signals ] = cfc_util_swsignals(signals,window_size,step);
        end

       % Estimate CFC
        for met_idx = 1:length(cfg.metrics)
            if strcmp(cfg.metrics{met_idx},'ESC')
                esc(:,lo_idx,hi_idx,:) = cfc.est.corr(signals.theta,signals.gamma_amp);
            elseif strcmp(cfg.metrics{met_idx},'NESC')
                nesc(:,lo_idx,hi_idx,:) = ncfc.est.corr(cos(signals.theta_phase),signals.gamma_amp);
            elseif strcmp(cfg.metrics{met_idx},'AEC')
                aec(:,lo_idx,hi_idx,:) = cfc.est.corr(signals.theta_amp,signals.gamma_amp);
            elseif strcmp(cfg.metrics{met_idx},'PLV')
                plv(:,lo_idx,hi_idx,:) = cfc.est.plv(signals.theta_phase,signals.gamma_amp_phase);
            elseif strcmp(cfg.metrics{met_idx},'GLM')
                tmp = cfc.est.glm(signals.theta_phase,signals.gamma_amp);
                glm(:,lo_idx,hi_idx,:) = tmp.fnorm;
            elseif strcmp(cfg.metrics{met_idx},'GLMDV')
                tmp = cfc.est.glmdv(signals.theta_phase,signals.gamma_amp);
                glmdv(:,lo_idx,hi_idx,:) = tmp.fnorm;
            elseif strcmp(cfg.metrics{met_idx},'MI')
                mi(:,lo_idx,hi_idx,:) = cfc.est.mi(signals.theta_phase,signals.gamma_amp);
            elseif strcmp(cfg.metrics{met_idx},'MI_NORM')
                mi_norm(:,lo_idx,hi_idx,:) = cfc.est.minorm(signals.theta_phase,signals.gamma_amp);
            elseif strcmp(cfg.metrics{met_idx},'MI_TORT')
                 mi_tort(:,lo_idx,hi_idx,:) = cfc.est.mitort(signals.theta_phase,signals.gamma_amp);
            else
                 fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
            end
        end
    end
end

fprintf('\n');
%% Create output struct

cmg = struct('lo_freqs',    lo_freqs,...
             'hi_freqs',    hi_freqs,...
             'signal',      signals.signal,...
             'time_vect',   signals.time_vect);

% add metrics
for met_idx = 1:length(cfg.metrics)
    if strcmp(cfg.metrics{met_idx},'ESC')
        cmg.esc = esc;
    elseif strcmp(cfg.metrics{met_idx},'NESC')
        cmg.nesc = nesc;
    elseif strcmp(cfg.metrics{met_idx},'AEC')
        cmg.aec = aec;
    elseif strcmp(cfg.metrics{met_idx},'PLV')
        cmg.plv = plv;
    elseif strcmp(cfg.metrics{met_idx},'GLM')
        cmg.glm = glm;
    elseif strcmp(cfg.metrics{met_idx},'GLMDV')
        cmg.glmdv = glmdv;
    elseif strcmp(cfg.metrics{met_idx},'MI')
        cmg.mi = mi;
    elseif strcmp(cfg.metrics{met_idx},'MI_NORM')
        cmg.mi_norm = mi_norm;
    elseif strcmp(cfg.metrics{met_idx},'MI_TORT')
        cmg.mi_tort = mi_tort;
    end
end

end
