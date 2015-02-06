function out = pacman_sw(obj)
%
% Phase amplitude coupling measures are numerous
%
% Must be passed a struct with:
%     sr: sampling rate
%     window_size: length in seconds for sliding window
%     step: step size between windows in seconds
% and EITHER:
%     signal: 1d broadband signal to be analysed
%     lo_bounds: frequency range for modulating signal (eg [5 7])
%     hi_bounds: frequency range for modulated signal (eg [50 60])
% OR
%     theta: 1d modulating signal
%     gamma: 1d modulated signal
%     lo_bounds: frequecy range for modulating signal
%
% and optionally:
%     true_timecourse: 1d signal indicating where pac exists
%     zero_pad: the number of samples to pad the time series when filtering
% Function computing a range of PAC related measures
%
%     Modulation Index - Canolty et al 2006
%     Phase Locking Value - Lachaux 1999
%     GLM PAC - Penny et al 2008
%     Envelope-to-signal correlation - Bruns & Eckhorn 2004
%     Normalised envelope-to-signal correlation  -Penny et al 2008

%% Housekeeping

% Make sure we haven't been given contradicting arguments
if isfield(obj,'signal') && isfield(obj,'theta')
    error('Please define either a signal and frequency bands or separate modulating and modulated signals');
end

if ~isfield(obj,'pad')
    obj.pad = 100;
end

% Generate time vector
if isfield(obj,'signal')
    [~,nsamples] = size(obj.signal);
else
    [~,nsamples] = size(obj.theta);
end
time_vect = 0:1/obj.sr:(nsamples-1) * (1/obj.sr);

% If we are passed one signal with frequency bands, perform filtering to
% create the modulating and modulated signal
if isfield(obj,'signal')    
    % Compute frequencies of interest
    pad_signal = [zeros(1,obj.pad), obj.signal, zeros(1,obj.pad)];
    theta = eegfilt_silent(pad_signal,obj.sr,obj.lo_bounds(1),obj.lo_bounds(2),0,[],0,'fir1');
    theta = theta(obj.pad:end-obj.pad-1);
    gamma = eegfilt_silent(pad_signal,obj.sr,obj.hi_bounds(1),obj.hi_bounds(2),0,[],0,'fir1');
    gamma = gamma(obj.pad:end-obj.pad-1);
else
    theta = obj.theta;
    gamma = obj.gamma;
end

%% Compute signals

% Compute amplitude time series
theta_amp = abs(hilbert(theta));
gamma_amp = abs(hilbert(gamma));

% Compute phase time series
theta_phase = angle(hilbert(theta));
% gamma_phase = angle(hilbert(gamma)); %% Unused

% Compute phase of gamma amplitude
gamma_amp_phase = angle(hilbert(gamma_amp));

% Compute theta-filtered gamma amplitude
%if isfield(obj,'signal')
    pad_gamma = [zeros(1,obj.pad), gamma_amp, zeros(1,obj.pad)];
    pad_gamma = eegfilt_silent(pad_gamma,obj.sr,obj.lo_bounds(1),obj.lo_bounds(2),0,[],0,'fir1');
    gamma_amp_theta = pad_gamma(obj.pad:end-obj.pad-1);
    %else
    % We cannot compute this for the two signal approach as we don't know
    % the exact epected theta bounds a priori
%    gamma_amp_theta = nan;
%end

signals = struct('theta', theta, ...
                 'gamma', gamma, ...
                 'gamma_amp',gamma_amp, ...
                 'theta_amp',theta_amp, ...
                 'theta_phase', theta_phase, ...
                 'gamma_amp_phase', gamma_amp_phase, ...
                 'gamma_amp_theta', gamma_amp_theta, ...
                 'time_vect', time_vect);

%% Set up sliding window

% Convert window and step to samples
window_size = floor( obj.window_size / ( 1 / obj.sr ) );
step = floor( obj.step / ( 1 / obj.sr ) );

% Find number of windows
nwindows = floor ( (nsamples - window_size) / step);

%% Compute the metrics

results = [];
results.nwindows = nwindows;

ft_progress('init', 'textbar', 'Please wait...')

for idx = 1:nwindows
    
    ft_progress(idx/nwindows, 'Processing event %d from %d', idx, nwindows)
    
    start_idx = (idx-1)*step + 1;
    end_idx = (idx-1)*step + window_size + 1;
    
    results.time_vect(idx) = (time_vect(start_idx)+time_vect(end_idx)) / 2;
    
    if isfield(obj,'true_timecourse')
        %results.true_timecourse(idx) = (obj.true_timecourse(start_idx)+obj.true_timecourse(end_idx)) / 2;
        results.true_timecourse(idx) = (obj.true_timecourse(floor((start_idx+end_idx)/2)));
    end
    
    glm = glm_estimator(theta_phase(start_idx:end_idx),gamma_amp(start_idx:end_idx));
    results.glm(:,idx) = glm.r2;
    results.plv(idx) = plv_estimator(theta_phase(start_idx:end_idx),gamma_amp_phase(start_idx:end_idx));
    results.esc(idx) = esc_estimator(theta(start_idx:end_idx),gamma_amp(start_idx:end_idx));
    results.aec(idx) = aec_estimator(theta_amp(start_idx:end_idx),gamma_amp(start_idx:end_idx));
    results.nesc(idx) = nesc_estimator(theta_phase(start_idx:end_idx),gamma_amp(start_idx:end_idx));
    results.mi(idx) = abs(est_canolty(theta_phase(start_idx:end_idx),gamma_amp(start_idx:end_idx),obj.sr));
    results.psi(idx) = psi_estimator(theta_phase(start_idx:end_idx),gamma_amp_theta(start_idx:end_idx),mean(obj.lo_bounds));
    
    if isfield(obj,'signal')    
        results.voytek(idx) = voytek_estimator(theta(start_idx:end_idx),gamma_amp_theta(start_idx:end_idx));
    else
        results.voytek(idx) = 0;
    end

    % metrics from mw
    S = struct( 'data_high', gamma(start_idx:end_idx), ...
                'data_low', theta(start_idx:end_idx), 'tres', 1/obj.sr,...
                'high_band',obj.hi_bounds, 'low_band',obj.lo_bounds);
    tmp = cross_freq_coupling(S);
    results.mark.plv(idx) = tmp.plv;
    results.mark.plv_cross(idx) = tmp.plv_cross;
    results.mark.pac_plv(idx) = tmp.pac_plv;
    results.mark.gs_plv(idx) = tmp.gs_plv;
    results.mark.gs_pac_plv(idx) = tmp.gs_pac_plv;

    
end

out = struct('signals',signals,'results',results');
