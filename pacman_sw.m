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
    [nsamples,~] = size(obj.signal);
else
    [nsamples,~] = size(obj.theta);
end
time_vect = [0:1/obj.sr:(nsamples-1) * (1/obj.sr)]';

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
pad_gamma = [zeros(obj.pad,1); gamma_amp; zeros(obj.pad,1)];
pad_gamma = eegfilt_silent(pad_gamma',obj.sr,obj.lo_bounds(1),obj.lo_bounds(2),0,[],0,'fir1');
gamma_amp_theta = pad_gamma(obj.pad:end-obj.pad-1)';

signals = struct('theta', theta, ...
                 'gamma', gamma, ...
                 'gamma_amp',gamma_amp, ...
                 'theta_amp',theta_amp, ...
                 'theta_phase', theta_phase, ...
                 'gamma_amp_phase', gamma_amp_phase, ...
                 'gamma_amp_theta', gamma_amp_theta, ...
                 'time_vect', time_vect);
if isfield(obj,'true_timecourse')
    signals.true_timecourse = obj.true_timecourse;
else
    signals.true_timecourse = zeros(size(theta,1),size(theta,2));
end

%% Set up sliding window

% Convert window and step to samples
window_size = floor( obj.window_size / ( 1 / obj.sr ) );
step = floor( obj.step / ( 1 / obj.sr ) );

% Find number of windows
nwindows = floor ( (nsamples - window_size) / step);

% Make sliding window data
fields = fieldnames(signals);
for i = 1:numel(fields)
    signals.(fields{i}) = make_sw_data(signals.(fields{i}),window_size,step);
end

%% Compute the metrics
metrics = {'esc','nesc','aec','plv','voytek','mi','glm'};

results = [];
if sum(ismember(metrics,'esc')) == 1
    results.esc = esc_estimator(signals.theta,signals.gamma_amp);
end
if sum(ismember(metrics,'nesc')) == 1
    results.nesc = nesc_estimator(signals.theta,signals.gamma_amp);
end
if sum(ismember(metrics,'aec')) == 1
    results.aec = aec_estimator(signals.theta,signals.gamma_amp);
end
if sum(ismember(metrics,'plv')) == 1
    results.plv = plv_estimator(signals.theta_phase,signals.gamma_amp_phase);
end
if sum(ismember(metrics,'mi')) == 1
    results.mi = mi_estimator(signals.theta_phase,signals.gamma_amp);
end
if sum(ismember(metrics,'glm')) == 1
    results.glm = glm_estimator(signals.theta_phase,signals.gamma_amp);
end
if sum(ismember(metrics,'voytek')) == 1
    results.voytek = voytek_estimator(signals.theta_phase,signals.gamma_amp);
end

results.nwindows = nwindows;
out.signals = signals;
out.results = results;
out.results.time_vect = mean(out.signals.time_vect,1);
out.results.true_timecourse = mean(out.signals.true_timecourse,1);
