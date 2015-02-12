function signals = make_pac_signals(signal,sr,hi_bounds,lo_bounds,time_vect,true_timecourse)

[nsamples,~] = size(signal);

if nargin < 6
    true_timecourse = zeros(nsamples,1);
end

if nargin < 5
    time_vect = (0:1/sr:(nsamples-1) * (1/sr))';
end

pad = 150;

    % Compute frequencies of interest
pad_signal = [zeros(pad,1); signal; zeros(pad,1)];
theta = eegfilt_silent(pad_signal',sr,lo_bounds(1),lo_bounds(2),0,[],0,'fir1');
theta = theta(pad:end-pad-1)';
gamma = eegfilt_silent(pad_signal',sr,hi_bounds(1),hi_bounds(2),0,[],0,'fir1');
gamma = gamma(pad:end-pad-1)';

%% Compute signals

% Compute amplitude time series
theta_amp = abs(hilbert(theta));
gamma_amp = abs(hilbert(gamma));

% Compute phase time series
theta_phase = angle(hilbert(theta));

% Compute phase of gamma amplitude
gamma_amp_phase = angle(hilbert(gamma_amp));

% Compute theta-filtered gamma amplitude
pad_gamma = [zeros(pad,1); gamma_amp; zeros(pad,1)];
pad_gamma = eegfilt_silent(pad_gamma',sr,lo_bounds(1),lo_bounds(2),0,[],0,'fir1');
gamma_amp_theta = pad_gamma(pad:end-pad-1)';

%% TODO: This is wasteful, there is a lot in memory and only the filtering is
%time consuming. Eventually we want to be able to only save the filtered time
%courses and compute the simple stuff when it is loaded back in
signals = struct('theta',               theta, ...
                 'gamma',               gamma, ...
                 'gamma_amp',           gamma_amp, ...
                 'theta_amp',           theta_amp, ...
                 'theta_phase',         theta_phase, ...
                 'gamma_amp_phase',     gamma_amp_phase, ...
                 'gamma_amp_theta',     gamma_amp_theta, ...
                 'time_vect',           time_vect, ...
                 'true_timecourse',     true_timecourse, ...
                 'signal',              signal, ...
                 'lo_bounds',           lo_bounds, ...
                 'hi_bounds',           hi_bounds, ...
                 'sr',                  sr);
