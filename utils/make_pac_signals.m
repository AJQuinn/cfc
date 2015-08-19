function signals = make_pac_signals(signal,sr,hi_bounds,lo_bounds,time_vect,true_timecourse,hi_trans,lo_trans)
%% Create the ingredients for CFC metric estimation.
%
% signal can be a 1 or 2d array [channels x samples]. There may only be one or
% two channels within signal. If one channel is passed in both the modulating
% and modulated time series and generated from it. If two channels are passed
% in the modulating time series will be extracted from the first channel and
% the modulated time series will be extracted from the second.
%
% If the lo_bounds variable is one number, the modulating phase will be
% estimated with an interpolation method.

[nchannels,nsamples] = size(signal);
if nchannels > 2
    error('%s channels were passed in, two is the maximum', char(nchannels))
end

if length(lo_bounds) == 2
    % We only need these if we're filtering for the modulating time series
    if (nargin < 8 || isempty(lo_trans))
        % Use a constant transition band
        lo_trans = [lo_bounds(1)-5 lo_bounds(2)+5];
    end
end
if nargin < 7 || isempty(hi_trans)
    hi_trans = [hi_bounds(1)-7.5 hi_bounds(2)+7.5];
end

if nargin < 6 || isempty(true_timecourse)
    true_timecourse = zeros(nsamples,1);
end

if nargin < 5 || isempty(time_vect)
    time_vect = (0:1/sr:(nsamples-1) * (1/sr));
end

%% Compute frequencies of interest

if length(lo_bounds) == 2
    theta_cfg.order = 512;
    theta_cfg.sample_rate = sr;
    theta_cfg.centre_freq = (lo_bounds(1)+lo_bounds(2))/2;
    theta_cfg.pass_width = lo_bounds(2)-lo_bounds(1);
    theta_cfg.trans_width = lo_trans(2) - lo_trans(1);
    theta_cfg.method = 'twopass';

    theta = fir_filter_data(signal(1,:),theta_cfg);
else
    theta = [];
end

gamma_cfg.order = 512;
gamma_cfg.sample_rate = sr;
gamma_cfg.centre_freq = (hi_bounds(1)+hi_bounds(2))/2;
gamma_cfg.pass_width = hi_bounds(2)-hi_bounds(1);
gamma_cfg.trans_width = hi_trans(2) - hi_trans(1);
gamma_cfg.method = 'twopass';

if nchannels == 1
    gamma = fir_filter_data(signal(1,:),gamma_cfg);
elseif nchannels == 2
    gamma = fir_filter_data(signal(2,:),gamma_cfg);
end

%% Compute signals

% Compute amplitude time series
gamma_amp = abs(hilbert(gamma));

% Compute phase time series
if length(lo_bounds) == 1
    % We only have one theta frequency, extract phase manually
    theta_phase = theta_waveform(signal(1,:),lo_bounds,sr);
else
    % We have filter bounds, extract phase with hilbert
    theta_phase = angle(hilbert(theta));
end

% Compute phase of gamma amplitude
gamma_amp_phase = angle(hilbert(gamma_amp));

% Compute theta-filtered gamma amplitude
if length(lo_bounds) == 2
    gamma_amp_theta = fir_filter_data(gamma_amp(1,:),theta_cfg);
else
    gamma_amp_theta = theta_waveform(signal(1,:),lo_bounds,sr);
end

%% TODO: This is wasteful, there is a lot in memory and only the filtering is
%time consuming. Eventually we want to be able to only save the filtered time
%courses and compute the simple stuff when it is loaded back in
signals = struct('theta',               theta, ...
                 'gamma',               gamma, ...
                 'gamma_amp',           gamma_amp, ...
                 'theta_phase',         theta_phase, ...
                 'gamma_amp_phase',     gamma_amp_phase, ...
                 'gamma_amp_theta',     gamma_amp_theta, ...
                 'time_vect',           time_vect, ...
                 'true_timecourse',     true_timecourse, ...
                 'signal',              signal, ...
                 'lo_bounds',           lo_bounds, ...
                 'hi_bounds',           hi_bounds, ...
                 'sr',                  sr);
