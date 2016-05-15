
function [signals,gamma_cfg,theta_cfg] = cfc_util_basesignals(signal,sr,hi_bounds,lo_bounds,time_vect,true_timecourse,hi_trans,lo_trans,edge_width)
%% Create the ingredients for CFC metric estimation.
%
% signal is an array [channels x samples x realisations]. There may only be one
% or two channels within signal. If one channel is passed in both the
% modulating and modulated time series and generated from it. If two channels
% are passed in the modulating time series will be extracted from the first
% channel and the modulated time series will be extracted from the second.
%
% If the lo_bounds variable is one number, the modulating phase will be
% estimated with an interpolation method.

[nchannels,nsamples,nrealisations] = size(signal);
if nchannels > 2
    error('%s channels were passed in, two is the maximum', char(nchannels))
end

if nargin < 9 || isempty(edge_width)
    edge_width = 0;
end

if length(lo_bounds) == 2
    % We only need these if we're filtering for the modulating time series
    if (nargin < 8 || isempty(lo_trans))
        % Use a constant transition band
        lo_trans = [lo_bounds(1)-2 lo_bounds(2)+2];
    end
end
if nargin < 7 || isempty(hi_trans)
    hi_trans = [hi_bounds(1)-7.5 hi_bounds(2)+7.5];
end

if nargin < 6 || isempty(true_timecourse)
    true_timecourse = zeros(1,nsamples);
end

if nargin < 5 || isempty(time_vect)
    time_vect = (0:1/sr:(nsamples-1) * (1/sr));
end

if nsamples < sr*4
    % Should probably expose this at the top
    order = 128;
else
    order = 256;
end

%% preallocated
theta = zeros(nchannels,nsamples,nrealisations);
gamma = zeros(nchannels,nsamples,nrealisations);
theta_phase = zeros(nchannels,nsamples,nrealisations);
gamma_amp = zeros(nchannels,nsamples,nrealisations);
theta_amp = zeros(nchannels,nsamples,nrealisations);
gamma_amp_phase = zeros(nchannels,nsamples,nrealisations);
gamma_amp_theta = zeros(nchannels,nsamples,nrealisations);

%% Compute frequencies of interest

for idx = 1:nrealisations

    if length(lo_bounds) == 2
        theta_cfg.order = order;
        theta_cfg.sample_rate = sr;
        theta_cfg.centre_freq = (lo_bounds(1)+lo_bounds(2))/2;
        theta_cfg.pass_width = lo_bounds(2)-lo_bounds(1);
        theta_cfg.trans_width = lo_trans(2) - lo_trans(1);
        theta_cfg.method = 'twopass';

        theta(1,:,idx) = cfc_filt_fir(signal(1,:,idx),theta_cfg);
    else
        theta = [];
    end

    gamma_cfg.order = order;
    gamma_cfg.sample_rate = sr;
    gamma_cfg.centre_freq = (hi_bounds(1)+hi_bounds(2))/2;
    gamma_cfg.pass_width = hi_bounds(2)-hi_bounds(1);
    gamma_cfg.trans_width = hi_trans(2) - hi_trans(1);
    gamma_cfg.method = 'twopass';

    if nchannels == 1
        gamma(1,:,idx) = cfc_filt_fir(signal(1,:,idx),gamma_cfg);
    elseif nchannels == 2
        gamma(1,:,idx) = cfc_filt_fir(signal(2,:,idx),gamma_cfg);
    end

    %% Compute signals

    % Compute amplitude time series
    gamma_amp(1,:,idx) = abs(hilbert(gamma(1,:,idx)));
    theta_amp(1,:,idx) = abs(hilbert(theta(1,:,idx)));

    % Compute phase time series
    if length(lo_bounds) == 1
        % We only have one theta frequency, extract phase manually
        theta_phase(1,:,idx) = cfc_util_thetawaveform(signal(1,:,idx),lo_bounds,sr);
    else
        % We have filter bounds, extract phase with hilbert
        theta_phase(1,:,idx) = angle(hilbert(theta(1,:,idx)));
    end

    % Compute phase of gamma amplitude
    gamma_amp_phase(1,:,idx) = angle(hilbert(gamma_amp(1,:,idx)));

    % Compute theta-filtered gamma amplitude
    if length(lo_bounds) == 2
        gamma_amp_theta(1,:,idx) = cfc_filt_fir(gamma_amp(1,:,idx),theta_cfg);
    else
        gamma_amp_theta(1,:,idx) = cfc_util_thetawaveform(signal(1,:,idx),lo_bounds,sr);
    end

end

%% TODO: This is wasteful, there is a lot in memory and only the filtering is
%time consuming. Eventually we want to be able to only save the filtered time
%courses and compute the simple stuff when it is loaded back in
edge_width_samp = fix(sr * edge_width);
return_inds = (1+edge_width_samp):nsamples-edge_width_samp;

signals = struct('theta',               theta(:,return_inds,:), ...
                 'gamma',               gamma(:,return_inds,:), ...
                 'theta_amp',           theta_amp(:,return_inds,:),...
                 'gamma_amp',           gamma_amp(:,return_inds,:), ...
                 'theta_phase',         theta_phase(:,return_inds,:), ...
                 'gamma_amp_phase',     gamma_amp_phase(:,return_inds,:), ...
                 'gamma_amp_theta',     gamma_amp_theta(:,return_inds,:), ...
                 'time_vect',           time_vect(:,return_inds), ...
                 'true_timecourse',     true_timecourse(:,return_inds,:), ...
                 'signal',              signal(:,return_inds,:), ...
                 'window_len',          time_vect(end)-time_vect(1),...
                 'lo_bounds',           lo_bounds, ...
                 'hi_bounds',           hi_bounds, ...
                 'sr',                  sr);
