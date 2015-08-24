function out = pacman_sw(obj,verbose)
%
% Phase amplitude coupling measures are numerous
%
% Must be passed a struct with:
%     sr: sampling rate
%     window_size: length in seconds for sliding window
%     step: step size between windows in seconds
%     signal: 1d broadband signal to be analysed
%     lo_bounds: frequency range for modulating signal (eg [5 7])
%     hi_bounds: frequency range for modulated signal (eg [50 60])
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

if nargin < 2
    verbose = false;
end

% Make sure we haven't been given contradicting arguments
%if isfield(obj,'signal') && isfield(obj,'theta')
%    error('Please define either a signal and frequency bands or separate modulating and modulated signals');
%end

if ~isfield(obj,'pad')
    obj.pad = 100;
end
if ~isfield(obj,'true_timecourse')
    obj.true_timecourse = zeros(size(obj.signal,1),size(obj.signal,2));
end

% Generate time vector
if isfield(obj,'signal')
    [nsamples,~] = size(obj.signal);
end
time_vect = (0:1/obj.sr:(nsamples-1) * (1/obj.sr))';

if verbose;disp('Computing Signals...');end;

signals = cfc_util_basesignals(obj.signal,obj.sr,obj.hi_bounds,obj.lo_bounds,time_vect,obj.true_timecourse);

%% Set up sliding window

if verbose;disp('Creating sliding windows...');end

% Convert window and step to samples
window_size = floor( obj.window_size / ( 1 / obj.sr ) );
step = floor( obj.step / ( 1 / obj.sr ) );

% Find number of windows
nwindows = floor ( (nsamples - window_size) / step);

% Make sliding window data
skip_field = {'hi_bounds','lo_bounds','sr'};
fields = fieldnames(signals);
for i = 1:numel(fields)
    if strmatch(fields{i},skip_field)
        continue
     else
        signals.(fields{i}) = cfc_util_swsignals(signals.(fields{i}),window_size,step);
    end
end

%% Compute the metrics
metrics = {'esc','nesc','aec','plv','voytek','mi','glm'};

results = [];
if sum(ismember(metrics,'esc')) == 1
    if verbose;disp('Computing ESC...');end
    results.esc = esc_estimator(signals.theta,signals.gamma_amp);
end
if sum(ismember(metrics,'nesc')) == 1
    if verbose;disp('Computing NESC...');end
    results.nesc = nesc_estimator(signals.theta_phase,signals.gamma_amp);
end
if sum(ismember(metrics,'aec')) == 1
    if verbose;disp('Computing AEC');end
    results.aec = cfc_est_aec(signals.theta_amp,signals.gamma_amp);
end
if sum(ismember(metrics,'plv')) == 1
    if verbose;disp('Computing PLV');end
    results.plv = plv_estimator(signals.theta_phase,signals.gamma_amp_phase);
end
if sum(ismember(metrics,'mi')) == 1
    if verbose;disp('Computing MI...');end
    results.mi = mi_estimator(signals.theta_phase,signals.gamma_amp);
end
if sum(ismember(metrics,'glm')) == 1
    if verbose;disp('Computing GLM...');end
    results.glm = glm_estimator(signals.theta_phase,signals.gamma_amp);
end
if sum(ismember(metrics,'voytek')) == 1
    if verbose;disp('Computing voytek...');end
    results.voytek = voytek_estimator(signals.theta_phase,signals.gamma_amp);
end

results.nwindows = nwindows;
out.signals = signals;
out.results = results;
out.results.time_vect = mean(out.signals.time_vect,1);
out.results.true_timecourse = out.signals.true_timecourse(floor(window_size/2),:);
