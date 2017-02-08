function ret = peak_lock( signal, time_vect, sample_rate, filt_cfg, winsize, peaklock, edges )
%function ret = peak_lock( signal, time_vect, sample_rate, normalise, resolution )
%
% Function for epoching a continuous time-series around either peaks or troughs
% within a specific frequency band
%

if nargin < 3 || isempty(sample_rate)
    sample_rate = ceil(1/diff(time_vect(1:2)));
end

if nargin < 5 || isempty(winsize)
    winsize = fix(sample_rate/3);
end

if nargin < 6 || isempty(peaklock)
    disp('hey');
    peaklock = true;
end

if nargin < 7 || isempty( edges )
    edges = fix(sample_rate/2);
end

% If we're passed the full filter specification, use that, otherwise define a
% filt_specification based on the frequency info provided. This should be a low
% and high value to filter around i.e. filt_cfg = [8 13]
if ~isstruct( filt_cfg ) && numel( filt_cfg ) == 2

    freqs = filt_cfg;

    % pks are troughs...
    filt_cfg = [];
    filt_cfg.sample_rate = sample_rate;
    filt_cfg.centre_freq = mean(freqs);
    filt_cfg.pass_width = diff(freqs);
    filt_cfg.trans_width = diff(freqs) + 6;
    filt_cfg.method = 'twopass';
    filt_cfg.order = [];

    % Filter signal
    [filt_signal,~] = cfc.filter.fir(signal,filt_cfg,[],[]);
    filt_signal = filt_signal(edges:end-edges);

    edges_applied = true;

elseif ~isstruct( filt_cfg )
    % If this isn't a filter or just two values, use it as reference
    if size(signal) == size(filt_cfg)
        filt_signal = filt_cfg;
        clear filt_cfg;
    else
        disp('Sizes dont match');
    end
        edges_applied = false;

else

    % Filter signal
    [filt_signal,~] = cfc.filter.fir(signal,filt_cfg,[],[]);
    filt_signal = filt_signal(edges:end-edges);
    edges_applied = true;

end

% Find peaks
% if peaklock == true
%     disp('peak');
%     locs = findpeaks(filt_signal);
%     pks = filt_signal(locs);
% else
%     % or troughs
%     disp('trough');
%     locs = findpeaks(-filt_signal);
%     pks = filt_signal(locs);
% end

% Find peaks
if peaklock == true
    [locs,pks] = cfc.peak.peakfinder(filt_signal,[],[],[],false);
else
    % or troughs
    [locs,pks] = cfc.peak.peakfinder(-filt_signal,[],[],[],false);
    pks = -pks;
end

% account for filter edge removal
%locs = locs + edges;



% Find epochs around each peak or trough
trls = zeros(2,length(locs));
trls(1,:) = locs-winsize;
trls(2,:) = locs+winsize;
% sanity check
for ind = 1:sum(sum(trls<=0))
    if trls(1,1) <= 0;
        trls = trls(:,2:end);
        pks = pks(2:end);
        locs = locs(2:end);
    end
end
for ind = 1:sum(sum(trls>=length(signal)-sample_rate))
    if trls(2,end) >= (length(signal)-sample_rate);
        trls = trls(:,1:end-1);
        pks = pks(1:end-1);
        locs = locs(1:end-1);
    end
end

if edges_applied == true
    % add edges back on
    locs = locs + edges;
    trls = trls + edges;
end

% only take peaks with some sensible magnitude
zmag = zscore(pks);
% locs = locs( zmag > -1 );
% pks = pks( zmag > -1 );

% average waveform
waveform = zeros(size(trls,2),2*winsize + 1);
for idx = 1:size(trls,2)
    waveform(idx,:) = signal(trls(1,idx):trls(2,idx));
end
avg_waveform = mean(waveform,1);

% Define output structure
ret = [];
ret.avg_waveform = avg_waveform;
ret.waveform = waveform;
ret.trls = trls;
ret.winsize = winsize;
ret.avg_time = linspace(-(winsize/sample_rate),(winsize/sample_rate),size(avg_waveform,2));
ret.pks = pks;
ret.locs = locs;

