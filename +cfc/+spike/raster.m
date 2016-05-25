function raster = raster(spike_times, spike_allocs, sample_rate, time_range)
%function out = cfc.spike.raster(spike_times, spike_allocs, sample_rate, time_range)
%
% Plot up a raster plot for a given time range
%
% spike_times: vectors of spike timings (in samples)
% spike_allocs: vector indicating which cluster each spike belongs to
% sample_rate
% time_range: times in seconds to display ie [ 0 100 ]
% clusters_idx: optional vector of clusters to use ie [ 3 5 11 12 ]

nclusters = double(max(spike_allocs));
nsamples = double(spike_times(end));

% Get time vectors and points of interest
time_vect = 0:1/sample_rate:nsamples/sample_rate;
time_idx = time_vect > time_range(1) & time_vect < time_range(2);

% Convert timings to samples
time_range(1) = floor(time_range(1) / (1 / sample_rate));
time_range(2) = floor(time_range(2) / (1 / sample_rate));

spike_idx = spike_times > time_range(1) & spike_times < time_range(2);

raster = nan(nclusters,sum(time_idx));

for idx = 1:nclusters
    % Get cluster of interest
    alloc_idx = spike_allocs == idx;

    % Get time point of each interesting spike in samples
    spk = spike_times(alloc_idx);


    raster(idx,spike_times(alloc_idx) & spike_idx) = 1 ;
end

size(time_vect(time_idx))


figure
for idx = 1:nclusters
    plot(time_vect(time_idx),raster(idx,:)+idx/2,...
          'k.','MarkerSize',3);hold on
end

