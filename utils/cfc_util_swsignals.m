function [ sw_signals ] = cfc_util_swsignals(signals,window_size,step)
%
% Function for splitting the data within a cfc signals struct into separate %
% sliding windows based on a window length and step size (must be set in
% samples rather than ms)
%
% This turns signals struct containing [nchannels x nsamples] data into an
% [nchannels x windowsize x nwindows] matrix

[nchannels,nsamples] = size(signals.signal); % there should only be one channel here
nwindows = fix ( (nsamples - window_size) / step);

fields = fieldnames(signals);

% We don't want to split these variables
skip_field = {'hi_bounds','hi_bandwidth','hi_steps',...
          'lo_bounds','lo_bandwidth','lo_steps','sr',...
          'sample_rate','switching_freq'};

%% Main loop
for i = 1:numel(fields)
    if strmatch(fields{i},skip_field)
        % Just add the variable without manipulating it
        sw_signals.(fields{i}) = signals.(fields{i});
    else
        % Make sliding window version of variable
        sw_data = zeros(nchannels,window_size,nwindows);
        data = signals.(fields{i});
        for idx = 1:nwindows
            start_idx = (idx-1)*step + 1;
            end_idx = (idx-1)*step + window_size;

            sw_data(1,:,idx) = squeeze(data(:,start_idx:end_idx));
        end

        sw_signals.(fields{i}) = sw_data;

    end
end

