function [ sw_signals ] = swsignals(signals,window_size,window_step)
%
% Function for splitting the data within a cfc signals struct into separate
% sliding windows based on a window length and step size (must be set in
% samples rather than ms)
%
% This turns signals struct containing [nchannels x nsamples] data into an
% [nchannels x windowsize x nwindows] matrix

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
        data = signals.(fields{i});

        sw_signals.(fields{i}) = make_sw(data,window_size,window_step);

    end
end
