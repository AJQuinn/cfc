function [ out ] = getslice( signals, start, stop )
%% Extract one window from a cfc_signals struct
%
% signals: struct
%   cfc signals structure
% start: double
%   start of the slice to extract
% stop: double
%   end of the slice to extract


idx = signals.time_vect > start & signals.time_vect < stop;

fields = fieldnames(signals);

% We don't want to split these variables
skip_field = {'hi_bounds','hi_bandwidth','hi_steps',...
          'lo_bounds','lo_bandwidth','lo_steps','sr',...
          'sample_rate','switching_freq'};

out = [];
for i = 1:numel(fields)
    if strmatch(fields{i},skip_field)
        % Just add the variable without manipulating it
        out.(fields{i}) = signals.(fields{i});
    else
        tmp = signals.(fields{i});
        out.(fields{i}) = tmp(:,start:stop,:);
    end
end
