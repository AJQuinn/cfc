function [ out ] = cfc_util_stacktrials( signals, window_size, mode )
%% Switches a cfc_signals struct between [channels x samples x trials] and
% [channels x window_size*trials x nwindows]. Useful for task evoked CFC estimation
%
% When using 'stack' mode, this may reduce the number of samples in the
% dataset if the final window would be less than the specified window size
%
% signals: struct
%   cfc signals structure
% window_size: double
%   length of window in milliseconds
% mode: string
%   flag to indicate which transform to do
%   'stack': [channels x samples x trials] to [channels x window_size*trials x nwindows]
%   'unstack': [channels x window_size*trials x nwindows] to [channels x samples x trials]



% Convert window and step to samples
window_size = floor( window_size / ( 1000 / signals.sr ) );

if strcmp( mode, 'stack');
    [nchannels, nsamples, ntrials] = size(signals.signal);
    % Find number of windows
    nwindows = fix(nsamples / window_size);
    nsamples = window_size*nwindows;
elseif strcmp( mode, 'unstack' )
    [nchannels, tmp, nwindows] = size(signals.signal);
    ntrials = tmp / window_size;
    %ntrials = signals.ntrials;
    nsamples = ((nwindows-1)*window_size) + window_size;
    %nsamples = signals.nsamples;
else
    error('mode not recognised, please use "stack" or "unstack"');
end

fields = fieldnames(signals);

% We don't want to split these variables
skip_field = {'hi_bounds','hi_bandwidth','hi_steps',...
          'lo_bounds','lo_bandwidth','lo_steps','sr',...
          'sample_rate','switching_freq','nwindows',...
          'nsamples','ntrials','nchannels','window_size',...
          'time_vect','true_timecourse','window_len'};

out = [];
for i = 1:numel(fields)
    if strmatch(fields{i},skip_field)
        % Just add the variable without manipulating it
        out.(fields{i}) = signals.(fields{i});
    else
        tmp = signals.(fields{i});
        if strcmp(mode, 'stack')
            %[channels x samples x trials] to [channels x window_size*trials x nwindows]
            data = zeros(nchannels,window_size*ntrials,nwindows);
            for win_idx = 1:nwindows
                start_idx = (win_idx-1)*window_size + 1;
                end_idx = (win_idx-1)*window_size + window_size;
                data(:,:,win_idx) = reshape( tmp(:,start_idx:end_idx,:), nchannels, window_size*ntrials);
            end
            out.(fields{i}) = data;
        else
            %[channels x window_size*trials x nwindows] to [channels x samples x trials]
            data = zeros(nchannels,nsamples,ntrials);
            tmp = reshape( tmp, nchannels, window_size, ntrials, nwindows);
            for win_idx = 1:nwindows
                start_idx = (win_idx-1)*window_size + 1;
                end_idx = (win_idx-1)*window_size + window_size;
                data(:,start_idx:end_idx,:) = squeeze(tmp(:,:,:,win_idx));
            end
            out.(fields{i}) = data;
        end
    end
end

if strcmp(mode,'unstack')
    % We might have lost a sample
    out.time_vect = out.time_vect(1:size(out.theta,2));
end

out.nwindows = nwindows;
out.nsamples = nsamples;
out.ntrials = ntrials;
out.nchannels = nchannels;
out.window_size = window_size;

end
