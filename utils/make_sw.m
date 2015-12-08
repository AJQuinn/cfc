function out = make_sw( data, window_size, window_step )
%function out = make_sw( data, window_size, window_step )
%
% split a vector into a samples by windows matrix. window_size and window_step
% are given in samples.

    [nchannels,nsamples] = size(data); % there should only be one channel here
    nwindows = fix ( (nsamples - window_size) / window_step);

    out = zeros(window_size,nwindows);

    for idx = 1:nwindows
        start_idx = (idx-1)*window_step + 1;
        end_idx   = (idx-1)*window_step + window_size;

        out(:,idx) = squeeze(data(:,start_idx:end_idx));
    end
