function validate_cmg( signal,cfg )



%% Housekeeping
if ndims(signal) == 2
    [nchannels,nsamples] = size(signal);
    ntrials = 1;
else
    [nchannels,nsamples,ntrials] = size(signal);
end
if ~isfield(cfg,'filter_method');
    cfg.filter_method = 'twopass';
end
if ~isfield(cfg,'true_timecourse')
    cfg.true_timecourse = zeros(size(signal,1),size(signal,2));
end
if ~isfield(cfg,'pretrig')
    cfg.pretrig = 0;
end
if ~isfield(cfg,'full_estimate')
    cfg.full_estimate = 0;
end
% Generate time vector
time_vect = (0:1/cfg.sr:(nsamples-1) * (1/cfg.sr)) - cfg.pretrig;

n_lo_steps = (cfg.lo_bounds(2)-cfg.lo_bounds(1))/cfg.lo_step + 1;
lo_freqs = ones(2,n_lo_steps) .* repmat(cfg.lo_bounds(1):cfg.lo_step:cfg.lo_bounds(2),2,1);
lo_freqs(1,:) = lo_freqs(1,:) - cfg.lo_bandwidth/2;
lo_freqs(2,:) = lo_freqs(2,:) + cfg.lo_bandwidth/2;

n_hi_steps = floor((cfg.hi_bounds(2)-cfg.hi_bounds(1))/cfg.hi_step + 1);
hi_freqs = ones(2,n_hi_steps,n_lo_steps) .* repmat(cfg.hi_bounds(1):cfg.hi_step:cfg.hi_bounds(2),[2,1,n_lo_steps]);

if strcmp(cfg.hi_bandwidth,'adaptive')
    for idx = 1:n_lo_steps
        hi_bandwidth = lo_freqs(1,idx)+2;
        hi_freqs(1,:,idx) = hi_freqs(1,:,idx) - ones(1,n_hi_steps,1)*hi_bandwidth;
        hi_freqs(2,:,idx) = hi_freqs(2,:,idx) + ones(1,n_hi_steps,1)*hi_bandwidth;
    end
else
    size(hi_freqs)
    size(ones(n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2)
     hi_freqs(1,:,:) = hi_freqs(1,:,:) - ones(1,n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2;
     hi_freqs(2,:,:) = hi_freqs(2,:,:) + ones(1,n_hi_steps,n_lo_steps)*cfg.hi_bandwidth/2;

end


figure;subplot(111);hold on
for lo_idx = 1:n_lo_steps
    %subplot(5,6,lo_idx);hold on;grid on
    %set(gca,'YScale','log')
    title([num2str(lo_freqs(1,lo_idx)) ' ' num2str(lo_freqs(2,lo_idx))])
    %for hi_idx = 1:n_hi_steps

        if hi_freqs(1,hi_idx,lo_idx) <= lo_freqs(2,lo_idx) && ...
           cfg.full_estimate == 0
            % Skip estimation if the hi and lo bands overlap
            continue
        end

        if lo_freqs(2,lo_idx) < 0 ||  hi_freqs(1,hi_idx,lo_idx) < 0
            % skip estimation is any low frequency is below zero
            continue
        end

        %% Create PAC signal -
        lo = [lo_freqs(1,lo_idx), lo_freqs(2,lo_idx)];

        signals = cfc.util.basesignals(signal, ...
                                   cfg.sr, ...
                                   [hi_freqs(1,hi_idx,lo_idx) hi_freqs(2,hi_idx,lo_idx)], ...
                                   lo,...
                                   time_vect, ...
                                   cfg.true_timecourse,[],[],[],cfg.filter_method);

       if hi_idx == 1
           [pxx,f] = pwelch(signals.theta,[],[],[],cfg.sr);
           plot(f,pxx)
       end
       %[pxx,f] = pwelch(signals.gamma,[],[],[],cfg.sr);
       %plot(f,pxx)

end

