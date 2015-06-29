function cmg = estimate_comodulogram( obj )
%%
%
% Estimate the co-modulogram for a signal
%
% Must be passed a struct with:
%     sr: sampling rate
%     signal: 1d broadband signal to be analysed
%     obj.lo_bounds: frequency range for modulating signal (eg [5 7])
%     obj.lo_step: frequency steps for modulating signal (eg 1)
%     obj.lo_bandwidth: filter bandwidth for modulating signal (eg 2)
%     obj.hi_bounds: frequency range for modulated signal (eg [50 60])
%     obj.hi_step: frequency steps for modulated signal (eg 10)
%     obj.hi_bandwidth: bandwidth for modulated signal, either int or 'adaptive'
%
% and optionally:
%     true_timecourse: 1d signal indicating where pac exists
%     zero_pad: the number of samples to pad the time series when filtering
%     window_size: length in seconds for sliding window
%     step: step size between windows in seconds

%% Housekeeping

if nargin < 2
    verbose = false;
end

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

if isfield(obj,'window_size') && isfield(obj,'step')
    % Convert window and step to samples
    window_size = floor( obj.window_size / ( 1 / obj.sr ) );
    step = floor( obj.step / ( 1 / obj.sr ) );

    % Find number of windows
    nwindows = floor ( (nsamples - window_size) / step);
end


%% Create frequency ranges
%
% Create a [2 x nsteps] matrix of low frequency bounds to scan. lo_freqs(1,:)
% and lo_freqs(2,:) contain the high and low bounds to filter on based on the
% center frequencies in lo_bounds, lo_steps and lo_bandwidth. High frequency
% bandlimits are created in a [2, n_hi_steps, n_lo_steps] matrix. This is to
% allow for an adaptive bandwidth based on the lowest low frequency (if
% requested).

n_lo_steps = (obj.lo_bounds(2)-obj.lo_bounds(1))/obj.lo_step + 1;
lo_freqs = ones(2,n_lo_steps) .* repmat(obj.lo_bounds(1):obj.lo_step:obj.lo_bounds(2),2,1);
lo_freqs(1,:) = lo_freqs(1,:) - obj.lo_bandwidth/2;
lo_freqs(2,:) = lo_freqs(2,:) + obj.lo_bandwidth/2;

n_hi_steps = floor((obj.hi_bounds(2)-obj.hi_bounds(1))/obj.hi_step + 1);
hi_freqs = ones(2,n_hi_steps,n_lo_steps) .* repmat(obj.hi_bounds(1):obj.hi_step:obj.hi_bounds(2),[2,1,n_lo_steps]);

if strcmp(obj.hi_bandwidth,'adaptive')
    for idx = 1:n_lo_steps
        hi_bandwidth = lo_freqs(1,idx)/2;
        hi_freqs(1,:,idx) = hi_freqs(1,:,idx) - ones(1,n_hi_steps,1)*hi_bandwidth;
        hi_freqs(2,:,idx) = hi_freqs(2,:,idx) + ones(1,n_hi_steps,1)*hi_bandwidth;
    end
else
    size(hi_freqs)
    size(ones(n_hi_steps,n_lo_steps)*obj.hi_bandwidth/2)
     hi_freqs(1,:,:) = hi_freqs(1,:,:) - ones(1,n_hi_steps,n_lo_steps)*obj.hi_bandwidth/2;
     hi_freqs(2,:,:) = hi_freqs(2,:,:) + ones(1,n_hi_steps,n_lo_steps)*obj.hi_bandwidth/2;

end

% Preallocate metrics

esc = zeros(n_lo_steps,n_hi_steps);
nesc = zeros(n_lo_steps,n_hi_steps);
plv = zeros(n_lo_steps,n_hi_steps);
mi = zeros(n_lo_steps,n_hi_steps);
glm = zeros(n_lo_steps,n_hi_steps);
aec = zeros(n_lo_steps,n_hi_steps);


%%%%%%%%%%%%%%%%%%
%% Begin Main Loop
%%%%%%%%%%%%%%%%%%
msg = '';

for lo_idx = 1:n_lo_steps
    for hi_idx = 1:n_hi_steps

        fprintf(repmat(char(8),1,length(msg)));
        msg = sprintf('Computing low - %.2f:%.2f high - %.2f:%.2f', ...
            lo_freqs(1,lo_idx),lo_freqs(2,lo_idx),hi_freqs(1,hi_idx,lo_idx),hi_freqs(2,hi_idx,lo_idx));
        fprintf(msg);

        %% Create PAC signal -

        signals = make_pac_signals(obj.signal, ...
                                   obj.sr, ...
                                   [hi_freqs(1,hi_idx,lo_idx) hi_freqs(2,hi_idx,lo_idx)], ...
                                   [lo_freqs(1,lo_idx), lo_freqs(2,lo_idx)],...
                                   time_vect, ...
                                   obj.true_timecourse);



        if isfield(obj,'window_size')
            % Make sliding window data
            skip_field = {'hi_bounds','hi_bandwidth','hi_steps',...
                          'lo_bounds','lo_bandwidth','lo_steps','sr'};
            fields = fieldnames(signals);
            for i = 1:numel(fields)
                if strmatch(fields{i},skip_field)
                    continue
                 else
                    signals.(fields{i}) = make_sw_data(signals.(fields{i}),window_size,step);
                end
            end
        end
    
        % Estimate CFC
        esc(lo_idx,hi_idx) = esc_estimator(signals.theta,signals.gamma_amp);
        nesc(lo_idx,hi_idx) = nesc_estimator(signals.theta_phase,signals.gamma_amp);
        plv(lo_idx,hi_idx) = plv_estimator(signals.theta_phase,signals.gamma_amp_phase);
        glm(lo_idx,hi_idx) = glm_estimator(signals.theta_phase,signals.gamma_amp);
        mi(lo_idx,hi_idx) = mi_estimator(signals.theta_phase,signals.gamma_amp);
        aec(lo_idx,hi_idx) = aec_estimator(signals.theta_amp,signals.gamma_amp);

    end
end
fprintf('\n');
%% Create output struct

cmg = struct('lo_freqs',    lo_freqs,...
             'hi_freqs',    hi_freqs,...
             'esc',         esc,...
             'nesc',        nesc,...
             'plv',         plv,...
             'glm',         glm,...
             'mi',          mi,...
             'aec',         aec);
