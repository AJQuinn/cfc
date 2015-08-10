function obj = generate_pac_signal(S)

% S.seconds
% S.sample_rate
% S.modulating_freq
% S.modulating_amp
% S.modulated_freq
% S.modulated_amp
% S.phase_lag
% S.method
% S.noise_ratio
% S.theta

%if ~isfield('S','phase_lag')
%    S.phase_lag = pi/2; % Lock to peak of theta wave
%end

obj.time_vect = [0:1/S.sample_rate:S.seconds];
obj.sample_rate = S.sample_rate;
if isfield(S, 'switching_freq')
    obj.switching_freq = S.switching_freq;
end

obj

if S.method == 'aq'

        base_gamma = .1;

        % Generate coupled signals
        modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);
        modulated_ts = S.modulating_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        mod = modulating_ts;
        mod(mod < base_gamma) = base_gamma;
        % Make modulation time-varying
        state_switching = square(sin(2*pi*obj.switching_freq*obj.time_vect));
        %mod = mod .*state_switching;

        % Create phase scramble of original modulating signal
        s = phaseran(modulating_ts,1)';
        s(s < base_gamma) = base_gamma;

        % Scramble scrambled series in chunks
        chunk_size = 50;
        nchunks = floor(size(s,2) / chunk_size);
        s_chunk = reshape(s(1,1:chunk_size*nchunks),[chunk_size,nchunks]);
        s_chunk = s_chunk(:,randperm(size(s_chunk,2)));
        s(1,1:nchunks*chunk_size) = reshape(s_chunk,[1,chunk_size*nchunks]);

        % Where modulating signal is zero we substitute in the phase scrambled
        % modulating signal
        mod(state_switching < 1) = s(1,1:size(find(state_switching < 1),2));
        obj.mod = mod;

        modulated_ts = modulated_ts.*mod;

elseif S.method == 'mw'

        % Generate coupled signals
        Ttheta=1/S.modulating_freq;
        Tgamma=1/S.modulated_freq; %secs
        modulated_ts = S.modulated_amp * sin(2*pi*obj.time_vect/Tgamma)/2;
        if isfield(S,'theta')
            modulating_ts = sin(S.theta)';
        else
            modulating_ts= S.modulating_amp* sin(2*pi*obj.time_vect/Ttheta)/2;
        end
        % Vary coupling over time
        state_switching = (square(2*pi*obj.time_vect/obj.switching_freq)+1)/2;

        k=3;c=3;
        modulating_ts_lag=sin(2*pi*obj.time_vect/Ttheta + S.phase_lag)/2;
        amp_gamma=state_switching.*k./(1 + exp(-c*(modulating_ts_lag)));
        modulated_ts=amp_gamma.*modulated_ts;

        % Normalise variance of the modulated signal when not modulated
        amptmp=std(modulated_ts(find(modulated_ts~=0)));
        modulated_ts(find(modulated_ts==0))=randn(size(find(modulated_ts==0)))*amptmp;

elseif S.method == 'hippo'

        % Generate coupled signals
        Ttheta=1/S.modulating_freq;
        Tgamma=1/S.modulated_freq; %secs
        modulated_ts = S.modulated_amp * sin(2*pi*obj.time_vect/Tgamma)/2;
        modulating_ts= S.modulating_amp* load_pac_theta(S.sample_rate);

        % Vary coupling over time
        state_switching = (square(2*pi*obj.time_vect/obj.switching_freq)+1)/2;

        k=3;c=3;
        modulating_ts_lag=sin(2*pi*obj.time_vect/Ttheta + S.phase_lag)/2;
        amp_gamma=state_switching.*k./(1 + exp(-c*(modulating_ts_lag)));
        modulated_ts=amp_gamma.*modulated_ts;

        % Normalise variance of the modulated signal when not modulated
        amptmp=std(modulated_ts(find(modulated_ts~=0)));
        modulated_ts(find(modulated_ts==0))=randn(size(find(modulated_ts==0)))*amptmp;


elseif S.method == 'ab'

        % Generate coupled signals
        windows = 100;
        win_vect = 0:1/S.sample_rate:obj.time_vect(end)/windows;
        winsize = length(win_vect);

        lfwave = sin(angle(hilbert(sin(2*pi*S.modulating_freq*win_vect))));
        hfwave = sin(2*pi*S.modulated_freq*win_vect);
        hfmod  = (2*(lfwave+1)/2+0.5).*cos(angle(hilbert(hfwave)));

        m1 = normalise(gausswin(winsize).* lfwave');
        m2 = normalise(gausswin(winsize).* hfmod');

        % Vary coupling over time
        modulating_ts = []; modulated_ts = []; a = [];

        for n = 1:100
                if randi(4) == 1
                    modulating_ts = [modulating_ts; m1];
                    modulated_ts = [modulated_ts; m2];
                    a  = [a;  gausswin(winsize)];
                else
                    modulating_ts = [modulating_ts; zeros(winsize,1)];
                    modulated_ts = [modulated_ts; zeros(winsize,1)];
                    a  = [a;  zeros(winsize,1)];
                end
        end

        state_switching = a>.1;

end

noise = scale_signal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
                            S.noise_level,...
                            modulating_ts);
obj.modulating_ts = modulating_ts + noise;
noise = scale_signal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
                            S.noise_level,...
                            modulated_ts);
obj.modulated_ts = modulated_ts + noise;

obj.signal = obj.modulating_ts - obj.modulated_ts;

obj.state_switching = state_switching;
obj.time_vect = obj.time_vect;


