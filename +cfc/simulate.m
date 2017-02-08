function obj = simulate(S)
%%function obj = cfc.simulate(S)
% Generate a simulated signal with a known cross frequency coupling
%
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

% Housekeeping
obj.time_vect = [0:1/S.sample_rate:S.seconds];
obj.sr = S.sample_rate;

if isfield(S, 'switching_freq')
    obj.switching_freq = S.switching_freq;
    state_switching = square(sin(2*pi*S.switching_freq*obj.time_vect));
else
    state_switching = zeros(size(obj.time_vect));
end

if ~isfield(S,'randseed')
    rng('shuffle');
end

%% Development note
%
% Each function below should populate obj with a modulated and a modulating
% time-series along with any relevant state switching. These will be combined
% along with any scaled noise at the end.
%
% obj.modulating_ts
% obj.modulated_ts
%
% Any additional arguments required for the function should be commented here

if strcmp(S.method,'none')
        % two uncoupled signals, tonic low and high frequency oscillations

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        % ensure this is tonic
        state_switching = ones(size(modulating_ts,1),size(modulating_ts,2));

        obj.modulating_ts = obj.modulating_ts.*state_switching;

elseif strcmp(S.method,'basic')
        % A simple phase amplitude coupling between low and high frequency signals

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'sinsq')
        % phase-amplitude coupling between an asymmetrical (peak-to-trough) low
        % and a symmetrical high signal
        %
        % the degree of nonlinearity is set with S.nonlin_pow

        if ~isfield(S,'nonlin_pow')
            error('Please make sure S.nonlin_pow is defined to use this simulation')
        end

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        modulating_ts = (sin(2*pi*S.modulating_freq*obj.time_vect + pi)+1).^S.nonlin_pow;
        modulating_ts = -2*modulating_ts ./ max(modulating_ts) + 1;
        obj.modulating_ts = S.modulating_amp*modulating_ts;

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'sinsq_asym')
        % phase-amplitude coupling between a asymmetrical (peak-to-trough and
        % ascending-to-descending are asymmetrical) low and a symmetrical high
        % signal
        %
        % the degree of nonlinearity is set with S.nonlin_pow

        if ~isfield(S,'nonlin_pow')
            error('Please make sure S.nonlin_pow is defined to use this simulation')
        end


        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        % Generate two signals, one sinusoidal and one non-sinusoidal. These
        % are mixed within the waveform to generate the full waveform
        modulating_ts1 = (sin(2*pi*S.modulating_freq*obj.time_vect + pi)+1).^S.nonlin_pow;
        modulating_ts1 = -2*modulating_ts1 ./ max(modulating_ts1) + 1;
        modulating_ts2 = (sin(2*pi*S.modulating_freq*obj.time_vect + pi)+1).^1;
        modulating_ts2 = -2*modulating_ts2 ./ max(modulating_ts2) + 1;

        % Mix signals based on sinusoidal phase
        ang = angle(hilbert(sin(2*pi*S.modulating_freq*obj.time_vect + pi)));
        obj.modulating_ts(ang > 0) = modulating_ts1(ang>0);
        obj.modulating_ts(ang < 0) = modulating_ts2(ang<0);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'sawtooth')
        % Geneates a phase-amplitude coupling between a sawtooth and a high
        % frequency oscillation

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sawtooth(2*pi*S.modulating_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'square')
        % Geneates a phase-amplitude coupling between a sawtooth and a high
        % frequency oscillation

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sawtooth(2*pi*S.modulating_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'timedil')
        % Geneates a phase-amplitude coupling between a non-sinusoidal low
        % frequecy signal and a high frequency oscillation.
        % The non-sinusoidal signal is generated by perturbing the phase within
        % individual cycles to create wider peaks than troughs

        template_seconds = 1/S.modulating_freq;
        template_vect = linspace(0,template_seconds,template_seconds*S.sample_rate);

        delta = 1/S.sample_rate;
        dil = .25;
        ddelta = (-dil/(S.modulating_freq*1.5)) * sin( 2*pi*S.modulating_freq*template_vect);

        for idx = 1:length(template_vect)
            t_vect(idx) = (delta*idx) + ddelta(idx);
        end

        dev = cos( 2*pi*S.modulating_freq*t_vect );

        cont = true;
        obj.modulating_ts = dev;
        while cont == true
            obj.modulating_ts = [obj.modulating_ts dev];
            if length(obj.modulating_ts) + length(dev) > length(obj.time_vect);
                cont = false;
            end
        end

        obj.modulating_ts = padarray(obj.modulating_ts',size(obj.time_vect,2)-size(obj.modulating_ts,2),'post')';

        % Generate modulated time-series
        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'asymtimedil')
        % Geneates a phase-amplitude coupling between a non-sinusoidal low
        % frequecy signal and a high frequency oscillation.
        % The non-sinusoidal signal is generated by perturbing the phase within
        % individual cycles to create faster ascending edges than descending
        % and wider peaks than troughs

        template_seconds = 1/S.modulating_freq;
        template_vect = linspace(0,template_seconds,template_seconds*S.sample_rate);

        delta = 1/S.sample_rate;
        dil = .25;
        ddelta = (-dil/(S.modulating_freq*1.5)) * sin( 2*pi*S.modulating_freq*template_vect);

        for idx = 1:length(template_vect)
            t_vect(idx) = (delta*idx) + ddelta(idx);
        end

        dev = cos( 2*pi*S.modulating_freq*t_vect );
        typ = cos( 2*pi*S.modulating_freq*template_vect );

        % The deviant waveform is a mixture of the dilated waveform and a sinusoid
        dev = [ typ(1:fix(length(t_vect)/2)) dev((fix(length(t_vect)/2)+1):end) ];

        obj.state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        cont = true;
        obj.modulating_ts = dev;
        while cont == true
            obj.modulating_ts = [obj.modulating_ts dev];
            if length(obj.modulating_ts) + length(dev) > length(obj.time_vect);
                cont = false;
            end
        end

        obj.modulating_ts = padarray(obj.modulating_ts',size(obj.time_vect,2)-size(obj.modulating_ts,2),'post')';

        % Generate modulated time-series
        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'nonstat_mean')
        % A simple phase amplitude coupling between low and high frequency
        % signals in which the mean of the low-frequency signal drifts.
        %
        % Requires additional definition of
        % S.mean_freq

        if ~isfield(S,'mean_freq')
            disp('Using default mean frequency of .5Hz');
            S.mean_freq = .5;
        end

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);
        obj.modulating_ts = obj.modulating_ts + .25*sin(2*pi*S.mean_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);


elseif strcmp(S.method,'modal')
        % Generate a phase-amplitude-coupling signal between an autoregressive
        % oscillator and a high frequency signal. The low frequency signal has
        % a bandwidth introducing more realistic features such as envelope
        % dynamics.

        % Define AR oscillator
        r = .98;
        wr = 2*pi*S.modulating_freq/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        obj.modulating_ts = randn(1,length(time_vect));

        % Generate signal from model and parameters
        for idx = 4:length(time_vect)
            for ilag = 2:3
                obj.modulating_ts(idx) = obj.modulating_ts(idx) - squeeze(a1(ilag))*obj.modulating_ts(idx-ilag+1);
            end
        end

        % Normalise
        obj.modulating_ts = ( obj.modulating_ts - mean(obj.modulating_ts) ) / std(obj.modulating_ts);

        % Generate modulated signal
        obj.modulated_ts = S.modulated_amp.*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        % Modulate high frequency amplitude by magnitude of low frequency
        t = obj.modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (obj.modulated_ts.*t.*state_switching);

elseif strcmp(S.method,'basic_modalnoise')
        % Generate a simple phase-amplitude-coupling signal alongside an
        % unrelated oscillator close to the low frequency signal.

        % Geneate noise oscillator
        r = .98;
        wr = 2*pi*(S.modulating_freq*1.57)/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        obj.noise_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                obj.noise_ts(idx) = obj.noise_ts(idx) - squeeze(a1(ilag))*obj.noise_ts(idx-ilag+1);
            end
        end

        obj.noise_ts = ( obj.noise_ts - mean(obj.noise_ts) ) / std(obj.noise_ts);

        % Generate standard signals
        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect) + obj.noise_ts;

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'modal_modalnoise')
        % Generate a phase-amplitude-coupling signal between a low-frequency
        % oscillatory and a standard high frequency oscillation alongside an
        % unrelated oscillator close to the low frequency signal.

        r = .96;
        wr = 2*pi*(S.modulating_freq*2.57)/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        obj.noise_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                obj.noise_ts(idx) = obj.noise_ts(idx) - squeeze(a1(ilag))*obj.noise_ts(idx-ilag+1);
            end
        end

        obj.noise_ts = ( obj.noise_ts - mean(obj.noise_ts) ) / std(obj.noise_ts);

        r = .99;
        wr = 2*pi*S.modulating_freq/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        obj.modulating_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                obj.modulating_ts(idx) = obj.modulating_ts(idx) - squeeze(a1(ilag))*obj.modulating_ts(idx-ilag+1);
            end
        end

        obj.modulating_ts = ( obj.modulating_ts - mean(obj.modulating_ts) ) / std(obj.modulating_ts);
        obj.modulating_ts = obj.modulating_ts + obj.noise_ts;

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        t = obj.modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (obj.modulated_ts.*t.*state_switching);

elseif strcmp(S.method,'asymmodal')
        warning('Depreciated Option! - asymmodal');

        r = .98;
        wr = 2*pi*S.modulating_freq/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        modulating_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                modulating_ts(idx) = modulating_ts(idx) - squeeze(a1(ilag))*modulating_ts(idx-ilag+1);
            end
        end

        modulating_ts = ( modulating_ts - mean(modulating_ts) ) / std(modulating_ts);

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        amp = abs(hilbert(modulating_ts));
        ph = angle(hilbert(modulating_ts)*exp(-j*pi/2));
        di = diff(ph);
        di = di < -1;

        inds = find(di);
        ph2 = zeros(size(modulating_ts));
        for idx = 1:length(inds)-1

            start = inds(idx)+1;
            stop = inds(idx+1)+1;
            wavelen = stop-start;

            % Create asymmetrical modulation
            %mod = ( 1*linspace(1,.5,wavelen).^3 ).* ...
            %      1.*( abs(linspace(-1,1,wavelen).^1.5) ).* ...
            %      3.*sin(linspace(0,pi*2*1,wavelen));
            mod = ( 1*linspace(1,.5,wavelen).^2 ).* ...
                    1.*( abs(linspace(-1,1,wavelen).^1.7) ).* ...
                    4.*sin(linspace(0,pi*2*1,wavelen));

            % Modulate the phase dynamics
            ph2(start:stop-1) = -ph(start:stop-1) + mod;

        end
        modulating_ts = -sgolayfilt((amp.*(cos(ph2)))',2,21)';

        noise = cfc.util.scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        modulated_ts = S.modulated_amp.*sin(2*pi*S.modulated_freq*obj.time_vect);

        t = modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (modulated_ts.*t);
        obj.signal = obj.modulated_ts + modulating_ts;

        noise = cfc.util.scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            obj.signal);

        obj.signal = zscore(obj.signal + noise);

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

elseif strcmp(S.method,'kramer2008')

    % Stolen from supplementary materials of:
    % Kramer, Tort & Kopell (2008) Sharp edge artifacts and spurious
    % coupling in EEG frequency comodulation measures.
    % http://dx.doi.org/10.1016/j.jneumeth.2008.01.020

    f=S.modulating_freq;
    dt = 1/S.sample_rate;
    x1 = -cos(2*pi*(0:dt*f:.5));
    x2 = cos(2*pi*(0:dt*f:.5));
    theta = [x1 x2];                    %Define the low frequency (6 Hz) signal.
    state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    obj.signal = [];
    for k=1:length(theta):length(obj.time_vect)

        if state_switching(k) == 1
            fii = f*.8;
            r = ceil(rand()*(length(x1)/5));                %Define the duration of the sharp edge.
            x3i = ((0:r)/r);                      %Create the sharp edge.
            x3ii = (cos(2*pi*(0:dt*fii:.5))+1.0);%Define the taper of the sharp edge.
            x3 = [x3i x3ii] * .5;                      %Create the tapered sharp edge.

            pos = ceil((rand()*10));%+length(x1));        %Insert the sharp edge into the sinusoid.
            seed = theta;
            tmp = seed(pos:pos+length(x3)-1)+x3;
            tmp = tmp / 2;
            seed(pos:pos+length(x3)-1) = tmp;
        else
            r = rand*.25;
            seed = (theta*(1+r)) - (-1 - -(1-r));
        end
        obj.signal = [obj.signal seed];
    end

    noise = cfc.util.scalesignal(randn(size(obj.signal,1),size(obj.signal,2)),...
            S.noise_level,...
            obj.signal);
        obj.signal = obj.signal + noise;

    obj.signal = obj.signal + noise;      %Add noise to the signal.
    obj.signal = obj.signal(1:length(obj.time_vect));
    obj.signal = ( obj.signal - mean(obj.signal) );
    obj.state_switching = state_switching;
    return

elseif strcmp(S.method,'mw')

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

        signal = obj.modulating_ts - obj.modulated_ts;

elseif strcmp(S.method, 'hippo')

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

        signal = obj.modulating_ts - obj.modulated_ts;

elseif strcmp(S.method,'ab')

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

        signal = obj.modulating_ts - obj.modulated_ts;

end


%% Add noise scaled to the modulating time-series and return
obj.noise = cfc.util.scalesignal(randn(size(obj.modulated_ts,1),size(obj.modulated_ts,2)),...
                             S.noise_level,...
                             obj.modulated_ts);

obj.signal = obj.modulated_ts + obj.modulating_ts + obj.noise;

obj.state_switching = state_switching;

