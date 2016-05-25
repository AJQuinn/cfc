function obj = cfc_simulate(S)

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


% TODO: Gamma amplitude not properly normalised??

obj.time_vect = [0:1/S.sample_rate:S.seconds];
obj.sr = S.sample_rate;
if isfield(S, 'switching_freq')
    obj.switching_freq = S.switching_freq;
end

if strcmp(S.method,'none')

        modulated_ts = sin(2*pi*S.modulated_freq*obj.time_vect);
        modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        sb_double = S.modulated_amp*(sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) + ...
                                          (2*pi*S.modulated_freq*obj.time_vect)) - ...
                                     sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) - ...
                                          (2*pi*S.modulated_freq*obj.time_vect)));

        state_switching = zeros(size(modulating_ts,1),size(modulating_ts,2));

        noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            modulated_ts);
        obj.modulated_ts = (modulated_ts + noise + sb_double.*state_switching);

        obj.signal = obj.modulated_ts + modulating_ts;

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

        return
elseif strcmp(S.method,'aq')

        modulated_ts = sin(2*pi*S.modulated_freq*obj.time_vect);
        modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        sb_double = S.modulated_amp*(sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) + ...
                                          (2*pi*S.modulated_freq*obj.time_vect)) - ...
                                     sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) - ...
                                          (2*pi*S.modulated_freq*obj.time_vect)));

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            modulated_ts);
        obj.modulated_ts = (modulated_ts + noise + sb_double.*state_switching);

        obj.signal = obj.modulated_ts + modulating_ts;

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

        return

elseif strcmp(S.method,'nonsinusoidal')

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

    noise = cfc_util_scalesignal(randn(size(obj.signal,1),size(obj.signal,2)),...
            S.noise_level,...
            obj.signal);
        obj.signal = obj.signal + noise;

    obj.signal = obj.signal + noise;      %Add noise to the signal.
    obj.signal = obj.signal(1:length(obj.time_vect));
    obj.signal = ( obj.signal - mean(obj.signal) );
    obj.state_switching = state_switching;
    return

elseif strcmp(S.method,'sawtooth')

    obj.signal = S.modulating_amp * cos( 2*pi*S.modulating_freq*obj.time_vect + .01);
    sawt = S.modulating_amp *  sawtooth( 2*pi*S.modulating_freq*obj.time_vect - .01);

    obj.state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    obj.signal(logical(obj.state_switching)) = sawt(logical(obj.state_switching));

    noise = cfc_util_scalesignal(randn(size(obj.signal,1),size(obj.signal,2)),...
            S.noise_level,...
            obj.signal);
    obj.signal = obj.signal + noise;

    return


elseif strcmp(S.method,'square')

    obj.signal = S.modulating_amp * cos( 2*pi*S.modulating_freq*obj.time_vect + .01);
    sq = S.modulating_amp *  square( 2*pi*S.modulating_freq*obj.time_vect - .01);

    obj.state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    obj.signal(logical(obj.state_switching)) = sq(logical(obj.state_switching));

    noise = cfc_util_scalesignal(randn(size(obj.signal,1),size(obj.signal,2)),...
            S.noise_level,...
            obj.signal);
    obj.signal = obj.signal + noise;

    return
elseif strcmp(S.method,'timedil')

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

    obj.state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    cont = true;
    obj.signal = dev;
    while cont == true
        if obj.state_switching(length(obj.signal)) == 1
            obj.signal = [obj.signal dev];
        else
            obj.signal = [obj.signal typ];
        end
        if length(obj.signal) + length(dev) > length(obj.time_vect);
            cont = false;
        end
    end

    obj.signal = padarray(obj.signal',size(obj.time_vect,2)-size(obj.signal,2),'post')';

    noise = cfc_util_scalesignal(randn(size(obj.time_vect,1),size(obj.time_vect,2)),...
            S.noise_level,...
            obj.signal);
    obj.signal = obj.signal + noise;

    return


elseif strcmp(S.method,'asymtimedil')

    template_seconds = 1/S.modulating_freq;
    template_vect = linspace(0,template_seconds,template_seconds*S.sample_rate);

    delta = 1/S.sample_rate;
    dil = .4;
    ddelta = (-dil/(S.modulating_freq*1.5)) * sin( 2*pi*S.modulating_freq*template_vect);

    for idx = 1:length(template_vect)
        t_vect(idx) = (delta*idx) + ddelta(idx);
    end

    dev = cos( 2*pi*S.modulating_freq*t_vect );
    typ = cos( 2*pi*S.modulating_freq*template_vect );

    dev = [ typ(1:fix(length(t_vect)/2)) dev((fix(length(t_vect)/2)+1):end) ];

    obj.state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    cont = true;
    obj.signal = dev;
    while cont == true
        if obj.state_switching(length(obj.signal)) == 1
            obj.signal = [obj.signal dev];
        else
            obj.signal = [obj.signal typ];
        end
        if length(obj.signal) + length(dev) > length(obj.time_vect);
            cont = false;
        end
    end

    obj.signal = padarray(obj.signal',size(obj.time_vect,2)-size(obj.signal,2),'post')';

    noise = cfc_util_scalesignal(randn(size(obj.time_vect,1),size(obj.time_vect,2)),...
            S.noise_level,...
            obj.signal);
    obj.signal = obj.signal + noise;


    return

elseif strcmp(S.method,'nonstat_mean')

    obj.signal = S.modulating_amp * sin( 2*pi*S.modulating_freq*obj.time_vect );
    state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

    mean_term = 1.5*max(sin(2*pi*1*obj.time_vect),0);
    obj.signal(state_switching == 1) = obj.signal(state_switching == 1) + ...
                                mean_term(state_switching == 1);

    noise = cfc_util_scalesignal(randn(size(obj.signal,1),size(obj.signal,2)),...
            S.noise_level,...
            obj.signal);

    obj.signal = obj.signal + noise;
    obj.state_switching = state_switching;
    return

elseif strcmp(S.method,'modal')

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

        noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        modulated_ts = sin(2*pi*S.modulated_freq*obj.time_vect);
        noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            modulated_ts);

        t = modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (modulated_ts.*t);
        obj.signal = obj.modulated_ts + modulating_ts;

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

        return

elseif strcmp(S.method,'noisemodal')

        r = .98;
        wr = 2*pi*(S.modulating_freq*1.57)/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        noise_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                noise_ts(idx) = noise_ts(idx) - squeeze(a1(ilag))*noise_ts(idx-ilag+1);
            end
        end

        noise_ts = ( noise_ts - mean(noise_ts) ) / std(noise_ts);
        modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

                sb_double = S.modulated_amp*(sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) + ...
                                          (2*pi*S.modulated_freq*obj.time_vect)) - ...
                                     sin((2*pi*S.modulating_freq*obj.time_vect+S.phase_lag) - ...
                                          (2*pi*S.modulated_freq*obj.time_vect)));

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        modulated_ts = sin(2*pi*S.modulated_freq*obj.time_vect);
        noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            modulated_ts);
        obj.modulated_ts = (modulated_ts + noise + sb_double.*state_switching);

        obj.signal = obj.modulated_ts + noise_ts + modulating_ts;

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

        return

elseif strcmp(S.method,'modalnoisemodal')

        r = .92;
        wr = 2*pi*(S.modulating_freq*2.23)/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        noise_ts = randn(1,length(time_vect));

        for idx = 4:length(time_vect)
            for ilag = 2:3
                noise_ts(idx) = noise_ts(idx) - squeeze(a1(ilag))*noise_ts(idx-ilag+1);
            end
        end

        noise_ts = ( noise_ts - mean(noise_ts) ) / std(noise_ts);

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

        noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
            S.noise_level,...
            modulating_ts);
        obj.modulating_ts = modulating_ts + noise;
        modulated_ts = sin(2*pi*S.modulated_freq*obj.time_vect);
        noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
            S.noise_level,...
            modulated_ts);

        t = modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (modulated_ts.*t);
        obj.signal = obj.modulated_ts + noise_ts + modulating_ts;

        obj.state_switching = state_switching;
        obj.time_vect = obj.time_vect;

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

noise = cfc_util_scalesignal(randn(size(modulating_ts,1),size(modulating_ts,2)),...
                            S.noise_level,...
                            modulating_ts);
obj.modulating_ts = modulating_ts + noise;
noise = cfc_util_scalesignal(randn(size(modulated_ts,1),size(modulated_ts,2)),...
                            S.noise_level,...
                            modulated_ts);
obj.modulated_ts = modulated_ts + noise;

obj.signal = signal;

obj.state_switching = state_switching;
obj.time_vect = obj.time_vect;


