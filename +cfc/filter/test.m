

sample_rate = 512;
seconds = 15;
X = randn(1,sample_rate*seconds);
time_vect = linspace(0,seconds,length(X));

X = X + sin(2*pi*10*time_vect + pi);

%% Simple filter

theta_cfg.order = sample_rate/2;
theta_cfg.sample_rate = sample_rate;
theta_cfg.centre_freq = 10;
theta_cfg.pass_width = 2;
theta_cfg.trans_width = 5;
theta_cfg.method = 'twopass';

% check responses
cfc.filt.explore(theta_cfg);


%% Run at several different trans bands

figure;subplot(211); hold on

plot(time_vect(1:sample_rate/2),X(1:sample_rate/2))

subplot(212); hold on

tb = [3 5 7 9];

tmp_cfg = theta_cfg;

for idx = 1:length(tb)
    tmp_cfg.trans_width = tb(idx);

    [theta,D{idx}] = cfc.filt.fir(X,tmp_cfg);

    subplot(211)
    plot(time_vect(1:sample_rate/2),theta(1:sample_rate/2))

    subplot(212)
    [h,w] = freqz(D{idx});
    plot(w,abs(h))

end

legend({'Original','4','6','8','10'})

%% Run at several different centre freqs

figure;subplot(211); hold on

plot(time_vect(1:sample_rate/2),X(1:sample_rate/2))

subplot(212); hold on

tmp_cfg = theta_cfg;

tb = [8 10 12 14 16];

for idx = 1:length(tb)
    tmp_cfg.centre_freq = tb(idx);

    [theta,D{idx}] = cfc.filt.fir(X,tmp_cfg);

    subplot(211)
    plot(time_vect(1:sample_rate/2),theta(1:sample_rate/2))

    subplot(212)
    [h,w] = freqz(D{idx});
    plot(w,abs(h))

end

legend({'Original','8','10','12','14','16'})

%% Run at several different centre freqs - ONE PASS

figure;subplot(211); hold on

plot(time_vect(1:sample_rate/2),X(1:sample_rate/2))

subplot(212); hold on

tmp_cfg = theta_cfg;
tmp_cfg.method = 'onepass';

tb = [6 8 10 12 14 16];

for idx = 1:length(tb)
    tmp_cfg.centre_freq = tb(idx);

    [theta,D{idx}] = cfc.filt.fir(X,tmp_cfg);

    subplot(211)
    plot(time_vect(1:sample_rate/2),theta(1:sample_rate/2))

    subplot(212)
    [h,w] = freqz(D{idx});
    plot(w,abs(h))

end

subplot(211)
legend({'Original','6','8','10','12','14','16'})


%% Simulate a signal

S = struct('seconds',           40, ...
           'sample_rate',       sample_rate,...
           'modulating_freq',   8,...
           'modulating_amp',    1,...
           'modulated_freq',    60,...
           'modulated_amp',     .25,...
           'noise_ratio',       .2,...
           'phase_lag',         pi/2,...
           'noise_level',       -10,...
           'switching_freq',    .25,...
           'method',            'aq');
signal = cfc_simulate(S);

figure;subplot(211); hold on

plot(time_vect(1:sample_rate/2),X(1:sample_rate/2))

subplot(212); hold on

gamma_cfg = theta_cfg;
gamma_cfg.centre_freq = 60;
gamma_cfg.pass_width = 10;
gamma_cfg.trans_width = [];

tb = [5 10 15 20 25];

for idx = 1:length(tb)
    gamma_cfg.pass_width = tb(idx);
    gamma_cfg.trans_width = tb(idx) + 3;


    [theta,D{idx}] = cfc.filt.fir(X,gamma_cfg);

    subplot(211)
    plot(time_vect(1:sample_rate/2),theta(1:sample_rate/2))

    subplot(212)
    [h,w] = freqz(D{idx});
    plot(w,abs(h))

end

legend({'Original','8','10','12','14','16'})
