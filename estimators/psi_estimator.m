function [psi] = psi_estimator(modulating_signal,modulated_signal,theta_freq)

% Function for estimating the phase slope index between two signals using the function retrieved from: http://doc.ml.tu-berlin.de/causality/
% Modulating signal is the theta phase and the modulated signal is the phase of the gamma amplitude


srate = 400; % ASSUMED FOR NOW!!!
z = 2*pi*1i/srate;

df = (srate / 2) / length(modulating_signal);

% Compute Coherence
a = modulating_signal .* exp(-z*(1:length(modulating_signal))*theta_freq);
b = modulated_signal .* exp(-z*(1:length(modulated_signal))*theta_freq);

Sa = a * conj(a)';
Sb = b * conj(b)';
Sab = a * conj(b)';

C_f = Sab / sqrt(Sa * Sb);

% Compute Coherence at new freq
a = modulating_signal .* exp(-z*(1:length(modulating_signal))*theta_freq+df);
b = modulated_signal .* exp(-z*(1:length(modulated_signal))*theta_freq+df);

Sa = a * conj(a)';
Sb = b * conj(b)';
Sab = a * conj(b)';

C_df = Sab / sqrt(Sa * Sb);

% Compute PSI
psi = squeeze( imag( conj(C_f) * C_df ) );