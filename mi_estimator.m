function [mean_vector] = canolty_estimator(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006

% Get the phase time series from the modulating signal
phase = angle(hilbert(modulating_signal));

% Get the amplitude time series from the modualted signal
amp = abs(hilbert(modulated_signal));

% Combine into complex form
Z = amp .* exp(1i * phase);

% Compute mean vector
mean_vector = abs(mean(squeeze(Z)));
