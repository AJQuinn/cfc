function [mean_vector] = canolty_estimator(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006
% Modulating signal is the theta phase and modualted signal is the gamma amplitude


% Combine into complex form
Z = modulating_signal .* exp(1i * modulated_signal);

% Compute mean vector
mean_vector = abs(mean(squeeze(Z)));
