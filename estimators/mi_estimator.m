function [mean_vector] = mi_estimator(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006

% Combine into complex form
Z = modulated_signal .* exp(1i * modulating_signal);

% Compute mean vector
mean_vector = abs(mean(squeeze(Z)));
