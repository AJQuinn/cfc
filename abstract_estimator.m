function [plv] = plv_estimator(modulating_signal,modulated_signal)

% Function for estimating the Phase locking value as defined in Lachaux 1999
% Modulating signal is the theta phase and modulated signal is the phase of the gamma amplitude

[~,N] = size(modulating_signal);

phase_diff = modulating_signal - modulated_signal;

plv = sum(exp(phase_diff)) / N;

