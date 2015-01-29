function [nesc] = nesc_estimator(modulating_signal,modulated_signal)

% Function for estimating the normalised Envelope-to-signal correlation as defined in Penny et al 2008
% Modulating signal is the theta filtere time series and modulated_signal is the gamma amplitude

nesc = corr(cos(modulating_signal)',modulated_signal');
