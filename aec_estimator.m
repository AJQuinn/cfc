function [aec] = aec_estimator(modulating_signal,modulated_signal)

% Function for estimating the Amplitude envelope-to-signal correlation as defined in Bruns and Eckhorn 2004
% Modulating signal is the theta filtere time series and modulated_signal is the gamma amplitude

aec = corr(modulating_signal',modulated_signal');
