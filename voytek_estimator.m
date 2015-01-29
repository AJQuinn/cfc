function [pac] = voytek_estimator(modulating_signal,modulated_signal)

% Function for estimating phase amplitude coupling as defined in Voytek et al 
% The modulated signal here should be the amplitude time-series of the high frequency component filtered to the same frequencies as the modulating signal.

% Correlate the signals
pac = corr(modulating_signal', modulated_signal');
