function [plv] = plv_estimator(modulating_signal,modulated_signal)

% Function for estimating the phase amplitude coupling as defined in Lachaux 1999
% Modulating signal is the theta phase and the modulated signal is the phase of the gamma amplitude

[~,N] = size(modulating_signal);

phase_diff = modulating_signal - modulated_signal;

plv = abs( sum(exp(1j*phase_diff)) / N );

% Fisher transform - from mw's code
%plv=0.5*log((1+plv)/(1-plv));