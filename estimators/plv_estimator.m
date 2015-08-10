function [plv] = plv_estimator(modulating_signal,modulated_signal)
%
% Function for estimating the Phase locking value as defined in Lachaux 1999
% Modulating signal is the theta phase and modulated signal is the phase of the gamma amplitude
%
% Input data should be in arrays of [nchannels, nsamples]
% Output is a vector containing an estimate per epoch

[nchannels,nsamples] = size(modulating_signal);

phase_diff = modulating_signal - modulated_signal;

plv =  abs(sum(exp(1i.*phase_diff),2) ./ nsamples)';
