function [plv] = cfc_est_plv(modulating_signal,modulated_signal)
%
% Function for estimating the Phase locking value as defined in Lachaux 1999
% Modulating signal is the theta phase and modulated signal is the phase of the gamma amplitude
%
% Input data should be in arrays of [1 x nsamples x nrealisations]
% Output is a vector containing an estimate per epoch

if ndims(modulating_signal) == 2
    ax = 2;
else
    ax = 2;
end
nsamples = size(modulating_signal,2);

modulating_signal = squeeze(modulating_signal);
modulated_signal = squeeze(modulated_signal);

phase_diff = modulating_signal - modulated_signal;

plv =  abs(sum(exp(1i.*phase_diff),ax) ./ nsamples)';
