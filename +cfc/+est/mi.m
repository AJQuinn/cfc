function [mean_vector] = mi(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006
%
% Modulating signal is [1 x samples x realisations]
% Modulated signal is [1 x samples x realisations]
%
% A realisation can be either an epoch or a sliding window

if ndims(modulating_signal) == 2
    ax = 2;
else
    ax = 1;
end

% Combine into complex form
Z = modulated_signal .* exp(1i * modulating_signal);

% Compute mean vector
mean_vector = abs(mean(squeeze(Z),ax));
