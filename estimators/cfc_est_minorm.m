function [mean_vector] = cfc_est_minorm(modulating_signal,modulated_signal)
%
% Function estimating the mean vector as defined in Canolty et al 2006. The MI
% is returned as a z-value on a null distribution created by phase shifted
% surrogates. See supplimentary materials of Canolty et al 2006.
%
% Modulating signal is [1 x samples x realisations]
% Modulated signal is [1 x samples x realisations]
%
% A realisation can be either an epoch or a sliding window

if ndims(modulating_signal) == 2
    ax = 2;
    [~,nsamples] = size(modulating_signal);
    nrealisations = 1;
else
    ax = 1;
    modulating_signal = squeeze(modulating_signal);
    modulated_signal = squeeze(modulated_signal);
    [nsamples,nrealisations] = size(modulating_signal);
end

%modulating_signal = modulating_signal - mean(modulating_signal);

% Combine into complex form
Z = modulated_signal .* exp(1i * modulating_signal);

% Compute mean vector
mean_vector = squeeze(mean(Z,ax));

% Compute surrogates
nsurrogates = 200;
surrogate_Z = nan(nsurrogates,nrealisations);
skips = randi([fix(nsamples/10),9*fix(nsamples/10)],nsurrogates,1);
for idx = 1:nsurrogates
    tmp = circshift(modulated_signal,skips(idx),ax);
    surrogate_Z(idx,:) = abs(mean( tmp.*exp(1i * modulating_signal),ax));
end

% Normalise
Z_length = (abs(mean_vector)-mean(surrogate_Z,1))./std(surrogate_Z,[],1);
mean_vector = abs(Z_length.*exp(1i*angle(mean_vector)));
