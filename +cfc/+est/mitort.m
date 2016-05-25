function [mi] = mitort(modulating_signal,modulated_signal)
% Function for estimating the modulation index in Tort et al 2008 PNAS
% Modulating signal is the theta phase time series and modulated_signal is the gamma amplitude
%
% Modulating signal is [1 x samples x realisations]
% Modulated signal is [1 x samples x realisations]
%
% A realisation can be either an epoch or a sliding window

if ndims(modulating_signal) == 2
    modulating_signal = modulating_signal';
    modulated_signal = modulated_signal';
else
    modulating_signal = squeeze(modulating_signal);
    modulated_signal = squeeze(modulated_signal);
end

mi = zeros(1,size(modulating_signal,2));
nbin = 21;


for itrial = 1:size(modulating_signal,2)

    [~,bins] = hist(modulating_signal,nbin);
    mean_amp = zeros(1,20);
    for idx = 1:20
        inds = bins(idx) < modulating_signal(:,itrial) & bins(idx+1) > modulating_signal(:,itrial);
        mean_amp(idx) = mean(modulated_signal(inds,itrial));
    end

    mean_amp(isnan(mean_amp)) = 0;

    % https://github.com/cineguerrilha/Neurodynamics/blob/master/16ch/Comodulation/ModIndex_v1.m
    mi(itrial) = (log(nbin)-(-sum((mean_amp/nansum(mean_amp,2)).*log((mean_amp/nansum(mean_amp,2))))))/log(nbin);

end
