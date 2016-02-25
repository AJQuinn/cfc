function cfc_util_inspectsig( signals, freq_max, freq_scale )
%%function cfc_util_inspectsig( signals )
%
% This function will compute a range of CFC relevant metrics from a dataset
% which has been prepared with cfc_util_basesignal and produce a summary plot.
%
% The left of the plot will contain the power spectrum with the selected
% modulatingand modulating frequencyranges highlighted.
%
% The right of the plot displays:
% Top row - Modulating phase distribution
%         - Modulated power distribution
%         - Average power by modulating phase
% Bottom row - Phase trajectory in complex plan
%            - Amplitude trajectory in complex plan (may not be meaningful)
%            - Amplitude modulated phase trajectory (Canolty et al 2006)
%

if nargin < 3 || isempty(freq_scale)
    freq_scale = 'dB'
end

if nargin < 2 || isempty(freq_max)
    freq_max = signals.sr / 2;
end

winlen = signals.time_vect(end) - signals.time_vect(1);

% If we have been given windowed or epoched data, stack everything up
% This is a bit of a hack for the pwelch below, really should estimate this per
% epoch to avoid edges
stacked = 0;
if ndims ( signals.signal ) == 3
    signals = cfc_util_stacktrials(signals,signals.window_len *1000,'stack');
    stacked = 1;
end

figure('position',[100 100 1024, 512]);

%% Plot full data power spectrum
subplot(3,4,[1 2 5 6]);hold on;grid on;
[pow,f] = pwelch(signals.signal,[],[],[],signals.sr);
[pow_lo,~] = pwelch(signals.theta,[],[],[],signals.sr);
[pow_hi,~] = pwelch(signals.gamma,[],[],[],signals.sr);

if strcmp(freq_scale,'log')
    pow = log(pow);
    pow_lo = log(pow_lo);
    pow_hi = log(pow_hi);
elseif strcmp(freq_scale,'dB')
    pow = 20*log(pow);
    pow_lo = 20*log(pow_lo);
    pow_hi = 20*log(pow_hi);
end

plot(f,pow);
% add in frequencies of interest

yval = [min(pow) min(pow)];
plot(signals.lo_bounds,yval,'r','linewidth',10);
plot(signals.hi_bounds,yval,'g','linewidth',10);
axis('tight');
legend({'Signal','Modulating Range','Modulated Range'});
title('Power Spectrum and Frequency ranges')
xlabel('Frequency (Hz)')
ylabel('Power');
xlim([0 freq_max]);
yvals = ylim;

subplot(3,4,[9 10]);hold on; grid on;
plot(f,pow_lo);
plot(f,pow_hi);
legend({'Modulating Signal','Modulated Signal'});
title('Filtered Signal Power Spectra')
xlabel('Frequency (Hz)')
ylabel('Power');
xlim([0 freq_max])
ylim(yvals);




%% Phase and power distributions
subplot(343);hold on;grid on;
[counts,bins] = hist(signals.theta_phase,20);
hist(signals.theta_phase,20);
plot([-pi pi],[mean(counts) mean(counts)],'k--','linewidth',2);
xlabel('Modulating Phase');
xlim([-pi pi]);
title('Phase Distribution');

subplot(3,4,7);hold on;grid on;
hist(signals.gamma_amp,20);
xlabel('Modulated Power');
ylabel('Amplitude Envelope');
title('Power Distribution');

%% Phase on unit circle
subplot(3,4,4);grid on;hold on;
plot(sin(signals.theta_phase),cos(signals.theta_phase))
plot([0 mean(sin(signals.theta_phase))],[0 mean(cos(signals.theta_phase))],'r','linewidth',5)
xlabel('Real (cos)');xlim([-1.2 1.2])
ylabel('Imag (sin)');ylim([-1.2 1.2])
title('Modulating phase in z-plane')

%% Amp on unit circle
subplot(3,4,8);grid on;hold on;
plot(sin(signals.gamma_amp),cos(signals.gamma_amp))
plot([0 mean(sin(signals.gamma_amp))],[0 mean(cos(signals.gamma_amp))],'r','linewidth',5)
xlabel('Real (cos)');xlim([-1.2 1.2])
ylabel('Imag (sin)');ylim([-1.2 1.2])
title('Modulated amplitude in z-plane')


%% Gamma amp by theta phase
subplot(3,4,11)
[~,bins] = hist(signals.theta_phase,20);
counts = zeros(1,20);
for idx = 1:19
    inds = bins(idx) < signals.theta_phase & bins(idx+1) > signals.theta_phase;
    counts(idx) = mean(signals.gamma_amp(inds));
end
bar(bins,counts,pi/2);
xlim([-pi pi]);
xlabel('Modulating Phase');
ylabel('Modulated Power');
title('Phase-Amplitude Distribution');


%% Phase * Amp on unit circle
subplot(3,4,12);hold on;grid on;
canolty = signals.gamma_amp .* exp(1i * signals.theta_phase);
plot(real(canolty),imag(canolty))
plot([0 mean(real(canolty))],[0 mean(imag(canolty))],'r','linewidth',3)
xlabel('Real');
ylabel('Imag');
title('Canolty MI')

%% Make more detailed phase plot

% Undo stacking if we have stacked
if stacked == 1
    signals = cfc_util_stacktrials(signals,signals.window_len*1000,'unstack');
end

%%
figure('position',[1100 100 512, 768]);

% Plot the phase, separate lines for each epoch
subplot(311);hold on;grid on
plot(squeeze(signals.time_vect),squeeze(signals.theta_phase))
ylim([-pi pi]);xlim([signals.time_vect(1) signals.time_vect(end)])
xlabel('Time (seconds')
ylabel('Phase (rads)')

% Plot the mean phase, should be (at least close to) zero for al time points
% TODO: perhaps could smooth this to reduce noise/paranoia
subplot(312);hold on; grid on;
plot([signals.time_vect(1) signals.time_vect(end)],[0 0],'k','linewidth',2);
plot(signals.time_vect,squeeze(mean(signals.theta_phase,3)))
ylim([-pi pi]);xlim([signals.time_vect(1) signals.time_vect(end)]);
xlabel('Time (seconds)');
ylabel('Average phase');

% Plot the phase slips with an approximate threshold based on the median of the differential
subplot(313);hold on; grid on
slipdat = diff(unwrap(signals.theta_phase));
medslip = mean(median(slipdat));
medmedslip = mean(std(slipdat));
plot(signals.time_vect(2:end),squeeze(slipdat))
plot([signals.time_vect(1) signals.time_vect(end)],[medslip + 3*medmedslip,medslip + 3*medmedslip],'r--')
plot([signals.time_vect(1) signals.time_vect(end)],[medslip - 3*medmedslip,medslip - 3*medmedslip],'r--')
legend({'Phase differential','Approx upper bound','Approx lower bound'});
xlabel('Time (seconds)');
ylabel('Phase progression');
xlim([signals.time_vect(1) signals.time_vect(end)])
