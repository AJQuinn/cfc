function plot_comodulogram(cmg,spectrum,freq_vect)

if nargin < 2

    figure;
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        cmg.mi');
    colormap('hot');
    colorbar();
else

    lo_start = mean(cmg.lo_freqs(:,1));
    lo_stop = mean(cmg.lo_freqs(:,end));
    hi_start = mean(cmg.hi_freqs(:,1,1));
    hi_stop = mean(cmg.hi_freqs(:,end,end));

    lo_freq = freq_vect(lo_start < freq_vect & freq_vect < lo_stop);
    lo_pow = spectrum(lo_start < freq_vect & freq_vect < lo_stop);

    hi_freq = freq_vect(hi_start < freq_vect & freq_vect < hi_stop);
    hi_pow = spectrum(hi_start < freq_vect & freq_vect < hi_stop);

    figure;
    axes('position',[.3,.1,.6,.2])
    plot(lo_freq,lo_pow)
    xlim([lo_start lo_stop]);
    %ylim([0 max(lo_pow)*1.2]);
    xlabel('Modulating Frequency (Hz)');

    axes('position',[.1,.3,.2,.6])
    plot(hi_freq,hi_pow)
    ylim([0 max(max(hi_pow))*1.2])
    view(-270, 90);
    set(gca, 'xdir', 'reverse');
    xlim([hi_start hi_stop]);
    xlabel('Modulated Frequency (Hz)');

    axes('position',[.3,.3,.6,.6])
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        cmg.mi');
    colormap('hot');
    %colorbar();
    colorbar('position',[.92,.3,.03,.6]);
    %ax = axes('position',[.92,.3,.1,.6]);
    %set(gca,'visible','off');
    %colorbar('peer',ax);

end
