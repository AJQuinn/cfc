function plot_group_peaks( group_peaks )
%% Function to plot the distribution of peaks across a group

    nppts = length(group_peaks);
    freq_vect = group_peaks{1}.freq_vect;
    
    % Get smoothed fft from each participant
    % not currently sanity-checking the smoothing parameters....
    ffts = arrayfun(@(i) group_peaks{i}.smo_fft, 1:nppts,'UniformOutput',false);
    group_fft = cell2mat(ffts');
    
    % Get top peak from each participant
    peak_freqs = arrayfun( @(i) group_peaks{1}.peak_frequencies(group_peaks{i}.peaks_by_amplitude(1)), 1:nppts);
    peak_amps = arrayfun( @(i) group_peaks{1}.peak_amplitudes(group_peaks{i}.peaks_by_amplitude(1)), 1:nppts);
    
    
    %% Start plotting
    figure('Position',[1, 1, 800, 500]);
  
    % Average spectrum and individuals
    axes('position',[.3,.3,.6,.6])
    plot(freq_vect,mean(group_fft)); hold on
    grid on;
    amp_range = ylim;
    xlim([freq_vect(1) freq_vect(end)]);
    plot(peak_freqs,peak_amps,'r*');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);

    legend({'Average Spectrum','Individual Peak'});
    
    
    % horizontal histogram
    h_bins = freq_vect(1):freq_vect(end)/24:freq_vect(end);
    axes('position',[.3,.1,.6,.2])
    [n1,xout1] = hist(peak_freqs,h_bins);
    bar(xout1,n1,1,'b');hold on;grid on;
    xlim([freq_vect(1) freq_vect(end)]);
    
    xlabel('Frequency (Hz)');
    
    % vertical histogram
    v_bins = amp_range(1):amp_range(end)/24:amp_range(end);
    axes('position',[.1,.3,.2,.6])
    [n1, xout1] = hist(peak_amps,v_bins);
    barh(xout1,n1,1,'b');hold on;grid on;
    ylim(amp_range);
    ylabel('Amplitude');

end