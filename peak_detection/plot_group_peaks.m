function plot_group_peaks( varargin )
%% Function to plot the distribution of peaks across a group
%
% input is a cell array containing quinnpeaks output structs for one
% participant in each cell. Multiple (max of 4) cell arrays can be passed
% in in order to display different groups or conditions in one plot.
%
% plot_group_peaks( group_peak_data );
%
% plot_group_peaks( condition_1_peaks, condition_2_peaks );
%

    nconditions = length(varargin);

    colour_cycle = {'b','r','g','c'};
    if nconditions > 3;
        error('That is a lot of conditions, perhaps too many...');
    end

    %% Get the figure started
    figure('Position',[1, 1, 800, 500]);

    % Average spectrum and individuals
    main_axis = axes('position',[.3,.3,.6,.6]);
    horz_axis = axes('position',[.3,.1,.6,.2]);
    vert_axis = axes('position',[.1,.3,.2,.6]);

    %% Main loop
    for cond_idx = 1:nconditions
        group_peaks = varargin{cond_idx};

        nppts = length(group_peaks);
        nbins = fix(sqrt(nppts));
        if nbins < 5; nbins = 5; end;

        freq_vect = group_peaks{1}.freq_vect;

        % Get smoothed fft from each participant
        % not currently sanity-checking the smoothing parameters....
        ffts = arrayfun(@(i) group_peaks{i}.smo_fft, 1:nppts,'UniformOutput',false);
        group_fft = cell2mat(ffts');

        % Get top peak from each participant
        peak_freqs = arrayfun( @(i) group_peaks{i}.peak_frequencies(group_peaks{i}.peaks_by_amplitude(1)), 1:nppts);
        peak_amps = arrayfun( @(i) group_peaks{i}.peak_amplitudes(group_peaks{i}.peaks_by_amplitude(1)), 1:nppts);


        %% Start plotting
        % Average spectrum and individuals

        plot(main_axis,freq_vect,mean(group_fft,1),colour_cycle{cond_idx});
        hold(main_axis,'on');
        grid(main_axis,'on');
        plot(main_axis,peak_freqs,peak_amps,[colour_cycle{cond_idx} '*']);

        % horizontal histogram
        [n1,xout1] = hist(peak_freqs,nbins);
        bar(horz_axis,xout1,n1,1,colour_cycle{cond_idx});
        hold(horz_axis,'on');
        grid(horz_axis,'on');

        % vertical histogram
        [n1, xout1] = hist(peak_amps,nbins);
        barh(vert_axis,xout1,n1,1,colour_cycle{cond_idx});
        hold(vert_axis,'on');
        grid(vert_axis,'on');

        if cond_idx == 1;
            leg = {'Condition 1 Average','Condition 1 Individual'};
        else
            leg = cat(2,leg,{['Condition ' num2str(cond_idx) ' Average'], ...
                            ['Condition ' num2str(cond_idx) ' Individual']});
        end

    end

    % Main axis formatting
    legend(main_axis,leg);
    amp_range = ylim(main_axis);
    xlim(main_axis,[freq_vect(1) freq_vect(end)]);
    set(main_axis,'YTickLabel',[]);
    set(main_axis,'XTickLabel',[]);

    % Horizontal axis formatting
    xlabel(horz_axis,'Frequency (Hz)');
    xlim(horz_axis,[freq_vect(1) freq_vect(end)]);
    ylabel(horz_axis,'Frequency');

    % Vertical axis formatting
    ylim(vert_axis,amp_range);
    ylabel(vert_axis,'Amplitude');
    xlabel(vert_axis,'Frequency');

end
