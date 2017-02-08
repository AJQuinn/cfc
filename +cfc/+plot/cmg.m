function cmg(cmg,metric_name,spectrum,freq_vect,log_high,col_levels)

if nargin < 2 || isempty(metric_name)
    if isfield(cmg,'mi_norm')
        metric_name = 'mi_norm';
    elseif isfield(cmg,'mi')
        metric_name = 'mi';
    elseif isfield(cmg,'esc')
        metric_name = 'ecs';
    else
        error('Please specify which metric to plot!')
    end
end

metric = cmg.(metric_name);

if ndims(metric) >= 3
    metric = squeeze(mean(metric,4));
end

if nargin <= 2

    %% Just plot the cmg without any fancy spectrum sidebars
    %col_levels = linspace(min(min(metric)),max(max(metric)),24);
    col_levels = linspace(max(max(metric))/3,max(max(metric)),24);
    %col_levels = 0:.25/24:.25;
    figure;
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        squeeze(metric)',36,...
        'linestyle','none');
    caxis([col_levels(1) col_levels(end)]);
    colormap(flipud(colormap('hot')));
    xlabel('Modulating Frequency (Hz)');
    ylabel('Modulated Frequency (Hz)');
    set(gca,'Color',[0.8 0.8 0.8]);

    colorbar();
else

    if nargin < 6 || isempty(col_levels)
        % Set colour levels based on the data
        col_levels = linspace(max(max(metric))/9,max(max(metric)),24);
        %col_levels = 0:.25/36:.25;
    end

    if nargin < 5 || isempty(log_high)
        % Default to not using a log scale for the high frequency
        log_high = false;
    end

    %% Extract the relevant parts of the spectrum
    lo_start = mean(cmg.lo_freqs(:,1));
    lo_stop = mean(cmg.lo_freqs(:,end));
    hi_start = mean(cmg.hi_freqs(:,1,1));
    hi_stop = mean(cmg.hi_freqs(:,end,end));

    lo_freq = freq_vect(lo_start < freq_vect & freq_vect < lo_stop);
    lo_pow = spectrum(lo_start < freq_vect & freq_vect < lo_stop);

    hi_freq = freq_vect(hi_start < freq_vect & freq_vect < hi_stop);
    hi_pow = spectrum(hi_start < freq_vect & freq_vect < hi_stop);

    %% Plot modulating frequency underneath the main comodulogram
    figure('position',[100 100 1024 768]);
    axes('position',[.3,.1,.6,.2]);
    plot(lo_freq,lo_pow)
    xlim([lo_start lo_stop]);
    xlabel('Modulating Frequency (Hz)');
    set(gca,'FontSize',24);grid on;

    %% Plot the modulated frequency spectrum vertically to the left of the
    %  comodulogram
    axes('position',[.1,.3,.2,.6]);
    plot(hi_freq,hi_pow)
    ylim([0 max(max(hi_pow))*1.2])
    view(-270, 90);
    set(gca, 'xdir', 'reverse');
    xlim([hi_start hi_stop]);
    xlabel('Modulated Frequency (Hz)');
    % Use log scale if requested
    if log_high == true; set(gca, 'YScale','log'); end
    set(gca,'FontSize',24);grid on;


    %% Plot the main comodulogram
    axes('position',[.3,.3,.6,.6])
    contourf(squeeze(mean(cmg.lo_freqs,1)),...
        squeeze(mean(cmg.hi_freqs(:,:,1),1)),...
        squeeze(metric)',18,...
        'linestyle','none');
    caxis([col_levels(1) col_levels(end)]);
    colormap(flipud(colormap('hot')));
    colorbar('position',[.92,.3,.03,.6]);
    set(gca,'XTickLabel','');set(gca,'YTickLabel','')
    set(gca,'Color',[0.8 0.8 0.8]);
    grid on
end

set(gca,'FontSize',24)
