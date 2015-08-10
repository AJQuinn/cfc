function plot_cfc( cfc_results )

time = mean(cfc_results.time_vect,2);
thresh = nan;

for met_idx = 1:length(cfc_results.cfg.metrics)
    if strcmp(cfc_results.cfg.metrics{met_idx},'ESC')
        metric_name = 'ESC';
        metric = cfc_results.esc;
        metric_nulls = cfc_results.esc_null;
        if isfield(cfc_results,'esc_thresh');thresh = cfc_results.esc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'NESC')
        metric_name = 'NESC';
        metric = cfc_results.nesc;
        metric_nulls = cfc_results.nesc_null;
        if isfield(cfc_results,'nesc_thresh');thresh = cfc_results.nesc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'AEC')
        metric_name = 'AEC';
        metric = cfc_results.aec;
        metric_nulls = cfc_results.aec_null;
        if isfield(cfc_results,'aec_thresh');thresh = cfc_results.aec_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'PLV')
        metric_name = 'PLV';
        metric = cfc_results.plv;
        metric_nulls = cfc_results.plv+null;
        if isfield(cfc_results,'plv_thresh');thresh = cfc_results.plv_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'GLM')
        metric_name = 'GLM';
        metric = cfc_results.glm;
        metric_nulls = cfc_results.glm_null;
        if isfield(cfc_results,'glm_thresh');thresh = cfc_results.glm_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'MI')
        metric_name = 'MI';
        metric = cfc_results.mi;
        metric_nulls = cfc_results.mi_null;
        if isfield(cfc_results,'mi_thresh');thresh = cfc_results.mi_thresh;end
    else
        fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
    end

    % Sliding windows
    figure;
    plot(time,metric)
    hold on
    title(sprintf('Cross frequency coupling between %dHz and %dHz',...
                        cfc_results.cfg.lo_freq, cfc_results.cfg.hi_freq));
    xlabel('Time (seconds)')
    ylabel(metric_name)

    if ~isnan(thresh)
        plot([time(1) time(end)], [thresh(3), thresh(3)])
        plot([time(1) time(end)], [thresh(2), thresh(2)])
        plot([time(1) time(end)], [thresh(1), thresh(1)])
        legend({metric_name,'Thresh .001','Thresh .01','Thresh .05'});
    else
        legend({metric_name});
    end

    % Histograms
    bins = 0:max(metric)/24:max(metric);
    figure
    subplot(211)
    [n1, xout1] = hist(metric_nulls,bins);
    bar(xout1,n1,'r'); grid; hold on;
    ylabel('Frequency');

    if ~isnan(thresh)
        ybounds = ylim;
        plot([thresh(3) thresh(3)],[ybounds(1) ybounds(2)]);
        plot([thresh(2) thresh(2)],[ybounds(1) ybounds(2)]);
        plot([thresh(1) thresh(1)],[ybounds(1) ybounds(2)]);
        legend({'Null Distribution','Thresh .001','Thresh .01','Thresh .05'});
    else
        legend({'Null Distribution'});
    end

    subplot(212)
    [n2, xout2] = hist(metric,bins);
    bar(xout2,n2,'g'); grid; hold on;

    if ~isnan(thresh)
        ybounds = ylim;
        plot([thresh(3) thresh(3)],[ybounds(1) ybounds(2)]);
        plot([thresh(2) thresh(2)],[ybounds(1) ybounds(2)]);
        plot([thresh(1) thresh(1)],[ybounds(1) ybounds(2)]);
        legend({'Observed Distribution','Thresh .001','Thresh .01','Thresh .05'});
    else
        legend({'Observed Distribution'});
    end

    xlabel(metric_name);
    ylabel('Frequency');
    title('Null and observed distributions');

    thresh = nan;

end

