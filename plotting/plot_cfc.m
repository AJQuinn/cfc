function plot_cfc( cfc_results )

time = mean(cfc_results.time_vect,2);
thresh = nan;

for met_idx = 1:length(cfc_results.cfg.metrics)
    if strcmp(cfc_results.cfg.metrics{met_idx},'ESC')
        metric_name = 'ESC';
        metric = cfc_results.esc;
        if isfield(cfc_results,'esc_thresh');thresh = cfc_results.esc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'NESC')
        metric_name = 'NESC';
        metric = cfc_results.nesc;
        if isfield(cfc_results,'nesc_thresh');thresh = cfc_results.nesc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'AEC')
        metric_name = 'AEC';
        metric = cfc_results.aec;
        if isfield(cfc_results,'aec_thresh');thresh = cfc_results.aec_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'PLV')
        metric_name = 'PLV';
        metric = cfc_results.plv;
        if isfield(cfc_results,'plv_thresh');thresh = cfc_results.plv_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'GLM')
        metric_name = 'GLM';
        metric = cfc_results.glm;
        if isfield(cfc_results,'glm_thresh');thresh = cfc_results.glm_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'MI')
        metric_name = 'MI';
        metric = cfc_results.mi;
        if isfield(cfc_results,'mi_thresh');thresh = cfc_results.mi_thresh;end
    else
        fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfg.metrics{met_idx});
    end

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

    thresh = nan;

end

