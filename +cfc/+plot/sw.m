function sw( cfc_results, metric, outpath )

if nargin < 3 || isempty( outpath )
    outpath = nan;
end

if nargin < 2 || isempty( metric )
    % TODO: This option is currently broken, won't make any difference
    metric = 'MI';
end


time = cfc_results.time_vect;
thresh = nan;

for met_idx = 1:length(cfc_results.cfg.metrics)
    if strcmp(cfc_results.cfg.metrics{met_idx},'ESC')
        metric_name = 'ESC';
        metric = cfc_results.esc;
        if isfield(cfc_results,'esc_null'); metric_nulls = cfc_results.esc_null;end
        if isfield(cfc_results,'esc_thresh');thresh = cfc_results.esc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'NESC')
        metric_name = 'NESC';
        metric = cfc_results.nesc;
        if isfield(cfc_results,'nesc_null'); metric_nulls = cfc_results.nesc_null;end
        if isfield(cfc_results,'nesc_thresh');thresh = cfc_results.nesc_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'AEC')
        metric_name = 'AEC';
        metric = cfc_results.aec;
        if isfield(cfc_results,'aec_null'); metric_nulls = cfc_results.aec_null;end
        if isfield(cfc_results,'aec_thresh');thresh = cfc_results.aec_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'PLV')
        metric_name = 'PLV';
        metric = cfc_results.plv;
        if isfield(cfc_results,'plv_null'); metric_nulls = cfc_results.plv_null;end
        if isfield(cfc_results,'plv_thresh');thresh = cfc_results.plv_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'GLM')
        metric_name = 'GLM';
        metric = cfc_results.glm;
        if isfield(cfc_results,'glm_null'); metric_nulls = cfc_results.glm_null;end
        if isfield(cfc_results,'glm_thresh');thresh = cfc_results.glm_thresh;end
    elseif strcmp(cfc_results.cfg.metrics{met_idx},'MI')
        metric_name = 'MI';
        metric = cfc_results.mi;
        if isfield(cfc_results,'mi_null'); metric_nulls = cfc_results.mi_null;end
        if isfield(cfc_results,'mi_thresh');thresh = cfc_results.mi_thresh;end
    else
        %fprintf('CFC Metric %s not recognised!\nPlease choose from:\nESC, NESC, AEC, PLV, GLM and MI',cfc_results.cfg.metrics{met_idx});
    end

    % Sliding windows
    figure;
    if exist('metric_nulls') == 1
        subplot(211)
    else
        subplot(111)
    end

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

    if ~isnan(outpath)
        saveas(gcf,[outpath '_sw_timecourse'],'epsc')
    end

    if exist('metric_nulls') == 1
        % Histograms
        bins = 0:max(metric)/24:max(metric);
        subplot(223)
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

        subplot(224)
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

        if ~isnan(outpath)
            saveas(gcf,[outpath '_sw_histograms'],'epsc')
        end

    end


end

