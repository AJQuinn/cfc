function animcmg( obj, outpath )
% Make a movie from a sliding window comodulogram


ts_max = max(max(obj.signal));
ts_min = min(min(obj.signal));
if abs(ts_max) > abs(ts_min)
    ts_ylim = abs(ts_max) * 1.1;
else
    ts_ylim = abs(ts_min) * 1.1;
end


met_max = max(max(max(obj.mi_norm))) * .45;
%met_max = mean(mean(mean(obj.mi_norm))) + ( std(std(std(obj.mi_norm)))*2 );


figure;
%ax1 = axes('position',[.1,.1,.8,.25]);
ax1 = axes('position',[.25,.1,.5,.25]);
ax2 = axes('position',[.25,.43,.5,.55]);

nwindows = size(obj.mi_norm,1);

for idx = 1:nwindows

    %window_idx = obj.time_vect(idx,:) >= obj.time_vect(idx,:) & ...
    %    time_vect <= obj.time_vect(end,idx);

    plot(ax1,obj.time_vect(idx,:),obj.signal(idx,:));
    xlim(ax1,[min(obj.time_vect(idx,:)) max(obj.time_vect(idx,:))]);
    ylim(ax1,[-ts_ylim ts_ylim]);
    xlabel(ax1,'Time (seconds)');

    contourf(mean(obj.lo_freqs),mean(mean(obj.hi_freqs,3)),...
        squeeze(obj.mi_norm(idx,:,:))','LineStyle','none',...
        'parent',ax2);
    caxis([0 met_max]);
    colormap('hot');
    ax3 = colorbar('position',[.85,.43,.03,.55]);
    xlabel(ax2,'Modulating Frequency (Hz)');
    ylabel(ax2,'Modulated Frequency (Hz)');
    ylabel(ax3,'Modulation Index');

    print(sprintf('%s_%04d', outpath, idx),'-dpng');

end

cmd = sprintf('convert -delay 15 -loop 0 %s %s', ...
    sprintf('%s_*', outpath), ...
    [outpath '.gif']);

system(cmd);

end
