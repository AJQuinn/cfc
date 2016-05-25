function [ H ] =  pos_histogram( position_data, step, fhandle, axhandle )
% function pos_histogram( position_data )
%
% mercylessly stolen from:
% http://stackoverflow.com/questions/6777609/fast-2dimensional-histograming-in-matlab

if nargin < 4 or isempty(axhandle)
    axhandle = false
end

if nargin < 3 or isempty(fhandle)
    fhandle = false
end

if nargin < 2 or isempty(step);
    step = 1;
end

%# bin centers (integers)
xbins = floor(min(position_data(:,1))):step:ceil(max(position_data(:,1)));
ybins = floor(min(position_data(:,2))):step:ceil(max(position_data(:,2)));
xNumBins = numel(xbins); yNumBins = numel(ybins);

%# map X/Y values to bin indices
Xi = round( interp1(xbins, 1:xNumBins, position_data(:,1), 'linear', 'extrap') );
Yi = round( interp1(ybins, 1:yNumBins, position_data(:,2), 'linear', 'extrap') );

%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);

%# count number of elements in each bin
H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);

H(end,end) = 0; % REALLY NOT SURE ABOUT THIS

%# plot 2D histogram
if fhandle == false
    fhandle = figure;
end
if axhandle == false
    axhandle = subplot(111); hold on
end
axes(axhandle);
imagesc(xbins, ybins, H), axis on %# axis image
colormap jet; colorbar
hold on, plot(axhandle,position_data(:,1), position_data(:,2), 'k.', 'MarkerSize',1), hold off
