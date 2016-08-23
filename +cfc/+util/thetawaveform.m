function [ waveform_phase ] = thetawaveform( signal, theta_freq, sr )
%% Brute force estimation of theta phase by extraction of maxima, minima and
% zero crossings. This avoids bandwidth problems.

    min_dist = floor(.9*sr / (theta_freq));

    signal = sgolayfilt(signal,3,2.*round((min_dist+1)/2)-1);
    
    % Find local maxima
    [max_pks,max_locs] = findpeaks(signal,'MINPEAKDISTANCE', min_dist);
    % Find local minima
    [min_pks,min_locs] = findpeaks(-signal,'MINPEAKDISTANCE', min_dist);

    % Find ascending zero crossings
    [asz_pks,asz_locs] = findpeaks(diff(sign(signal)),'MINPEAKDISTANCE', min_dist);
    % Find descending zero crossings
    [dsz_pks,dsz_locs] = findpeaks(diff(-sign(signal)),'MINPEAKDISTANCE', min_dist);

    waveform_phase = zeros(size(signal,1),size(signal,2));

    % Interpolate partial first cycle if necessary
    if min_locs(1) < max_locs(1) && ...
            min_locs(1) < asz_locs(1) && ...
            min_locs(1) < dsz_locs(1)
        % We start with a minima
        min_idx = 1;
        disp('minima');

    elseif max_locs(1) < min_locs(1) && ...
            max_locs(1) < asz_locs(1) && ...
            max_locs(1) < dsz_locs(1)
        % We start with a maxima
        x = [max_locs(1) dsz_locs(1) min_locs(1)];
        y = [pi          (3*pi)/2    2*pi];
        min_idx = 2;
        disp('maxima');
        waveform_phase(max_locs(1):min_locs(1)) = interp1(x,y,max_locs(1):min_locs(1));

    elseif asz_locs(1) < min_locs(1) && ...
            asz_locs(1) < max_locs(1) && ...
            asz_locs(1) < dsz_locs(1)
        x = [asz_locs(1) max_locs(1) dsz_locs(1) min_locs(1)];
        y = [pi/2        pi          (3*pi)/2    2*pi];
        min_idx = 2;
        % We start ascending
        waveform_phase(asz_locs(1):min_locs(1)) = interp1(x,y,asz_locs(1):min_locs(1));

        disp('ascending')
    elseif dsz_locs(1) < min_locs(1) && ...
            dsz_locs(1) < max_locs(1) && ...
            dsz_locs(1) < asz_locs(1)
        % We start descending
        x = [dsz_locs(1) min_locs(1)];
        y = [(3*pi)/2    2*pi];
        min_idx = 2;
        waveform_phase(dsz_locs(1):min_locs(1)) = interp1(x,y,dsz_locs(1):min_locs(1));
        disp('decending')
    end

    % Interpolate for each individual oscillation
    size(min_locs,1)
    for idx = min_idx:size(min_locs,1)-1
%         asz_idx = find(asz_locs > min_locs(idx),1);
%         max_idx = find(max_locs > min_locs(idx),1);
%         dsz_idx = find(dsz_locs > min_locs(idx),1);
%         min2_idx = find(min_locs > min_locs(idx),1);
        asz_idx = find(asz_locs > min_locs(idx),1);
        max_idx = find(max_locs > asz_idx,1);
        dsz_idx = find(dsz_locs > max_idx,1);
        min2_idx = find(min_locs > dsz_idx,1);

        %oscillation = waveform_phase(min_locs(idx):min_locs(idx+1));
        x = [min_locs(idx) asz_locs(asz_idx) max_locs(max_idx) dsz_locs(dsz_idx) min_locs(min2_idx)]
        y = [0 pi/2 pi (3*pi)/2 2*pi];
        try
            waveform_phase(min_locs(idx):min_locs(min2_idx)) = interp1(x,y,min_locs(idx):min_locs(min2_idx));
        catch
            % We might have a duplicate value
            [~, ind] = unique(x);
            x(setdiff(1:5, ind)) = x(setdiff(1:5, ind))+1;
            waveform_phase(min_locs(idx):min_locs(min2_idx)) = interp1(x,y,min_locs(idx):min_locs(min2_idx));
        end
           
        min_locs(idx)
        
    end

    waveform_phase = waveform_phase - pi;

end

