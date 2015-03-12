function perm_signals = permute_signal(signals, method)
% Function for scrambling a cfc signals object for the generation of null
% distributions

perm_signals = signals;
skip_field = {'hi_bounds','lo_bounds','sr',...
              'time_vect','true_timecourse','signal'};
fields = fieldnames(signals);



if strcmp(method,'window')
    % We want to mix up the windows
    nwindows = size(signals.theta,2);
    
    for idx = 1:numel(fields)
        if strmatch(fields{idx},skip_field)
            continue
        else
            perm_idx = randperm(nwindows);
            perm_signals.(fields{idx}) = signals.(fields{idx})(:,perm_idx);
        end
    end
    
end