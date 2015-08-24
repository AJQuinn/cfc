function [weighted_signals] = cfc_util_weightsignal(signals,Gamma)

nstates = size(Gamma,2);

weighted_signals = cell(nstates);

for idx = 1:nstates
    
    % Weight signals
    skip_field = {'hi_bounds','time_vect','hi_steps',...
        'lo_bounds','lo_bandwidth','lo_steps','sr'};
    fields = fieldnames(signals);
    for i = 1:numel(fields)
        if strmatch(fields{i},skip_field)
            tmp.(fields{i}) = signals.(fields{i});
        else
            tmp.(fields{i}) = signals.(fields{i}).*Gamma(:,idx);
        end
    end
    weighted_signals{idx}.signals = tmp;
end