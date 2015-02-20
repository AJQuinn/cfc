function trans_results = transform_cfc(results)

trans_results = results;
trans_results.tranformed = 1;

if isfield(results,'esc')
    % Fisher's z transformation
    trans_results.esc = 0.5*log((1+results.esc)./(1-results.esc));
end
if isfield(results,'plv')
    % Arcsine transform
    trans_results.plv = asin( 2*results.plv - 1);
end
if isfield(results,'mi')
    % Log transform
    trans_results.mi = log(results.mi);
end
if isfield(results,'glm')
    % Fisher's z transformation
    trans_results.glm = 0.5*log((1+results.glm)./(1-results.glm));
end
