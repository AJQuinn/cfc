function [ hmm,Gamma,vpath ] = cfc_est_hmm(modulating_signal,modulated_signal,timelag,order,exptimelag)

if nargin < 5
    exptimelag = [];
end

if nargin < 4
    order = 20;
end

if nargin < 3
    timelag = 25;
end

signal = [modulating_signal, modulated_signal];

options = struct('K',2);
options.cyc = 100;
options.tol = 0.00001;
%options.inittype = 'GMM';
options.initcyc = 100;
options.initrep = 5;
options.order = order;
options.timelag = timelag;
options.embeddedlags = 0;
options.orderoffset = 25;
options.covtype = 'diag';
if ~isempty(exptimelag)
    options.exptimelag = exptimelag;
end

niters = 1;
minfe = 0;

for idx = 1:niters
    
    [hmm_tmp, Gamma_tmp, ~, vpath_tmp, ~, ~, fehist] = hmmmar (signal,size(signal,1),options);
    
    if minfe > fehist(end) || minfe == 0
        hmm = hmm_tmp;
        Gamma = Gamma_tmp;
        vpath = vpath_tmp;
    end
    
end

hmm.embedding = formorders(hmm.train.order,...
                hmm.train.orderoffset,...
                hmm.train.timelag,...
                hmm.train.exptimelag);
