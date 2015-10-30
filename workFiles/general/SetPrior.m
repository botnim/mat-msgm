function prior = SetPrior(N,K,gP)

switch gP.prUpdate
    
    case 'rapid'
        % use prior, with rapid updating
        prior = ones(N,K) / K;  % uniform prior
    
    case 'uniform'
        % prior in the form of a histogram,
        % corresponds to posterior predictive w. uniform
        % distribution of theta
        
        % last column of prior counts the amount of 'evidence'
        % for each variable, as it may be different in general
        prior = [ones(N,K), K * ones(N,1)];
        
    case 'jeff'
        % prior in the form of a histogram,
        % corresponds to posterior predictive w. Jeffreys prior
        prior = [ones(N,K) / K, ones(N,1)];
        
    case 'none'
        % do not use prior
        prior = [];      
end