function prior = UpdatePrior(prior,L,gP)

%
% UpdatePrior(prior,L) - update the prior with information from the current
% labeling assignment 'L'
%
% 'prior' refers to the estimated (marignal) distribution of the variables,
% it is an NxK array, where prior(n,k) is the (estimated) probability that
% the variable 'n' has label 'k'.
%

if isempty(prior) || isempty(L)
    % running msmrf without prior estimation
    
    return
end


N = size(prior,1);          % number of variables
K = size(prior,2);          % number of labels (+1, in case of 'hist')

switch gP.prUpdate
    
    case 'rapid'
        % rapid adaptation of the prior,
        % prior is mostly affected by current labeling
        
        alpha = 1/(K+1);            % weight parameter for new information
        prior_ = sparse(1:N,L,alpha,N,K);   % new information
        prior = K/(K+1) * prior + prior_;   % update by weighted average
        
    case 'uniform'
        % maintain a histogram of 'evidence'
        % to turn the histogram into a prior, it must be normalized.
        % The reason that it is not normalized at each update is to keep
        % count of how much 'evidence' was observed. This is necessary
        % because different vertices 'observe' differnt amount of
        % 'evidence' and in this way each new evidence gets a weight which
        % is proportional to its part out of the total evidence observed by
        % the vertex.
        % This is in fact posterior predictive with uniform prior
        
        prior = prior + sparse(1:N,L,1,N,K);    % update 'evidence'
        prior(:,end) = prior(:,end) + 1;        % update count
        
    case 'jeff'
        % similar to 'hist',
        % This is in fact posterior predictive with Jeffreys prior
        
        prior = prior + sparse(1:N,L,1,N,K);    % update 'evidence'
        prior(:,end) = prior(:,end) + 1;        % update count
                      
end