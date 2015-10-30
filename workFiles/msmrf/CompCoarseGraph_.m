function [Uc, Ec, Pc, vCoarse, plongTab] ...
    = CompCoarseGraph_(U,E,P,L,vEdgeList,gP)

%
% CompCoarseGraph(U,E,P,L,vEdgeList,gp) -
% compute coarse representation of the graphical model by contracting edges
%
% OUPUT:
%
%   Uc, Ec, Pc  -   unary, adjacency and pairwise of coarse graph
%
%   vCoarse     -   mapping from coarse vertices to fine vertices
%                   see definition below
%
%   plongTab    -   prolongation table, [NC]x[K+2] matrix          
%                   NC - number of coarse vertices
%                   plongTab(:,1) - fine-index of the interpolant
%                   plongTab(:,2) - coarse-index of the interpolator
%                   plongTab(:,2+k) - interpolation rule (1 <= k <= K)
%                       e.g. if the coarse-vertex 'i' has label 'k' then
%                       fine-vertex 'j' gets label plongTab(i,j,2+k)
%


M = size(E,1);              % number of edges in vEdgeList
N = size(U,1);              % number of vertices in the fine graph
K = size(U,2);              % number of labels

vFine = zeros(N,1);         % table for mapping fine vertices to coarse
                            % vertices; if 'i' is an index of a fine vertex
                            % then abs(vFine(i)) is its coarse counterpart
                            % - if - 
                            % vFine(i) > 0  --> 'i' is interpolator
                            % vFine(i) < 0  --> 'i' is interpolant

vCoarse = zeros(N,1);       % table for mapping coarse vertices to fine
                            % vertices; if 'I' is an index of a coarse
                            % vertex then vCoarse(I) is its corresponding
                            % fine vertex
                            
                            
plongTab = zeros(N,K+2);    % prolongation table
                            % plongTab(j,i,vec) mean that the coarse vertex
                            % indexed by 'i' interpolates fine-vertex 'j'
                            % according to the interpolation rule of 'vec'

vSkpList = ones(M,1);       % list of skipped edges; edges that were not
                            % considered in the coarse graph. Need to re-
                            % iterate over those
                            
cntrV = 1;                  % index counter for coarse [V]ertices

Uc = zeros(N,K);


if gP.wDegEnt
    % the unary term of each vertex is normalized by
    % the vertex' degree, compute the degree for all the vertices
      
    N = size(U,1);                      % number of vertices
    vDeg = accumarray(E(:,1),1,[N,1]);  % number of edges for which the
                                        % vertex is on the left side
                                        
    vDeg = vDeg + accumarray(E(:,2),1,[N,1]);   % degree of each vertex
    U_ = bsxfun(@rdivide,U,vDeg);
end


%
% iterate over edge list
% -- building the coarse unary term
for iM = 1
    
    pp = vEdgeList(iM);     % (signed) index of edge (ii,jj) in P
                            % if pp > 0 then ii is the 'left' vertex in P
                            % if pp < 0 then jj is the 'left' vertex in P
    ii = E(abs(pp),1);      % interpolator of edge iM
    jj = E(abs(pp),2);      % interpolant of edge iM
    if (pp < 0)
        % the direction of interpolation is reversed
        
        tmp = ii;
        ii = jj;
        jj = tmp;
    end
    
    
    %
    % check 'status' of the iM-th edge
    if (vFine(ii) == 0) && (vFine(jj) == 0)
        % both endpoints are free, edge can be inserted
        
        % update index table
        vFine(ii) = cntrV;          % vertex vFine(ii) is the interpolator
                                    % of coarse vertex cntrV
        vFine(jj) = -1 * cntrV;     % vertex vFine(jj) is the interpolant
                                    % of coarse vertex cntrV
        vCoarse(cntrV) = ii;        % coarse vertex cntrV maps to fine
                                    % vertex ii
        
        % find the interpolation assignment by solving:
        % l_j = argmin(@l) { \phi_{i,j} (l_i,l)  +  \phi_j (l) }
        % ..prepare the functional above
        pairwise = squeeze(P(:,:,abs(pp)));
        if (pp > 0)                 % this *if* makes sure that
            pairwise = pairwise';   % jj is on 1st dim, ii on 2nd dim
        end
        pairwise = bsxfun(@plus,pairwise,U(jj,:)'); % CHANGE U_ -> U
        [v, idx] = min(pairwise,[],1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pairwise = squeeze(P(:,:,abs(pp)));
        if (pp > 0)                 % this *if* makes sure that
            pairwise = pairwise';   % jj is on 1st dim, ii on 2nd dim
        end
        pairwise = bsxfun(@plus,pairwise,U(jj,:)');
        v = pairwise(sub2ind(size(pairwise),idx,1:K));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(L)
            % labels are initialized, overrule values
            % in the prolongation table
            v(L(ii)) = pairwise(L(jj),L(ii));
            idx(L(ii)) = L(jj);
        end
        v = v + U(ii,:);      % unary term of the new coarse vertex
        
        % update the coarse unary term
        Uc(cntrV,:) = v;
        cntrV = cntrV + 1;
        
        % update prolongation table
        plongTab(jj,:) = [jj,vFine(ii),idx];
        
        % update skipped list
        vSkpList(abs(pp)) = 0;
        
    elseif gP.bigAgg && (vFine(ii) > 0) && (vFine(jj) == 0)
        % edge starts at an interpolator vertex, and ends at a free vertex
        
        % update index table
        vFine(jj) = -1 * vFine(ii);
        
        % find the interpolation assignment
        pairwise = squeeze(P(:,:,abs(pp)));
        if (pp > 0)                 % this *if* makes sure that
            pairwise = pairwise';   % jj is on 1st dim, ii on 2nd dim
        end
        pairwise = bsxfun(@plus,pairwise,U(jj,:)'); % CHANGE U_ -> U
        [v, idx] = min(pairwise,[],1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pairwise = squeeze(P(:,:,abs(pp)));
        if (pp > 0)                 % this *if* makes sure that
            pairwise = pairwise';   % jj is on 1st dim, ii on 2nd dim
        end
        pairwise = bsxfun(@plus,pairwise,U(jj,:)');
        v = pairwise(sub2ind(size(pairwise),idx,1:K));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(L)
            % labels are initialized, overrule values
            % in the prolongation table
            v(L(ii)) = pairwise(L(jj),L(ii));
            idx(L(ii)) = L(jj);
        end

        % update the coarse unary term
        Uc(vFine(ii),:) = Uc(vFine(ii),:) + v;
        
        % update prolongation table
        plongTab(jj,:) = [jj,vFine(ii),idx];
        
        % update skipped list
        vSkpList(abs(pp)) = 0;
               
    end
    
end


%
% consider vertices that were not accounted for
% by the previous loop
% -- complete the construction of the coarse unary term
idx = find(vFine == 0);             % indices of unaccnt'd fine vertices
cntrV_ = cntrV + length(idx) - 1;   % update counter of coarse vertices
vCoarse(cntrV:cntrV_) = idx;        % update vCoarse
Uc(cntrV:cntrV_,:) = U(idx,:);      % update Uc
vFine(idx) = cntrV : cntrV_;        % update vFine
vCoarse = vCoarse(1:cntrV_);        % trim excess 0's
Uc = Uc(1:cntrV_,:);                % trim excess 0's

vSkpList = find(vSkpList);
Mc = length(vSkpList);              % max number of edges in coarse graph
Pc = zeros(K,K,Mc);                 % coarse pairwise
Ec = zeros(Mc,2);                   % coarse adjacency matrix

mAdj = sparse([],[],[],...
    cntrV_,cntrV_,Mc);              % adjacency matrix of coarse graph

cntrE = 1;                          % index counter for coarse [E]dges

%
% iterate over skipped edges
% -- building the coarse pairwise/connectivity
for iS = 1 : Mc
    
    ii = E(vSkpList(iS),1); % left vertex of edge iS wrt fine graph
    jj = E(vSkpList(iS),2);	% right vertex of edge iS wrt fine graph

    ii_ = vFine(ii);        % index of left vertex wrt coarse graph
    jj_ = vFine(jj);        % index of right vertex wrt coarse graph
    
    % interpolated vertices are assigned with a
    % label according to their interpolators!
    if (ii_ < 0)
        % vertex ii_ is interpolated
        iiperm = plongTab(ii,3:end);
        ii_ = abs(ii_);
    else
        iiperm = 1 : K;       
    end  
    if (jj_ < 0)
        % vertex jj_ is interpolated
        jjperm = plongTab(jj,3:end);
        jj_ = abs(jj_);
    else
        jjperm = 1 : K;
    end
    
    % process the pairwise
    % such that interpolated labels are considered correctly
    pairwise = P(:,:,vSkpList(iS));
    pairwise = pairwise(iiperm,jjperm);
    
    % check if edge is a self-loop
    % --this can happen in big aggregator mode
    if (ii_ == jj_)
        % update the relevant unary term
        
        Uc(ii_,:) = Uc(ii_,:) + diag(pairwise)';
        continue                % stop processing the edge
    end
        
    % check if edge exists in coarse graph
    idx = mAdj(ii_,jj_);
    if (idx == 0)
        % edge has not been represented on coarse graph
        idx = cntrE;
        mAdj(ii_,jj_) = idx;
        mAdj(jj_,ii_) = -idx;   % flag that edge is represented in reverse
        Ec(idx,:) = [ii_, jj_];
        cntrE = cntrE + 1;
    end
    
    % update pairwise in Pc
    if (idx < 0)
        % the edge is represented in reverse
        
        pairwise = pairwise';
        idx = abs(idx);
    end
    Pc(:,:,idx) = Pc(:,:,idx) + pairwise;
    
end

Ec = Ec(1:cntrE-1,:);
Pc = Pc(:,:,1:cntrE-1);

% trim excess 0's in plongTab - 
% it wasn't trimmed until now in order to allow
% fast indexing in lines ~225-235
plongTab = plongTab(vFine < 0,:);


end