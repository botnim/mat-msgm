function x = qpboms(U,E,P,L,gP)

u = U';
p = [E(:,1)'; E(:,2)'; ...
        squeeze(P(1,1,:))'; ...         	% Eij(0,0)
        squeeze(P(1,2,:))'; ...          	% Eij(0,1)
        squeeze(P(2,1,:))'; ...       	% Eij(1,0)
        squeeze(P(2,2,:))'];    

    
% initial solution
x = QPBO_wrapper_mex(u, p, int32(ones(size(U,1),1)), 'q');

if any(x == -1)
    
    % new graph
    [U_, E_, P_, ~, vCR2Fine, ~] = ...
        MakeCRgraph(U,E,P,x(x ~= -1)+1,find(x ~= -1), []);

    if isequal(U_, U)
        % solve and return

        % coarsen graph
        [U, P] = ReparamGraph(U,E,P,gP);
        vEdgeList = ScoreEdges(U,E,P,[],L,gP); 
        [Uc, Ec, Pc, vCoarse, plongTab] = CompCoarseGraph(U,E,P,L,vEdgeList,gP);
        Lc = Restrict(L,vCoarse);
        
        % recursive call
        xc = qpboms(Uc,Ec,Pc,Lc,gP);

        % interpolate
        x_ = Prolong(double(xc+1),vCoarse,plongTab);
        x_(x_ == 1) = 0;
        x_(x_ == 2) = 1;
        x(vCR2Fine) = x_;

    else
        % recursive call

        x_ = qpboms(U_,E_,P_,[],gP);
        x(vCR2Fine) = x_;
    end
end

    