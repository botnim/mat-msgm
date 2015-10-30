function [L, prior, Q] = Solve(Uc,Ec,Pc,Lc,priorc,...
                            U,E,P,L,prior,...
                            vCoarse,plongTab,gP)

%
% Solve(Uc,Ec,Pc,Lc,priorc) - solve the graphical model, either in exact or
% approximately
%


%
% apply appropriate method -
% depending on the number of variables
if (size(Uc,1) == 1)
    % only 1 vertex, winner takes all
    
    [~, Lc] = min(Uc,[],2);
    
elseif (size(Uc,1) == 2)
    % 2 vertices remaining, find optimal assignment
    
    K = size(Uc,2);     % number of labels
    pairwise = Pc;
    ii = Ec(1);         % left vertex of pairwise
    jj = Ec(2);         % right vertex of pairwise
    
    % find optimal assignment
    minEng = Inf;
    Lc = zeros(2,1);
    for i = 1 : K
        for j = 1 : K    
            
            e = Uc(ii,i) + Uc(jj,j) + pairwise(i,j);
            if (e <= minEng)
                minEng = e;
                Lc([ii,jj]) = [i;j];
            end
        end
    end
    
else
    % graphical model has > 3 vertices,
    % provide an approximate solution
    
    if any(Lc)
        % there's an initial labeling
        Lc = RelaxGraph(Uc,Ec,Pc,Lc,gP);
        
    else
        % no initial labeling
        [~, Lc] = min(Uc,[],2);        % initialize - winner takes all
        Lc = RelaxGraph(Uc,Ec,Pc,Lc,gP);
    end
    
    priorc = UpdatePrior(priorc,Lc,gP);
    
end



%
% interpolate the solution
L = Prolong(Lc,vCoarse,plongTab);
if any(priorc)
    % update the prior, where relevant
    prior(vCoarse,:) = priorc;
end

%
% assertion - remove later
EnergyAssert(U,E,P,L,Uc,Ec,Pc,Lc);

%
% relax the graph and update the prior
L = RelaxGraph(U,E,P,L,gP);
prior = UpdatePrior(prior,L,gP);

Q.e = [Energy(Uc,Ec,Pc,Lc), Energy(U,E,P,L)];
Q.nv = [size(Uc,1), size(U,1)];

