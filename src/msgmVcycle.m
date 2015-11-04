function [L, prior, Q] = msgmVcycle(U,E,P,L,prior,gP,fineGraph)

%
% Vcycle(D,E,W,L,gP)    -   A multiscale scheme for finding minimum energy
%                           assignments of pairwise MRFs
%
% Assume N variables, M edges, and K labels
%
% INPUT:
%
%   U   -   unary term, NxK matrix.
%           U(i,li) is the cost of assigning label 'li' to vertex 'i'
%
%
%   E   -   adjacency matrix, NxN matrix.
%           Let i,j denote two vertices; if the edge (i,j) exists, then
%           E(i,j) is an index pointing to the explicit representation
%           of their pairwise term in matrix 'P'. E(i,j) is zero otherwise.
%
%   P   -   explicit representation of all pairwise terms, KxKxM matrix.
%           If the edge connecting vertices (i,j) is indexed by 'k', then
%           P(li,lj,k) is the cost of assigning label 'li' to 'i' and
%           'lj' to 'j'. Note that in general P(y,z,:)!=P(z,y,:)
%
%   L   -   initial labeling of the vertices, Nx1 vector.
%           For an unlabeled graph, L should be the empty matrix [].
%
%   gP  -   set of parameters. See setParams().
%
%   fineGraph - (true/false) was the original fine graph passed as input to
%               msmrf(), or is it a coarsened version of it?
%               This information is used to save on computational resources
%               by not applying 'ReparmGraph()' on the finest graph more
%               than once.
%
%
% OUTPUT:
%
%   L   -   a labeling assignment to the vertices, Nx1 vector.
%


%
% reparameterize the unary and pairwise
% terms of the graphical model U,E,P
if (~fineGraph)
    % finest graph was already reparameterized

    [U, P] = ReparamGraph(U,E,P,gP);
end



%
% relax the graph
gP.bFineGraph = fineGraph;
[L, Q.e, Q.nv] = RelaxGraph(U,E,P,L,gP);
prior = UpdatePrior(prior,L,gP);


%
% score the edges, a preprocessing step
% for CompCoarseGraph()
vEdgeList = ScoreEdges(U,E,P,prior,L,gP);      

%
% compute coarse graphical model,
% according to 'vEdgeList'
[Uc, Ec, Pc, vCoarse, plongTab] = msgmCoarsening(U,E,P,L,vEdgeList,gP);
Lc = msgmInherit(L,vCoarse);
priorc = [];


%
% check stopping condition
% TODO: check if crsRatioThrs is necessary, if not - can bring to first
% step if V-cycle
% TODO: set size(Uc,1) <= 2
if (size(Uc,1)/size(U,1) >= gP.crsRatioThrs) || ...
        (size(Uc,1) <= 500)
    % coarsening ratio is above threshold, stop the recursion
   
    [L, prior, Q] = Solve(Uc,Ec,Pc,Lc,priorc, ...
        U,E,P,L,prior,...
        vCoarse,plongTab,gP);
    
    return
end


%
% recursive call
[Lc, priorc, Qc] = msgmVcycle(Uc,Ec,Pc,Lc,priorc,gP,false);
Q.e = [Q.e, Qc.e];
Q.nv = [Q.nv, Qc.nv];

%
% interpolate solution
if not(gP.bPlongCR)

    L = Prolong(Lc,vCoarse,plongTab);
else

    [Ucr, Ecr, Pcr, Lcr, vCR2Fine, L] = MakeCRgraph(U,E,P,Lc,vCoarse,plongTab);
    L_ = ProlongP(L,Lcr,vCR2Fine);
    gP_ = gP;
    gP_.numRelax = 5;
    Lcr = RelaxGraph(Ucr,Ecr,Pcr,Lcr,gP_);
    L = ProlongP(L,Lcr,vCR2Fine);
    EnergyAssert(U,E,P,L,U,E,P,L_);
end



%
% assertion - remove later
% EnergyAssert(U,E,P,L,Uc,Ec,Pc,Lc);
    

%
% relax the graph
[L, e_, nv_] = RelaxGraph(U,E,P,L,gP);
Q.e = [Q.e, e_];
Q.nv = [Q.nv, nv_];
prior = UpdatePrior(prior,L,gP);


end