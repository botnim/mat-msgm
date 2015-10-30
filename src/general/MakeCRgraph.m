function [Ucr, Ecr, Pcr, Lcr, vCR2Fine, L] = ...
    MakeCRgraph(U,E,P,Lc,vCoarse, plongTab)

%
% MakeCRgraph(U,E,P,Lc,vCoarse) - build a graphical model for compatible
% relaxations; it will be used to interpolate a coarse labeling to a fine
% scale in the probabilistic prolongation scheme
%

M = size(E,1);
N = size(U,1);

L = zeros(N,1);
L(vCoarse) = Lc;            % interpolated coarse labels

vbFine = true(N,1);
vbFine(vCoarse) = false;    % indicator of fine nodes (that go to CR graph)

% % % interpolate fine vertices whose interpolator's label hadn't changed
% % % get indices of fine vertices that should be interpolated
% % vbPlongVertex = (plongTab(vbFine,2) == L(plongTab(vbFine,1)));
% % vTmpInd = find(vbFine);
% % % interpolate these vertices according to rule
% % L(vTmpInd(vbPlongVertex)) = plongTab(vTmpInd(vbPlongVertex),3);
% % vbFine(vTmpInd(vbPlongVertex)) = false; % these vertices don't go in CR graph

if ~any(vbFine)
    % all vertices have been labeled, terminate
    
    Ucr = [];
    Ecr = [];
    Pcr = [];
    vCR2Fine = [];
    Lcr = [];
    
else
    % prepare compatible relaxations graph
    
    % set maps Fine <-> CR, and CR data structures
    vCR2Fine = find(vbFine);
    vFine2CR = zeros(1,N);                      % vFine2CR holds a map of vertex inds
    vFine2CR(vCR2Fine) = 1 : length(vCR2Fine);	% from the fine graph to the CR graph

    Ucr = U(vCR2Fine,:);
    vEdgesInd = zeros(M,1);     % vector of edge indices, keeps track of which
                                % edges are present in the CR graph

    cntrE = 1;                  % counter for #edges in CR graph



    % loop over edges,
    % to consider all relevant factors on the CR graph
    for iE = 1 : M

        ii = E(iE,1);
        jj = E(iE,2);

        if (vbFine(ii) && vbFine(jj))
            % both vertices belong to the fine graph

            vEdgesInd(cntrE) = iE;
            cntrE = cntrE + 1;

        elseif vbFine(ii)
            % only the first vertex is in the fine graph

            pairwise = P(:,:,iE);
            ii_ = vFine2CR(ii);     % ii's index w.r.t CR graph
            Ucr(ii_,:) = Ucr(ii_,:) + pairwise(:,L(jj))';

        elseif vbFine(jj)
            % only the second vertex is in the fine graph

            pairwise = P(:,:,iE);
            jj_ = vFine2CR(jj);     % jj's index w.r.t CR graph
            Ucr(jj_,:) = Ucr(jj_,:) + pairwise(L(ii),:);

        end

    end

    vEdgesInd = vEdgesInd(1:cntrE-1);
    Ecr = E(vEdgesInd,:);
    Ecr = vFine2CR(Ecr);
    Pcr = P(:,:,vEdgesInd);
    [~, Lcr] = min(Ucr,[],2);   % !!!! labels must be initialized for SwapQPBO !!!!!
    
end




