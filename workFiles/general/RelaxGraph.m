function [L, e, v] = RelaxGraph(U,E,P,L,gP,varargin)

%
% RelaxGraph(U,E,P,L,gP) -  relax the graphical model given by U,E,P using
% the algorithm specified by gP. Assuming initial labeling given by L (can
% be an empty matrix)

if ~isempty(varargin)
    gP.numRelax = varargin{1};
else
    dTargetEng = -Inf;
end

e = [];
v = [];

if any(L) && ~isempty(U)
    % no point in relaxing if there's no initial labeling
    
    e = Energy(U,E,P,L);            % for assertion - remove later
    
    switch (gP.relaxType)
        
        case 'MQPBO'
            L = SwapMQPBO(U,E,P,L,gP.numRelax,dTargetEng,gP);
        
        case 'QPBO'
            L = SwapQPBO_(U,E,P,L,gP.numRelax,dTargetEng,gP);
            
        case 'LSA'
            if ((gP.bFineGraph) && (numel(L) == prod(gP.imSz)))
                L = reshape(L,gP.imSz);
            end
            L = lsa_tr(U,E,P,L);
            L = L(:);
            
    end
    e_ = Energy(U,E,P,L);           % for assertion - remove later
    
    e = round(e * 1e6) / 1e6;
    e_ = round(e_ * 1e6) / 1e6;
    assert(e_ <= e);
    
    e = e_;
    v = size(U,1);
    
end