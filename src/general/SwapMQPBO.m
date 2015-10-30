function L = SwapMQPBO(U,E,P,L,N,varargin)

% ab-Swap on top of QPBO(I) for non-submodular energies
%

if ~isempty(varargin)
    dTargetEng = varargin{1};
else
    dTargetEng = -Inf;
end

e = Energy(U,E,P,L);

% reparameterize from k to 2 labels
K = size(U,2);
if (K > 2)
    
    gP.CONST = 1000;
    [U_, E_, P_, const] = kto2(U,E,P,gP);
    L_ = labelskto2(L,size(U,1),K);
end


for itr = 1 : N

    L = swap(U_,E_,P_,L_);
    ne = Energy(U_,E_,P_,L_) + const;
    e = round(e*1e6) / 1e6;
    ne = round(ne*1e6) / 1e6;
    assert(ne <= e);        % energy does not increase...
    if (e == ne) || (ne <= dTargetEng)
%         e = ne;
        break;              % too small an improvement
    end
    e = ne;
end

L = labels2tok(L_,K);

function L = swap(U,E,P,L)
%
% swap labels a and b in curent solution l
%



% working 2 labels

% validate label order
a = 1;
b = 2;

sel = true(size(U,1),1);

u = U';
p = [E(:,1)'; E(:,2)'; ...
        squeeze(P(1,1,:))'; ...         	% Eij(0,0)
        squeeze(P(1,2,:))'; ...          	% Eij(0,1)
        squeeze(P(2,1,:))'; ...       	% Eij(1,0)
        squeeze(P(2,2,:))'];    

u(isinf(u)) = 1e3;
p(isinf(p)) = 1e3;


ig = int32(L(sel)==b);
x = QPBO_wrapper_mex(u, p, ig, 'i');
rl = a*ones(numel(x),1);
rl(x > 0) = b;
L(sel) = rl;
