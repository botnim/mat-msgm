function L = SwapQPBO(U,E,P,L,N,varargin)

% ab-Swap on top of QPBO(I) for non-submodular energies
%

if ~isempty(varargin)
    dTargetEng = varargin{1};
else
    dTargetEng = -Inf;
end


K = size(U,2);

e = Energy(U,E,P,L);

combs = combnk(1:K,2);
combs = combs(randperm(size(combs,1)),:);
for itr = 1 : N

    perm = randperm(size(combs,1));
    for i = 1 : length(perm)
        a = combs(perm(i),1);
        b = combs(perm(i),2);
        L = swap(U,E,P,L,a,b);
    end
    
    ne = Energy(U,E,P,L);
    e = round(e*1e6) / 1e6;
    ne = round(ne*1e6) / 1e6;
    assert(ne <= e);        % energy does not increase...
    if (e == ne) || (ne <= dTargetEng)
%         e = ne;
        break;              % too small an improvement
    end
    e = ne;
end


function L = swap(U,E,P,L,a,b)
%
% swap labels a and b in curent solution l
%

sel = (L==a) | (L==b); % participating nodes
if ~any(sel)
    return;
end
n = size(U,1);

% reduce the energy:
% For variables with nieghbors not eq to a/b unary term must be adapted to
% compensate for possible label change w.r.t non-participating nieghbor

M = size(E,1);
ii = E(:,1);
jj = E(:,2);
ww = (1 : M)';

% right direction
pn = sel(ii) & (~sel(jj));          % problematic neighbors
if any(pn)
    indtmp = sub2ind(size(P),a*ones(nnz(pn),1),L(jj(pn)),ww(pn));
    U(:,a) = U(:,a) + accumarray(ii(pn), P(indtmp), [n 1]);
    indtmp = sub2ind(size(P),b*ones(nnz(pn),1),L(jj(pn)),ww(pn));
    U(:,b) = U(:,b) + accumarray(ii(pn), P(indtmp), [n 1]);
end
% reverse direction
pn = sel(jj) & (~sel(ii));
if any(pn)
    indtmp = sub2ind(size(P),L(ii(pn)),a*ones(nnz(pn),1),ww(pn));
    U(:,a) = U(:,a) + accumarray(jj(pn), P(indtmp), [n 1]);
    indtmp = sub2ind(size(P),L(ii(pn)),b*ones(nnz(pn),1),ww(pn));
    U(:,b) = U(:,b) + accumarray(jj(pn), P(indtmp), [n 1]);
end

u = U(sel,[a b])'; % unary
% [ii, jj, ww] = find(E(sel,sel));

indsel = zeros(n,1);
indsel(sel) = 1 : nnz(sel);     % map from all vertices to their new index
                                % as 'selected' vertices
edges = sel(ii) & sel(jj);      % indices of edges in the proxy graph
ii_ = indsel(ii(edges));      	% left vertex of selected pairs
jj_ = indsel(jj(edges));      	% right vertex of selected pairs
ww_ = ww(edges);                % indices of edges of selected pairs


if isempty(ww_)
    p = zeros(6,0);
else
    p = [ii_'; jj_'; ...
        squeeze(P(a,a,ww_))'; ...         	% Eij(0,0)
        squeeze(P(a,b,ww_))'; ...          	% Eij(0,1)
        squeeze(P(b,a,ww_))'; ...       	% Eij(1,0)
        squeeze(P(b,b,ww_))'];              % Eij(1,1)
end
ig = int32(L(sel)==b);
x = QPBO_wrapper_mex(u, p, ig, 'i');
rl = a*ones(numel(x),1);
rl(x > 0) = b;
L(sel) = rl;
