function [l e] = swap_qpbo_gen(varargin)

% ab-Swap on top of QPBO(I) for non-submodular energies
%

if (numel(varargin) == 5)
    % covering case of running Shai's benchmark
    D = varargin{1};
    mV = varargin{2};
    mW = varargin{3};
    il = varargin{4};
    [E, W] = mat2gen(mW,mV);
    N = 100;
elseif (numel(varargin) == 6)
    % covering case of N specified by input
    D = varargin{1};
    E = varargin{2};
    W = varargin{3};
    il = varargin{4};
    N = varargin{6};
else
    % N not specified by input
    D = varargin{1};
    E = varargin{2};
    W = varargin{3};
    il = varargin{4};
    N = 2;
end


[n, nl] = size(D);

if nargin == 2
    l = randi(nl,n,1); % initial guess
else
    l = il;
end

e = Energy_gen(D,E,W,l);

combs = combnk(1:nl,2);

for itr=1:N
%     for a=1:nl
%         for b=(a+1):nl            
%             l = swap(D,E,W,l,a,b);             
%         end
%     end
    
    perm = randperm(size(combs,1));
    for i = 1 : length(perm)
        a = combs(perm(i),1);
        b = combs(perm(i),2);
        l = swap(D,E,W,l,a,b);
    end
    
    ne = Energy_gen(D,E,W,l);
    e = round(e*1e6) / 1e6;
    ne = round(ne*1e6) / 1e6;
    assert(ne <= e); % energy does not increase...  UNCOMMENT!
    if e-ne < 1e-6
        e = ne;
        break; % too small an improvement
    end
    e = ne;
end


e(1,2)=itr;


function l = swap(D,E,W,l,a,b)
%
% swap labels a and b in curent solution l
%

sel = (l==a) | (l==b); % participating nodes
if ~any(sel)
    return;
end
n = size(D,1);

% reduce the energy:



% For variables with nieghbors not eq to a/b unary term must be adapted to
% compensate for possible label change w.r.t non-participating nieghbor



[ii, jj, ww] = find(E+E');

% % pn = sel(ii) & (~sel(jj)); % problematic neighbors
% % Dc(:,a) = Dc(:,a) + accumarray(ii(pn), wij(pn).*va(l(jj(pn))), [n 1]);
% % Dc(:,b) = Dc(:,b) + accumarray(ii(pn), wij(pn).*vb(l(jj(pn))), [n 1]);

pn = sel(ii) & (~sel(jj)); % problematic neighbors
% right direction
iii = ii(ii<jj);
jjj = jj(ii<jj);
www = ww(ii<jj);
pn_ = pn(ii<jj);
if any(pn_)
    indtmp = sub2ind(size(W),www(pn_),a*ones(nnz(pn_),1),l(jjj(pn_)));
    D(:,a) = D(:,a) + accumarray(iii(pn_), W(indtmp), [n 1]);
    indtmp = sub2ind(size(W),www(pn_),b*ones(nnz(pn_),1),l(jjj(pn_)));
    D(:,b) = D(:,b) + accumarray(iii(pn_), W(indtmp), [n 1]);
end
% right direction
iii = ii(ii>jj);
jjj = jj(ii>jj);
www = ww(ii>jj);
pn_ = pn(ii>jj);
if any(pn_)
    indtmp = sub2ind(size(W),www(pn_),l(jjj(pn_)),a*ones(nnz(pn_),1));
    D(:,a) = D(:,a) + accumarray(iii(pn_), W(indtmp), [n 1]);
    indtmp = sub2ind(size(W),www(pn_),l(jjj(pn_)),b*ones(nnz(pn_),1));
    D(:,b) = D(:,b) + accumarray(iii(pn_), W(indtmp), [n 1]);
end

u = D(sel,[a b])'; % unary
[ii, jj, ww] = find(E(sel,sel));
W_ = W(ww,:,:);
if isempty(ww)
    p = zeros(6,0);
else
    p = [ii'; jj'; ...
        W_(:,a,a)'; ...           	% Eij(0,0)
        W_(:,a,b)'; ...            	% Eij(0,1)
        W_(:,b,a)'; ...           	% Eij(1,0)
        W_(:,b,b)'];                % Eij(1,1)
end
ig = int32(l(sel)==b);
x = QPBO_wrapper_mex(u, p, ig, 'i');
rl = a*ones(numel(x),1);
rl(x > 0) = b;
l(sel) = rl;
