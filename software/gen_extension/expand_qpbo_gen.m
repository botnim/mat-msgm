function [l e] = expand_qpbo_gen(varargin)
%
% a-Expand on top of QPBO(I) for non-submodular energies
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

l = l(:);

e = Energy_gen(D,E,W,l);

for itr=1:N
%     for a=1:nl    
%                 
%         l = expand(D,E,W,l,a);
%                
%     end
    perm = randperm(nl);
    for i = 1 : length(perm)
        a = perm(i);
        l = expand(D,E,W,l,a);
    end

    ne = Energy_gen(D,E,W,l);
    e = round(e*1e6) / 1e6;
    ne = round(ne*1e6) / 1e6;
    assert(ne <= e); % energy does not increase... .  UNCOMMENT!
    if e-ne < 1e-6
        e = ne;
        break; % too small an improvement
    end
    e = ne;
end


e(1,2)=itr;


function l = expand(D, E, W, l, a)
%
% expand label a in curent solution l
%

sel = (l~=a); % participating nodes

if ~any(sel)
    return;
end
n = size(D,1);

% reduce the energy:



% For variables with nieghbors eq to a unary term must be adapted to
% compensate for possible label change w.r.t non-participating nieghbor

[ii, jj, ww] = find(E+E');

dc = D( sub2ind(size(D), (1:n)', l(:)) );
dc = [dc(:) D(:,a)];


pn = sel(ii) & (~sel(jj)); % problematic neighbors
% right direction
iii = ii(ii<jj);
www = ww(ii<jj);
pn_ = pn(ii<jj);
if any(pn_)
    indtmp = sub2ind(size(W),www(pn_),l(iii(pn_)),a*ones(size(iii(pn_))));
    dc(:,1) = dc(:,1) + accumarray(iii(pn_), W(indtmp), [n 1]);
    indtmp = sub2ind(size(W),www(pn_),a*ones(size(iii(pn_))),a*ones(size(iii(pn_))));
    dc(:,2) = dc(:,2) + accumarray(iii(pn_), W(indtmp), [n 1]);
end
% reverse direction
iii = ii(ii>jj);
www = ww(ii>jj);
pn_ = pn(ii>jj);
if any(pn_)
    indtmp = sub2ind(size(W),www(pn_),a*ones(size(iii(pn_))),l(iii(pn_)));
    dc(:,1) = dc(:,1) + accumarray(iii(pn_), W(indtmp), [n 1]);
    indtmp = sub2ind(size(W),www(pn_),a*ones(size(iii(pn_))),a*ones(size(iii(pn_))));
    dc(:,2) = dc(:,2) + accumarray(iii(pn_), W(indtmp), [n 1]);
end

% indtmp = sub2ind(size(W),ww(pn),a*ones(size(ii(pn))),a*ones(size(ii(pn))));
% dc(:,2) = dc(:,2) + accumarray(ii(pn), W(indtmp), [n 1]);



u = dc(sel,:)'; % unary
[ii, jj, ww] = find(E(sel,sel));
assert(isequal(E(sel,sel),triu(E(sel,sel))));
W_ = W(ww,:,:);
SZ = size(W_);
rl = l(sel);

if isempty(ii)
    p = zeros(6,0);
else
    p = [ii'; jj'; ...
        W_( sub2ind(SZ,(1:SZ(1))',rl(ii),rl(jj)) )'; ...     	% Eij(0,0)
        W_( sub2ind(SZ,(1:SZ(1))',rl(ii),ones(size(jj))*a) )'; 	% Eij(0,1)
        W_( sub2ind(SZ,(1:SZ(1))',ones(size(ii))*a,rl(jj)) )'; 	% Eij(1,0)
        W_(:,a,a)'];                                            % Eij(1,1)
end
ig = zeros(nnz(sel),1,'int32');
% ENG = Energy_gen(D(sel,:),E(sel,sel),W,rl);
x = QPBO_wrapper_mex(u, p, ig, 'i');
rl(x>0) = a;
% NEW_ENG = Energy_gen(D(sel,:),E(sel,sel),W,rl);
% assert(NEW_ENG <= ENG);
l(sel) = rl;
