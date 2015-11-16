%
% toy problem for msmrf
%
% Run simulations of msmrf on a 4-connected grid of size 
% GRID_SIZE x  GRID_SIZE, with each variable taking label in 1,...,K.
% Repeat this for 'numTests' times.
%
% NOTE: the main purpose of this simulation is to test the coarsening
% scheme by counting the number of times that the multiscale framework was
% able to recover the global minimum (without relaxations!). For this
% purpose an exhaustive, brute-force search is carried out inside the loop.
% This brute-force search becomes a computational overload very quickly,
% when the number of labels increases, or more importantly, when increasing
% the grid size.
%

GRID_SIZE = 3;    	% 4-connected grid GRID_SIZE by GRID_SIZED
K = 2;              	% number of labels
numTests = 500;     	% number of tests

%
% initialization of some vars for adjacency matrix
sz = [GRID_SIZE, GRID_SIZE];
[ii, jj] = sparse_adj_matrix(sz, 1, 1);
sel = ii<jj;
ii = ii(sel);
jj = jj(sel);

%
% make data matrices
E = [ii, jj];
M = size(E,1);      % number of edges
N = GRID_SIZE^2;    % number of variables


% %
% % set global parameters
% gP = setParams;

% eMin = zeros(numTests,1);
% load(['eMin_N',num2str(N),'_K',num2str(K),'.mat']);
% eMin = round(eMin * 10) / 10;

vEnergy = zeros(numTests,1);

count  = 0;
eTotal = 0;
eMS = zeros(numTests,1);
eBL = zeros(numTests,1);
for i = 1 : numTests

    rng(i);
    
    U = randn(N,K);
    U = round(U, 1);

    P = randn(K,K,M);
    P = round(P, 1);
    
% %     ig = int32(randi(2,size(U,1),1));
% %     
% %     tic;
% %     u = U';
% %     p = [E(:,1)'; E(:,2)'; ...
% %             squeeze(P(1,1,:))'; ...         	% Eij(0,0)
% %             squeeze(P(1,2,:))'; ...          	% Eij(0,1)
% %             squeeze(P(2,1,:))'; ...       	% Eij(1,0)
% %             squeeze(P(2,2,:))'];    
% % 
% %     u(isinf(u)) = 1e3;
% %     p(isinf(p)) = 1e3;
% %     
% %     ig_ = ig;
% %     for j = 1 : 6
% %         
% %         x = QPBO_wrapper_mex(u, p, ig_, 'i');
% %         ig_ = x+1;
% %     end
% %     t1 = toc;
% % 
% %     tic;
% %     y = qpboms(U,E,P,ig,gP);
% %     t2 = toc;
% % 
% %     e1 = Energy(U,E,P,double(x+1));
% %     e2 = Energy(U,E,P,double(y+1));
% %     
% %     eng(i,:) = [e1, e2];
% %     time(i,:) = [t1, t2];
    
    
    % find energy for k-label representation
    gP = msgmParams;
    gP.imSz = [GRID_SIZE, GRID_SIZE];
    gP.relaxType = 'NONE';
    gP.numV = 1;
    gP.numRelax = 1;
    gP.bPlongCR = false;
    tic;

    [Lk, ~, ~] = msgm(U,E,P,[],gP);
    eMS(i) = round(Energy(U,E,P,Lk) * 1e2) / 1e2;
    tMS = toc;
    
    tic;
    gP = setParams;
    gP.imSz = [GRID_SIZE, GRID_SIZE];
    gP.relaxType = 'NONE';
    gP.numV = 1;
    gP.numRelax = 1;
    gP.bPlongCR = false;
    [Lk, ~, ~] = msmrf(U,E,P,[],gP);
    eBL(i) = round(Energy(U,E,P,Lk) * 1e2) / 1e2;
    tBL = toc;
    
%     gP = setParams;
%     gP.imSz = [GRID_SIZE, GRID_SIZE];
%     gP.relaxType = 'LSA';
%     gP.numV = 0;
%     gP.bPlongCR = false;
%     tic;
%     [Lk, ~, ~] = msgm_REAL(U,E,P,[],gP);
%     eSS = round(Energy(U,E,P,Lk) * 1e2) / 1e2;
%     tSS = toc;

%     disp('------');
%     disp(strcat('multiscale energy:  ',num2str(eMS)));
%     disp(strcat('baseline energy:    ',num2str(eBL)));
%     disp(strcat('singlescale energy: ',num2str(eSS)));
%     disp('------');
%     disp(strcat('multiscale time:    ',num2str(tMS)));
%     disp(strcat('baseline time:      ',num2str(tBL)));
%     disp(strcat('singlescale time:   ',num2str(tSS)));
    
    % ... same but with MQPBO
%     gP.MQPBO = true;
%     gP.numRelax = 10;
%     tic;
%     [L_, ~, ~] = msmrf(U,E,P,[],gP);
%     e_ = round(Energy(U,E,P,L_) * 1e2) / 1e2;
%     t_ = toc;
%     
%     tmp(i,:) = [ek, e_];
%     time(i,:) = [tk,t_];

% %     % find energy for 2-label representation
% %     tic;
% %     gP.b2lbl = true;
% %     gP.CONST = 100;
% %     gP.bEdgeThrs = 1;
% %     gP.binEdges = 5;
% %     gP.MQPBO = true;
% %     [U_, E_, P_, const] = kto2(U,E,P,gP);
% %     tic;
% %     [L2, ~, ~] = msmrf(U_,E_,P_,[],gP);
% %     e2 = round((Energy(U_,E_,P_,L2) + const) * 1e2) / 1e2;
% %     t2 = toc;
% %     L2_ = labels2tok(L2,K);  
% %     e2_ = round(Energy(U,E,P,L2_) * 1e2) / 1e2; 
% %     assert(e2 == e2_);
% %     
% %     
% %     % find energy for 2-label representation
% %     tic;
% %     gP.b2lbl = true;
% %     gP.CONST = 100;
% %     gP.bEdgeThrs = 1;
% %     gP.binEdges = 5;
% %     [U_, E_, P_, const] = kto2exp(U,E,P,gP);
% %     tic;
% %     [L2, ~, ~] = msmrf(U_,E_,P_,[],gP);
% %     e2exp = round((Energy(U_,E_,P_,L2) + const) * 1e2) / 1e2;
% %     t2exp = toc;
% %     L2_ = labels2tokexp(L2,K);  
% %     e2_ = round(Energy(U,E,P,L2_) * 1e2) / 1e2; 
% %     assert(e2exp == e2_);
    
% %     tmp(i,:) = [ek, e2, e2exp];
    
    a = 1;
    
    
%     % find minimum energy
%     eMin_ = realmax;
%     c = CounterInBase(K);
%     for j = 1 : K ^ size(U,1)
%         
%         L = zeros(1,size(U,1));
%         tmpL = c.next();
%         L(end-length(tmpL)+1:end) = tmpL;	% next chronological labeling
% 
%         L = L' + 1;
%         
%         tmpEnergy = round(Energy(U,E,P,L) * 1e2) / 1e2;
%         if (tmpEnergy < eMin_)
%             eMin_ = tmpEnergy;
%             vLmin = L;
%         end
%     end
%     
%     eMin(i) = eMin_;  
% 
%     if ~isequal(e,eMin(i))
%         assert(e > eMin(i));
%         count = count + 1;
%     end
    
end

