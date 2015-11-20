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

GRID_SIZE = 100;    	% 4-connected grid GRID_SIZE by GRID_SIZED
K = 2;              	% number of labels
numTests = 50;     	% number of tests

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
for i = 1 : numTests

    rng(i);
    
    gP = setParams;

    U = randn(N,K);
    U = round(U * 1e1) / 1e1;

    P = 2 * randn(K,K,M);
    P = round(P * 1e1) / 1e1;
        
    % find energy for k-label representation
    gP = setParams;
    gP.imSz = [GRID_SIZE, GRID_SIZE];
    gP.relaxType = 'LSA';
    gP.numV = 1;
    tic;
    [~, eMS(i), ~] = msmrf(U,E,P,[],gP);
    
    gP.relaxType = 'QPBO';
    [~, eQPBO(i), ~] = msmrf(U,E,P,[],gP);
%     
%     disp(strcat('multiscale energy:  ',num2str(eMS)));
%     disp(strcat('multiscale time:    ',num2str(tMS)));
%     disp(strcat('singlescale energy: ',num2str(eSS)));
%     disp(strcat('singlescale time:   ',num2str(tSS)));
    
end

