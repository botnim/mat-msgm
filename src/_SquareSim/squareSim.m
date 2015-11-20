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
K = 4;              	% number of labels
numTests = 10;           % number of tests

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


vEnergy = zeros(numTests,1);

count  = 0;
eTotal = 0;
eMS = zeros(numTests,1);
eBL = zeros(numTests,1);
for i = 1 : numTests

    disp(num2str(i));
    rng(i);
    
    G.u = round(randn(N,K), 1);
    G.p = round(randn(K,K,M), 1);
    G.adj = E;
    
    % msgm
    param = msgmParams;
    param.imSz = [GRID_SIZE, GRID_SIZE];
    param.optimization = 'LSA';
    param.numSwapIterations = 1;
    param.bSoftInterpolation = false;

    
    [~, eMS(i), tMS] = msgm(G, [], param);
    param.optimization = 'QPBO';
    [~, eBL(i), tBL] = msgm(G, [], param);
    
end

