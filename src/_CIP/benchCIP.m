function [NULL, flag] = benchCIP(gP, lP, TAG, ~)

%
% running simulations on CIP dataset


rng('shuffle');
NULL = rand(1);                    	% first random number, make sure shuffle is on

%addpath(genpath('../'));                    % main files of msmrf
addpath(genpath('../../datasets/CIP'));   	% chinese character dataset
%addpath(genpath('~/mosek'));
%addpath(genpath('../../software'));

%
% run simulation
% [~, factorGraph] = fge_read(lP.filename);       % load file
% [D, E, W] = fg2std(factorGraph);                % "standard" form
% save(['data/files/',lP.filename(1:end-3),'mat'],...
%     'D','E','W');
ld = load([lP.filename(1:end-3),'mat']);
U = ld.D;
E = ld.E;
P = ld.W;

E_ = E;

% read image size (for LSA-TR relaxation)
sName = lP.filename(1:end-4);
k_ = strfind(sName,'_');
szDim1 = str2double(sName(k_(end)+1:end));
szDim2 = str2double(sName(k_(end-1)+1:k_(end)-1));

%
% transfer data structures to new form
[ii, jj] = find(E);
E = [ii, jj];
P = permute(P,[2,3,1]);


% % Nowozin's energy
% load cip_results_map.mat;                       % load results
% e(1) = fge_energy(1,factorGraph,lP.sol);
% e(2) = Energy_gen(D,E,W,lP.sol(:));             % sanity check for D,E,W
%                                                 % --were they transferred
%                                                 % correctly?


% % Shai and Meirav's energy
% [U,P] = ReparamGraph(U,E,P,gP);
% V = ones(2) - eye(2);
% P = bsxfun(@minus,P,P(1,1,:));   % zero out diagonal elements of P
% vals = squeeze(P(1,2,:));
% W = sparse(E(:,1),E(:,2),vals,size(U,1),size(U,1));
% W = W + W';
% [cP, cD, cW] = buildEnergyPyramid(U, V, W, 200);
% sm_l = discreteMultiscaleOptimization...
%     (cP, cD, V, cW, @(u, v, w, il, varargin) swap_qpbo_gen(u, v, w, il,1));


%e = zeros(1,11);
%t = zeros(1,11);

gP_ = msgmParams();

% UAI output
% sfileout = fullfile(pwd,'uai',strcat(lP.filename(1:end-3),'uai'));
% mat2uai(U,E,P,sfileout);

% tic;
% gP_.relaxType = 'QPBO';
% [L, e(1), ~, ~] = msmrf(U,E,P,[],gP_);
% t(1) = toc;

%[~, L] = min(U,[],2);

%%%% parll - seq
%gP_.relaxType = 'LSA-TR';
%gP_.numRelax = 1;
%gP_.numV = 1;
%gP_.prUpdate = 'none';
%
%for i = 1 : 10
%	tic;
%	[~, e(i), ~, ~] = msmrf(U,E,P,L,gP_);
%	t(i) = toc;
%end

gP_.numV = 2;
gP_.relaxType = 'LSA';
gP_.imSz = [szDim2, szDim1];
gP_.bPlongCR = true;
tic;
[~, e, ~] = msgm(U,E,P,[],gP_);
t = toc;
%%%%%%%% end parll - seq




%%%%% IRAD
% gP_.relaxType = 'QPBO';
% gP_.numRelax = 1;
% gP_.numV = 1;
% gP_.prUpdate = 'none';
% % gP_.bPlongCR = true;
% [~, e(1), ~, ~] = msmrf(U,E,P,L,gP_);
% 
% gP_.numV = 2;
% [~, e(2), ~, ~] = msmrf(U,E,P,L,gP_);
%%% end IRAD

% ld = load(['results_/result_',num2str(TAG),'.mat']);
% for i = 1 : (10 * gP_.numV)
% 
%     ig = qpboms(U,E,P,ig,gP_);
%     ig = ig + 1;
% end

% e = Energy(U,E,P,sm_l{1});
e = round(e * 1e2) / 1e2;


%
% save results
TAG = lP.itr;
save([gP.dirName,'/result_',num2str(TAG),'.mat'],'e','t');

flag = 1;
