function [NULL, flag] = benchCIP(gP, lP, TAG, ~)

%
% running simulations on CIP dataset


rng('shuffle');
NULL = rand(1);                    	% first random number, make sure shuffle is on

addpath(genpath('../'));                    % main files of msmrf
addpath(genpath('../../datasets/CIP'));   	% chinese character dataset
addpath(genpath('~/mosek'));
addpath(genpath('../../software'));

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


% e = zeros(1,1);
% t = zeros(1,1);

gP_ = setParams;

%% IbyL

% Pruning classifiers
lambda          = 0.1;     
data            = load(['externals/IbyL/Learning/Classifiers/Stereo/FastPD/Pruning_1.mat']);
classifiers     = data.classifier;
numScale        = data.options.numScale;
solver          = int32(0);

% Parameters
options               	= struct;
options.numScale      	= numScale;
options.solver       	= solver;
options.iterMax     	= int32(10);
options.pruning         = false;
options.multigrid       = false;
options.exportFeatures  = false;

mrf = convert2IbyL(U,E,P,gP_,strcat(lP.filename(1:end-3),'mat'));

% Muliscale optimization
clc;
options.exportFeatures  = false;
options.numScale      	= numScale;
t = tic;
[r_MS,f_MS,t_MS,m_MS] = IbyL_mex(mrf,options);
t(1) = toc(t);
e(1) = r_MS.nrg / 1e10;

% IbyL optimization (pruning)
options.pruning         = true;
options.classifiers     = classifiers;
t = tic;
[r_PR,f_PR,t_PR,m_PR] = IbyL_mex(mrf,options);
t(2) = toc(t);
e(2) = r_PR.nrg / 1e10;
%%%%%%%%%% end IbyL

% tic;
% gP_.relaxType = 'QPBO';
% [L, e(1), ~, ~] = msmrf(U,E,P,[],gP_);
% t(1) = toc;

% for i = 1 : 10
% 	gP_.relaxType = 'LSA';
% 	gP_.imSz = [szDim2, szDim1];
% 	gP_.numV = 1;
% 	tic;
% 	[~, e(i), ~, ~] = msmrf(U,E,P,[],gP_);
% 	t(i) = toc;
% end

%gP_.relaxType = 'LSA';
%gP_.imSz = [szDim2, szDim1];
%gP_.numV = 10;
%gP_.prUpdate = 'none';
%gP_.bPlongCR = true;
%[L, e(11), ~, ~] = msmrf(U,E,P,[],gP_);


% ld = load(['results_/result_',num2str(TAG),'.mat']);
% for i = 1 : (10 * gP_.numV)
% 
%     ig = qpboms(U,E,P,ig,gP_);
%     ig = ig + 1;
% end

% e = Energy(U,E,P,sm_l{1});
e = round(e * 1e2) / 1e2;

e
t
disp('done.');


%
% save results
save([gP.dirName,'/result_',num2str(TAG),'.mat'],'e');

flag = 1;
