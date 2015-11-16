addpath(genpath('C:\DATA\Dropbox (Personal)\Weizmann\project\mat-msgm\datasets\CIP'));    	% dataset folder
addpath(genpath('C:\DATA\Dropbox (Personal)\Weizmann\matlab\WeizGrid'));         	% WeizGrid on UNIX
%addpath(genpath('~/matlab/msmrf/workFiles'));
%addpath(genpath('../../externals'));

%
% set simulation parameters
NAME = 'CIP';
load cip_results_map.mat;
lP(length(instance_files)).filename = [];
[lP.filename] = deal(instance_files{:});
[lP.sol] = deal(sol{:});
[lP.k] = deal(1);

lP = lP(9);

lP = lP(1:end);
itr = num2cell(1:numel(lP));
[lP.itr] = deal(itr{:});
% lP = lP([16,29,46,64,79,99]);
% lP = lP(1);

%
% set global parameters
gP = setParams;


%
% create folder for results
gP.path = [pwd,'/'];
gP.dirName = 'results';
if not(exist([gP.path,gP.dirName],'dir'))
    mkdir([gP.path,gP.dirName]);
end
save([gP.path,gP.dirName,'/INFO.mat'],'gP');


%
% run the code on the cluster
WGjob = WGexec('nparallels', 115, 'Name', NAME, 'Mail', 'oymeir@gmail.com', ...
    'WorkFunc', 'benchCIP', 'LocalDebug', ispc, 'GlobalParams', gP, ...
    'SubParams', lP, 'RngShuffle', true, 'WaitTillFinished', false);

%
% finish up
% wrapSim(gP);



