function lP = setLocalParams(sP)

%
% "flatten" the simulation parameters, such that they can be written in a
% single loop
%
%
% example:  if we have parameters (alpha,beta) and wish to run simulations
% for every combination of the parameters, this would usually be written
%
% 
%                 for iAlpha = 1 : length(alpha)
%                     for iBeta = 1 : length(beta)
%                         doSim(alpha(iAlpha),beta(iBeta));
%                     end
%                 end
%
%
% after "flattening" the loop, this becomes
%
%                 for iLoop = 1 : prod([length(alpha),length(beta)])
%                     doSim(lP(iLoop).alpha, lP(iLoop).beta);
%                 end
%
% input:
%   sP  -   sturcture, fields correspond to parameters
%           i.e.    sP.alpha = [1,2,3];  sP.beta = [{'a'}, {'b'}];
%
% output:
%   lP  -   structure with same fields as sP, keeping all possible combs
%           i.e., for input as above, the output is:
%
%                     lP(1)           alpha: 1
%                                     beta: {'a'}
% 
%                     lP(2)           alpha: 2
%                                     beta: {'a'}
% 
%                     lP(3)           alpha: 3
%                                     beta: {'a'}
% 
%                     lP(4)           alpha: 1
%                                     beta: {'b'}
% 
%                     lP(5)           alpha: 2
%                                     beta: {'b'}
% 
%                     lP(6)           alpha: 3
%                                     beta: {'b'}



%
% flatten parameters
cPars = struct2cell(sP);                % param values in cell array
ind = cell(1,length(cPars));
for i = 1 : length(cPars)
    ind{i} = 1 : length(cPars{i});
end
combs = cell(1,length(cPars));          % cell array for all parameters
[combs{:}] = ndgrid(ind{:});          	% create all combinations

%
% create parameter structure with proper fields
lP = struct;
parNames = fieldnames(sP);
lP(length(combs{1}(:))).(parNames{1}) = [];     % initialization of structure
for i = 1 : length(parNames)
    tmp = num2cell(cPars{i}(combs{i}(:)));
	[lP.(parNames{i})] = deal(tmp{:});
end
    
    
    
    




