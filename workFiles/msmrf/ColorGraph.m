function L = ColorGraph(Lc,vCoarse,plongTab)

%
% Prolong(Lc,plongTab) - given a labeling 'Lc' of the coarse graph,
% interpolate the assignment according to the 'plongTab'
%

Nc = length(Lc);
N = Nc + size(plongTab,1);          % total number of variables
L = zeros(N,1);
Lc = (1:Nc)';

L(vCoarse) = Lc;                    % map interpolators to their vertex

ii = plongTab(:,2);                 % interpolators
jj = plongTab(:,1);                 % interpolants

L(jj) = L(ii);

end