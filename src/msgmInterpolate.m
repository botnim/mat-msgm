function x = msgmInterpolate(xc, mapFineToCoarse, mapInterpolation)
% Prolong(Lc,plongTab) - given a labeling 'Lc' of the coarse graph,
% interpolate the assignment according to the 'plongTab'
%

    x = xc(mapFineToCoarse);

    inds = sub2ind(size(mapInterpolation), ...
            (1 : numel(mapFineToCoarse))', ...
            x);
    
    x = mapInterpolation(inds);

% N = length(Lc) + size(plongTab,1);  % total number of variables
% L = zeros(N,1);
% 
% L(vCoarse) = Lc;                    % map interpolators to their vertex
% 
% ii = plongTab(:,2);                 % interpolators
% jj = plongTab(:,1);                 % interpolants
% plongTab = plongTab(:,3:end);       % prolongation table
% 
% idx = sub2ind(size(plongTab),1:size(plongTab,1),Lc(ii)');
% L(jj) = plongTab(idx);

end