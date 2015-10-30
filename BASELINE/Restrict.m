function uc = Restrict(u,vCoarse)

%
% Restrict(u,vCoarse) - restrict the vector 'u', which represents either
% a labeling assignment or a prior of the vertices, to those entries which
% are represented on the coarse graph (given by vCoarse)
%


if isempty(u)
    uc = [];
else
    uc = u(vCoarse,:);
end