function eng = Energy(U,E,P,L)

%
% computes energy
%

L = L(:)';              % row vector, easier for indexing
eng = 0;

%
% pairwise term
if ~isempty(P)
    idx = sub2ind(size(P), ...
        L(E(:,1)), ...      % label of left vertex
        L(E(:,2)),...       % label of right vertex
        1 : size(E,1));     % index of edge
    eng = sum(P(idx));
end


%
% unary term
idx = sub2ind(size(U),...
    1 : size(U,1), ...  % index of vertex
    L);                 % label of vertex
eng = eng + sum(U(idx));
