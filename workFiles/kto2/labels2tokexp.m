function L_ = labels2tokexp(L,k)
%labelsk2to translates a boolean labeling given by reduction of the
%function 'kto2' to the "ordinary" k-labeling

L = reshape(L,k-1,[])' - 1;
L_ = zeros(size(L,1),1);

idx1 = ~any(L,2);   % all 0's map to 1st label
L_(idx1) = 1;

[~, vLbl] = max(L(~idx1,:),[],2);
L_(~idx1) = vLbl + 1;

end

