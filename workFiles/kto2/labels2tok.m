function L_ = labels2tok(L,k)
%labelsk2to translates a boolean labeling given by reduction of the
%function 'kto2' to the "ordinary" k-labeling

L = reshape(L,k-1,[])' - 1;
L_ = k - sum(L,2);

end

