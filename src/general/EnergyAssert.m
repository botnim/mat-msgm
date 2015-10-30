function EnergyAssert(U,E,P,L,Uc,Ec,Pc,Lc)

%
% asserts that the prolongated energy (of U,E,P given L) is not greater
% than the energy of the coarse graph (Uc, Ec, Pc, given Lc).
%
% If this rule is violated, it's a strong indication that the coarsening in
% CompCoarseGraph(..) has bugs.

e = Energy(Uc,Ec,Pc,Lc);
e_ = Energy(U,E,P,L);
e = round(e * 1e4) / 1e4;
e_ = round(e_ * 1e4) / 1e4;
if (e_ > e)
    error(['coarse: ',num2str(e),', fine: ',num2str(e_)]);
end

end