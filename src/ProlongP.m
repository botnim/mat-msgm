function L = ProlongP(L,Lcr,vCR2Fine)

% ProlongP(Lc,vCoarse,Lcr,vCR2Fine) - prolong a coarse labeling given by
% 'Lc' and integrate the compatiable relaxations labeling of the fine
% scale, given by 'Lcr'
%


L(vCR2Fine) = Lcr;
