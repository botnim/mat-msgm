function L = lsa_tr(U,E,P,L)

% L = lsa_tr(U,E,P)
% run LSA-TR algorithm on energy (2 labels only, at this point)
%

% reparameterization
PE = [E(:,1), E(:,2), ...
        squeeze(P(1,1,:)), ...         	% Eij(0,0)
        squeeze(P(1,2,:)), ...          	% Eij(0,1)
        squeeze(P(2,1,:)), ...             % Eij(1,0)
        squeeze(P(2,2,:))];   
[newUE, newSubPE, newSuperPE, newConst] = reparamEnergy(U', PE);

eng.UE = newUE;
eng.subPE = newSubPE;
eng.superPE = newSuperPE;
eng.constTerm = newConst;

% run LSA-TR
L = LSA_TR(eng, 0, L);

end

