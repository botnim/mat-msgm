% playing around with the following structure:
%
% o-o-o       2-1-3
%   |           |
%   o           4
%
% in particular, want to see if exact MAP recovery is possible, no matter
% what the unary and pairwise are set to be

K = 2;              % number of labels
numTests = 1000;    % number of tests

addpath(genpath('../'));
gP = setParams;

%
% build the graph
E = [1, 2; ...
    1, 3; ...
    1, 4];
% P = repmat(ones(2)-eye(2),[1,1,3]);
% U = [4, 0; ...
%     0, 1; ...
%     0, 1; ...
%     0, 1];
P = randi(5,K,K);
U = randi(5,4,K);


count  = 0;
eTotal = 0;
for i = 135 : numTests
    rng(i);

%     U = randn(N,K);
%     P = randn(K,K,M);
%     
%     U = round(U * 1e1) / 1e1;
%     P = round(P * 1e1) / 1e1;

    P = randi(5,K,K,3);
    U = randi(5,4,K);

    
    % find minimum energy
    eMin = realmax;
    c = CounterInBase(K);
    for j = 1 : K ^ size(U,1)
        
        L = zeros(1,size(U,1));
        tmpL = c.next();
        L(end-length(tmpL)+1:end) = tmpL;	% next chronological labeling

        L = L' + 1;
        
        tmpEnergy = round(Energy(U,E,P,L) * 1e2) / 1e2;
        if (tmpEnergy < eMin)
            eMin = tmpEnergy;
            vLmin = L;
        end
    end
    
    % find energy using msmrf
    L_ = Vcycle(U,E,P,[],gP);
    e = round(Energy(U,E,P,L_) * 1e2) / 1e2;
    eTotal = eTotal + e;
    

    if ~isequal(e,eMin)
        assert(e > eMin);
        count = count + 1;
    end
    
end

eTotal
count

