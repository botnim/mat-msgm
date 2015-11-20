function x = msgmSwap(G, x, param)
% msgmSwap(G, x, param) ab-Swap for "hard" energies with general pairwise potentials
%
    % possible a/b combinations
    ab_combs = combnk(1:G.numLabels, 2);

    % do param.numSwapIterations iterations
    for i = 1 : param.numSwapIterations

        % random ordering of ab_combs
        perm = randperm(size(ab_combs, 1));
        for j = 1 : size(ab_combs, 1)

            a = ab_combs(perm(j),1);
            b = ab_combs(perm(j),2);
            x = swap(G, x, a, b, param);
        end
    end
end


%% Swap-QPBO

function xnew = swap(G, x, a, b, param)
% optimize by allowing only swap move a <-> b

    % find variables whose label is a or b
    vb_ab = (x == a) | (x == b);
    xnew = x;
    
    if (any(vb_ab))
        
        % fix labels of variables whose labels is not a,b
        [Gcond, xcond] = msgmConditionalDist(G, x, ~vb_ab);
        
        % only (a <--> b) swap move is allowed
        Gcond.u = Gcond.u(:,[a,b]);
        Gcond.p = Gcond.p([a,b],[a,b],:);
        xcond = 1 * (xcond == a) + 2 * (xcond == b);
        
        % improve a,b initial guess
        if (strcmp(param.optimization, 'QPBO'))

            % prepare pairwise potentials for QPBO wrapper
        
            xab = msgmQPBO(Gcond, xcond);
        else
            
            xab = msgmLSA(Gcond, xcond, param);
        end

        % insert optimized a,b move (xab) into the initial labeling x
        xab = a * (xab == 1) + b * (xab == 2);
        xnew(vb_ab) = xab;

        % assert that energy does not increase
        assert(msgmEnergy(G, xnew) <= msgmEnergy(G, x));
    end
end
