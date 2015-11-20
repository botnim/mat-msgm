function param = msgmParams()
% msgmParams() setting parameters for msgm()

    param.numVcycles = 1;               % num of V-cycles
    param.numMinVars = 10;              % num of variables for V-cycle stopping criterion
    param.optimization = 'LSA';         % 'QPBO' or 'LSA', 'NONE' for skipping
    param.numSwapIterations = 1;        % num of 'SWAP' iterations (if applicable)
    param.numEntropyBins = 20;          % num of bins for conditional entorpy scores

    param.imSz = [];                    % used for LSA-euc mode, for optimizing grids

    param.bSoftInterpolation = true;    % boolean flag for using soft interpolation
end
% 
% param.bigAgg = 1;              % (1/0)     use big aggregators
%                             %           randomness;
%                             %           if 0  - do not bin
%                             %           k > 0 - use 'k' bins for binning
% param.wDegEnt = 1;
% 
% 
% param.bPlongCR = true;        % prolonb via Compatible Relaxations
% 
                            
                            
                            
% TO KEEP

%gP.bHardInterpolation = ~gP.bPlongCR;
                                                        
