function gP = setParams()

%
% setting some default parameters for Vcycle/msmrf

gP.numV = 5;                % (int)     number of V-cycles
gP.crsRatioThrs = 1;% WAS 0.85;  	% (0<k<1)   set coarsening ratio
gP.bigAgg = 1;              % (1/0)     use big aggregators
gP.numRelax = 1;            % (int)     number of relaxations
gP.binEdges = 20;           % (int)     bin edge-scores in ScoreEdges() for
                            %           randomness;
                            %           if 0  - do not bin
                            %           k > 0 - use 'k' bins for binning
gP.prUpdate = 'none';       %
gP.wDegEnt = 1;


gP.bPlongCR = true;        % prolonb via Compatible Relaxations
gP.bEdgeThrs = 0;           % (0 <= val <= 1) keep edges with <=  val of entropy
gP.bAggMode = 0;            % (0/1) 0 for regular edge contraction
                            %       1 for big aggregator identification
                                                        
