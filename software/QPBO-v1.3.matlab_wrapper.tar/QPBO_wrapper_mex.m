%
% A wrapper for QPBO library by Vladimir Kolmogorov.
%
% Usage: 
%   [x] = QPBO_wrapper_mex(UTerm, PTerm, ig, method );
%
% Inputs:
%  UTerm - 2xN matrix of unary terms (N-number of variables)
%          Each term (col) is a pair [Ei(0), Ei(1)].'
%  PTerm - 6xE matrix of pairwise terms (E-number of pair-wise relations (edges))
%          Each term (col) is of the form [i j Eij(0,0), Eij(0,1), Eij(1,0), Eij(1,1)].'
%          Where i,j are the indices of the variables (first variable is of index 1, last one is N)
%  ig    - int32 0-based labels for initial guess (for "improve" method)
%  method - What type of method to use for optimization (QPBO/Probe/QPBOI)
%           a character either 'q', 'p' or 'i'
%
%
% Compile using  
% >> compile_qpbo_wrapper
%
%
