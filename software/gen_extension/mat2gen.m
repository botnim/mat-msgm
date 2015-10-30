function [E, W_] = mat2gen(W,V)

%
% mat2gen(mW,mV)
%
% switch from a special form of pairwise term to a general form;
% given the connectivity weight matrix mW, and the label interaction
% matrix mV, output connectivity matrix mE and interaction matrix mW_ 
% in general form
%

%
%
W = W - diag(diag(W));

%
% make connectivity matrix
[ii, jj, wij] = find(triu(W));
E = sparse(ii,jj,1:length(ii),size(W,1),size(W,2));
% E = E + E';


%
% make pairwise matrix
W_ = repmat(reshape(V,[1,size(V)]),[length(wij),1,1]);
W_ = bsxfun(@times,wij,W_);



