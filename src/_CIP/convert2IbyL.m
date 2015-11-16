function mrf = convert2IbyL(U,E,P,gP,filename)
    
%     % comment out
%     filename = 'TST_test_0001_88_115.mat';
%     gP = setParams;
% 
%     % transform data structures to new form
%     % ld = load([lP.filename(1:end-3),'mat']);
%     ld = load(filename);
%     U = ld.D;
%     E = ld.E;
%     P = ld.W;
% 
%     E_ = E;
% 
%     [ii, jj] = find(E);
%     E = [ii, jj];
%     P = permute(P,[2,3,1]);

    % reparameterization
    [U, P] = ReparamGraph(U,E,P,gP);
    P = bsxfun(@minus,P,P(1,1,:));
    weight = P(1,2,:);
    weight = weight(:);

    % extract dimensions
    inds = strfind(filename,'_');
    indDot = strfind(filename,'.');
    dim2 = str2double(filename(inds(end-1)+1:inds(end)-1));
    dim1 = str2double(filename(inds(end)+1:indDot-1));

    % build the mrf struct
    [~, initialGuess] = min(U,[],2);
    img = repmat(reshape(initialGuess - 1,dim2,dim1)',[1,1,3]);
    mrf.numCol = int32(dim2);
    mrf.numLine = int32(dim1);
    mrf.numLabels = int32(2);
    mrf.numEdges = int32(size(P,3));
    mrf.order = int32(1);
    mrf.unary = int64(U * 1e10);
    mrf.active = true(size(U));
    mrf.edges = int32(E)';
    mrf.weight = int64(weight * 1e10)';
    mrf.dist = int64(ones(2) - eye(2));
    mrf.labels = int32(initialGuess - 1);
    mrf.img = double(img);
end
