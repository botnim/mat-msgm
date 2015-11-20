function wrapSimCIP()

sPath = pwd;
dirName = 'C:\Users\OM\Desktop\results';
numTests = 100;

load('DTF-names.mat');

N = 1;      % number of test elements

% load file
e = zeros(numTests,N);
t = zeros(numTests,N);
for i = 1 : numTests
    
%     if (i == 58) || (i == 5)
%         continue;
%     end
    sName = lP(i).filename(1:end-4);
    k_ = strfind(sName,'_');
    szDim1 = str2double(sName(k_(end)+1:end));
    szDim2 = str2double(sName(k_(end-1)+1:k_(end)-1));

    sFile = [dirName,'/','result_',num2str(i),'.mat'];
    if exist(sFile,'file')

        ld = load(sFile);
        e(i,:) = ld.e;
%         t(i,:) = ld.t;
%         L = ld.L;
%         
%         L = reshape(L,szDim2,szDim1) - 1;
%         L = permute(L,[2,1]);
%         imwrite(L,['imgs/DTF-',num2str(i),'.png']);
        
    end
end

E = e;
E = round(E*1e2) / 1e2;

load([sPath,'CIP_optimal.mat']);
BEST(e(:,1) == 0) = [];
E(E(:,1) == 0) = [];
for i = 1 : N
    disp([num2str(i),' :    ',num2str(nnz(E(:,i) <=  BEST))]);
end

% E = bsxfun(@minus,E,BEST);
E = bsxfun(@minus,E,E(:,2));
[~, ind] = sort(E(:,1),'ascend');
% ind = 1:100;

%
% plot energies
figure; hold on;
% plot(1:100,zeros(size(BEST)),'-ok','linewidth',2); 	% BEST
plot(1:100,E(ind,2),'-or','linewidth',2);         	% parallel
plot(1:100,E(ind,1),'-ob','linewidth',2);           % seq., short
% plot(1:100,E(ind,3),'-og','linewidth',2);           % seq., long
% plot(1:100,E(ind,4),':og','linewidth',2);           % seq., long

cLeg = [{'parallel'}, {'sequential'}];%, {'seq., short'}, {'seq., long'}];
legend(cLeg,'FontSize',14);

xlabel('Instance #','FontSize',16);
ylabel('Energy','FontSize',16);

title(sprintf('64-connected grid\nSampling Empirical Distribution of CIP'),...
    'FontSize',14);