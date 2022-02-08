function allPCs = getPCs(alignedResps)

%% do PCA on the 3d matrix (cells x time x trials)
otl = size(alignedResps,1);
[i,j,k] = find(isnan(alignedResps));
nanTrials = unique(i);
for t = fliplr(nanTrials')
    alignedResps(t,:,:) = [];
end

%reshape the matrix to cells x time x trials
M = permute(alignedResps,[3 2 1]);
[i,j,k] = find(isnan(M));
nanCells = unique(i);

for c = fliplr(nanCells')
    M(c,:,:) = [];
end

% % center the data
% Mcent = zeros(size(M,1), size(M,2),size(M,3));
% for iCell = 1:size(M,1)
%     Mslice = squeeze(M(iCell,:,:));
%     MsliceMean = mean(Mslice,1);
%     Mcent(iCell,:,:) = Mslice - MsliceMean;
% end
% Mcent = M;

Mresh = reshape(M,size(M,1), size(M,2)*size(M,3));
Mresh = bsxfun(@minus, Mresh, mean(Mresh,2));
[U,S,V] = svd(Mresh,'econ');
newM = S*(V');
newMresh = reshape(newM,size(M,1),size(M,2), size(M,3));



%% reshape back into the original dimensions
%output is a matrix of size trials x time x PCs

allPCs = permute(newMresh(:,:,:),[3 2 1]);
tidx = [1 nanTrials'];
if ~isempty(nanTrials)
    for f = 1:length(tidx)-1
        apcs{f,:} = cat(1,allPCs(tidx(f):tidx(f+1)-1,:,:),nan(1,size(alignedResps,2),size(alignedResps,3)));
    end
    apcs{length(apcs)+1,:} = allPCs(tidx(end):end,:,:);
    allPCs = cell2mat(apcs);
end
