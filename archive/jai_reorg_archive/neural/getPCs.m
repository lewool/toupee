function allPCs = getPCs(alignedResps)

%% do PCA on the 3d matrix (cells x time x trials)

for a = 1:length(alignedResps)
%reshape the matrix to cells x time x trials
M = permute(alignedResps{a},[3 2 1]);
[i,~] = find(isnan(M));
nanCells = unique(i);

for c = fliplr(nanCells')
    M(c,:,:) = [];
end

% center the data
Mcent = zeros(size(M,1), size(M,2),size(M,3));
for iCell = 1:size(M,1)
    Mslice = squeeze(M(iCell,:,:));
    MsliceMean = mean(Mslice,1);
    Mcent(iCell,:,:) = Mslice - MsliceMean;
end
% Mcent = M;

Mresh = reshape(Mcent,size(Mcent,1), size(Mcent,2)*size(Mcent,3));
[U,S,V] = svd(Mresh,'econ');
newM = S*abs(V');
newMresh = reshape(newM,size(M,1),size(M,2), size(M,3));

%% reshape back into the original dimensions
%output is a matrix of size trials x time x PCs

allPCs{a} = permute(newMresh(1:numPCs,:,:),[3 2 1]);

end