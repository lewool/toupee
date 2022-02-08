function neuralData = getCellPCs(neuralData)

%INPUT: cellResps = one 2D matrix (cells x time)
%OUTPUT: one 2D matrix (time x PCs)

if size(neuralData.cellResps,1) > size(neuralData.cellResps,2)
    X = neuralData.cellResps';
else
    X = neuralData.cellResps;
end
X = zscore(X);

[U,S,V] = svd(X,'econ');

[~, max_Idx] = max(abs(V(:,1)));
    if sign(V(max_Idx,1)) < 0
        V = -V;
        U = -U;
    end
    
neuralData.PCs = V;