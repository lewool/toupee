function [svdResps, svdVals] = getSVD(alignedResps, c)
%% overview

% computes the singular value decomposition for a given set of time-aligned
% cell responses.
% alignedResps: 3D response array of size nTrials x nTimepoints x nCells
% c: how many components do you want to take? default = 1
% 2020 Jan 27 LEW

%% check input args

if nargin < 2
    c = 1;
end

%% initialize matrices

U = nan(size(alignedResps,1));
V = nan(size(alignedResps,2));
S = nan(length(U),length(V));

V_all = nan(size(alignedResps,2),c,size(alignedResps,3));
U_all = nan(size(alignedResps,1),c,size(alignedResps,3));
S_all = nan(c,size(alignedResps,3));

svdResps = nan(size(alignedResps));

%% SVD of all cell responses


for iCell = 1:size(alignedResps,3)
    try
        [U,S,V] = svd(alignedResps(:,:,iCell));
    catch
        U = nan(size(alignedResps,1));
        V = nan(size(alignedResps,2));
        S = nan(length(U),length(V));
    end
    
    [~, max_Idx] = max(abs(V(:,1)));
    if sign(V(max_Idx,1)) < 0
        V = -V;
        U = -U;
    end
    
    U_all(:,:, iCell) = U(:,1:c);
    V_all(:,:,iCell) = V(:,1:c);
    S_all(:,iCell) = diag(S(1:c,1:c)); 
    
    svdResps(:,:,iCell) = U(:,1:c)*S(1:c,1:c)*V(:,1:c)';    
end

if size(U_all,2) == 1
    U_all = squeeze(U_all);
    V_all = squeeze(V_all);
end

%% store outputs

svdVals.U = U_all;
svdVals.V = V_all;
svdVals.S = S_all;



