function [B,W,K] = kFactorize(X, Y, varargin)

% this function produces a kernel matrix K from a Toeplitz matrix X, and a 
% neural data matrix Y, then factorizes it into two matrices, B and W,
% for reduced-rank regression.

% X has size T x p, where T is the number of
% observations/timepoints and p is the number of predictors

% Y has size T x n, where T is the number of
% observations/timepoints and n is the number of neurons

% K has size p x n
% B has size p x r, where r is the rank. By default, rank is size p
% (i.e., nPredictors) but if you call a third variable you can set the rank
% anywhere between 1 and p

% K can be fully reconstructed by B*W, or as a reduced-rank
% reconstruction of B(:,1:t)*W(1:t,:), where t is the number of ranks to
% include.

if nargin<3
    t = size(X, 2);
else
    t = varargin{1};
    if t > size(X, 2)
        t = size(X, 2);
        warning('Rank must not exceed number of predictors; setting rank = size(X,2)')
    elseif t < 1
        t = 1;
        warning('Rank must be greater than 0; setting rank = 1')
    end
end

% Define constants
rr = size(X, 2);

full_covariance = cov([X Y]);
SSXX = full_covariance(1:rr, 1:rr);
SSYX = full_covariance((rr+1):end, 1:rr);
SSXY = full_covariance(1:rr, (rr+1):end);
SSYY = full_covariance((rr+1):end, (rr+1):end);

G = inv(SSYY);

% Define the matrix of eigen-values
[VVt,~] = eigs(sqrtm(G)*SSYX*inv(SSXX)*SSXY*sqrtm(G),t);

% Define the decomposition and mean matrices
AAt = sqrtm(inv(G))*VVt;
BBt = VVt'*sqrtm(G)*SSYX*inv(SSXX);

B = BBt';
W = AAt';
K = VVt';

end