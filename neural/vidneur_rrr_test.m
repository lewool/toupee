%% set up data matrices

% set up neural data (time x neurons, aka long skinny)
ndm = neuralData.cellResps;

% centering
ndm_ms = ndm - mean(ndm,2);

% load vid data
vidData = load('G:\Data\Video\LEW031\2020-02-03\1\2020-02-03_1_LEW031_eye_proc.mat');

% eyeData struct has the video frames, so load this
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

% pick which frames for downsampling
whichFrames = interp1(eyeData.timeAligned, 1:length(eyeData.timeAligned), neuralData.respTimes,'nearest');
whichFrames(isnan(whichFrames)) = 1;

% fetch video PCs & masks
motSVD = double(vidData.motSVD_1(whichFrames,:)); %time x pcs
motMasks = double(vidData.motMask_reshape_1); % Y x X x pcs

%%
% plot pixel masks of the top 20 video PCs
figure;
for rr = 1:20
    subplot(4,5,rr)
    imagesc(motMasks(:,:,rr))
    caxis([-.02 .02])
end
colormap(BlueWhiteRed)

%% SVD of neural data

% decomposition
[nd_u, nd_s, nd_v] = svd(ndm_ms,'econ');

% recomposition/approximation
r = 16; %rank (arbritrary)
ndm_r = nd_u(:,1:r)*nd_s(1:r,1:r)*nd_v(:,1:r)';


%% reduced-rank regression
% https://stats.stackexchange.com/questions/152517/what-is-reduced-rank-regression-all-about

% first, regular regression
B = motSVD\ndm_ms; 


% then, decompose fitted values (Y_hat)
[u,s,v] = svd(motSVD*B,'econ');

% for each rank, project weights from full regression onto a subspace in
% Y_hat, determined by rank r (so choose a few)
for r = 1:50
    Br(:,:,r) = B*v(:,1:r)*v(:,1:r)';
    L(r) = norm(ndm_ms - motSVD*B,'fro') + norm(motSVD*B - motSVD*Br(:,:,r),'fro');
end


%% reduced-rank regression
% https://stats.stackexchange.com/questions/152517/what-is-reduced-rank-regression-all-about

% first, regular regression
B = X_train\s1'; 
Yhat = X_test*B;

% then, decompose fitted values (Y_hat)
[u,s,v] = svd(Yhat,'econ');

% for each rank, project weights from full regression onto a subspace in
% Y_hat, determined by rank r (so choose a few)
for r = 1:50
    Br(:,:,r) = B*v(:,1:r)*v(:,1:r)';
    L(r) = norm(ndm_ms - Yhat,'fro') + norm(Yhat - X_test*Br(:,:,r),'fro');
end

%% apply weights to pixel masks
pc = 5;
clear wma
for vpc = 1:500
    wma(:,:,vpc) = motMasks(:,:,vpc)*Br(vpc,pc,16);
end
subplot(1,5,pc);imagesc(sum(wma,3));caxis([-.0002 .0002])
%% plot 'eigenfaces' 
figure;
for n = 1:nc
    subplot(5,7,n)
    imagesc(sum(wma(:,:,:,n),3))
    caxis([-.0000002 .0000002])
end