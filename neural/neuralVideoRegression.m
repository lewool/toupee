

%% set up neural data (data matrix)

rm = neuralData.cellResps';

% whiten
rm_ms = rm - mean(rm,2);

% set up test/train neural matrices
trimTime = size(rm_ms,2) - mod(size(rm_ms,2),2);
trimNeurons = size(rm_ms,1) - mod(size(rm_ms,1),2);

F_train = rm_ms(1:2:trimNeurons,1:2:trimTime);
F_test = rm_ms(1:2:trimNeurons,2:2:trimTime);
G_train = rm_ms(2:2:trimNeurons,1:2:trimTime);
G_test = rm_ms(2:2:trimNeurons,2:2:trimTime);

%% set up video data (predictor matrix)

vidData = load('G:\Data\Video\LEW031\2020-02-03\1\2020-02-03_1_LEW031_eye_proc.mat');

% eyeData struct has the video frames, so load this
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

% pick which frames for downsampling
whichFrames = interp1(eyeData.timeAligned, 1:length(eyeData.timeAligned), neuralData.respTimes,'nearest');
whichFrames(isnan(whichFrames)) = 1;

% fetch video PCs & masks
motSVD = double(vidData.motSVD_1(whichFrames,:));
motMasks = double(vidData.motMask_reshape_1);

% split into test/train matrices
X_test = motSVD(2:2:trimTime,:);
X_train = motSVD(1:2:trimTime,:);


%% convert neural data to SVCs

% compute  using the covariance of F and G data matrices
cov = F_train * G_train';
[u,~,v] = svd(cov,'econ');

% construct U'*F_train and V'*G_train
UF = u' * F_train;
VG = v' * G_train;

% decompose the data matrices for RRR
[u_UF,s_UF,v_UF] = svd(UF','econ');
[u_VG,s_VG,v_VG] = svd(VG','econ');


%% do the regression

maxrank = 100;

for r = 1:maxrank
    
    %fetch the r-rank approximations of U'F and V'G
    UFr = u_UF(:,1:r)*s_UF(1:r,1:r)*v_UF(:,1:r)';
    VGr = u_VG(:,1:r)*s_VG(1:r,1:r)*v_VG(:,1:r)';
    
    %compute the weight matrices A and B to solve U'F = A*X and V'G = B*X
    %for the rank-r approximation
    A(:,:,r) = X_train\UFr;
    B(:,:,r) = X_train\VGr;
    
    %loss function (Frobenius norm)
    A_L(r) = norm(UF' - X_train*A(:,:,r),"fro");
    B_L(r) = norm(VG' - X_train*B(:,:,r),"fro");
end


%%
figure;
plot([A_L ;B_L]')

%%
figure
rr = 1;
for neuralPC = 1:20
    for pc = 1:size(B,1)
        wma(:,:,pc) = motMasks(:,:,pc)*B(pc,neuralPC,rr);
    end
subplot(4,5,neuralPC)
imagesc(sum(wma,3))
caxis([-.001 .001])
end


%% https://stats.stackexchange.com/questions/152517/what-is-reduced-rank-regression-all-about

%first find ordinary-least-squares coefficients
B_ols = X_train\UF';

%project this onto r principal axes of the Y data
for r = 1:60
    Br(:,:,r) = B_ols*v_UF(:,1:r)*v_UF(:,1:r)';
end

for r = 1:60
    Br_L(r) = norm(UF' - X_train*Br(:,:,r),"fro");
end