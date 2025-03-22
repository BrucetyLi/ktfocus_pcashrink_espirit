clear;clc;close;

addpath('data')
addpath(genpath('bin'))
%% Loading data
key = 'NIST39743_12x';
disp(key)
data_all = load(key);
kspace = squeeze(data_all.(char(fieldnames(data_all))));


%% Basic Settings
% iteration parameters
opt = optset('n');
opt.th = 0;

opt.T1=@(T2)1.125;
opt.esp = 6.312/1000;
opt.etl = 30;
opt.RFr.angle = 110;
opt.FA_array = [175/110,145/110,ones(1,28)];

opt.num_low_phase = 6;

% parameters for ktFocuss
opt.ktfocuss.lambda1 = 0.001;
opt.ktfocuss.lambda2 = 0.5;
opt.ktfocuss.inner_iter = 40;
opt.ktfocuss.outer_iter = 2;

% parameters for PCA
lambda3 = 0.001;      
iteration = 40;   

% espirit parameters
opt.es.ksize = [max(ceil(opt.num_low_phase/4),2),max(ceil(opt.num_low_phase/4),2)];
opt.es.eigThresh_1 = 0.02;
opt.es.eigThresh_2 = 0.95;

%% pretreatment
if length(size(kspace)) ==4
    [opt.ny,opt.nx,opt.nt,opt.nc] = size(kspace);
else
    [opt.ny,opt.nx,opt.sl,opt.nt,opt.nc] = size(kspace);
    selection = 1;  
    kspace = squeeze(kspace(:,:,selection,:,:));
end
opt.mask = abs(squeeze(kspace(:,:,:,1)))>0;

% scaling  
scale = max(max(max(sqrt(sum(abs(ifftc(ifftc(kspace,2),1)).^2,4)))));
 y = kspace/scale;

% Sense Coil Sensitivity
% opt.maps = Sense(y,opt.num_low_phase);

% espirit
calib = crop(squeeze(sum(y,3)/size(y,3)),[opt.num_low_phase,opt.num_low_phase,opt.nc]);
[k,SE] = dat2Kernel(calib,opt.es.ksize);
idx = find(SE >= SE(1)*opt.es.eigThresh_1, 1, 'last' );
[ME,WE] = kernelEig(k(:,:,:,1:idx),[opt.ny,opt.nx]);
opt.maps = ME(:,:,:,end).*repmat(WE(:,:,end)>opt.es.eigThresh_2,[1,1,opt.nc]);

%% ktfocuss 
row = opt.ny;
col = opt.nx;
echo = opt.nt;
channel = opt.nc;

mat = ifftc(y,2);
low = mat;
low(1:floor(0.5*row-opt.num_low_phase/2),:,:,:) = 0;
low(ceil(0.5*row+opt.num_low_phase/2+1):end,:,:,:) = 0;
A = @(x,mask)  fftc(x,1).*mask;
AT = @(x,mask) ifftc(x.*mask,1);
channel_kt = zeros(size(mat));
for j = 1:channel
disp(['focuss channel  : ',num2str(j)]);
Y = squeeze(mat(:,:,:,j));
Low_Y = squeeze(low(:,:,:,j));
FOCUSS = KTFOCUSS(A,AT,Y,Low_Y,opt.mask,opt.ktfocuss.lambda2,opt.ktfocuss.lambda1, opt.ktfocuss.inner_iter,opt.ktfocuss.outer_iter);
channel_kt(:,:,:,j) = FOCUSS;
disp('*****************************************************')
end

%% epgfit for ktfocuss
% P = reshape(channel_kt,[row*col,echo,channel]);
% Z = sqrt(sum(abs(P).^2,3));
% 
% Z(abs(Z(:,1))<=opt.th,:) = 0;
% T2 = zeros(row*col,1);
% B1 = zeros(row*col,1);
% I = zeros(row*col,1);
% parfor j= 1:row*col
%     temp = Z(j,:);
%      if abs(temp(1))==0
%         T2(j)=0; B1(j)=1;I(j)=0;
%     else
%     [T2(j),B1(j),I(j)] = IniFit(squeeze(Z(j,:)),opt);
%      end
% end
% Tm = [T2,B1];
% 
% TIZ = [Tm,I*scale,Z*scale];
% save(fullfile('results','es','ktfocuss',key),'TIZ')

%% pca_shrink   
focuss_results =  PImerge(opt.maps,channel_kt);
[V,D]= eig(cov(reshape(focuss_results,[row*col,echo])));
W = initial(focuss_results,V);
T = @(x,option) forward(x,opt.mask,V,opt.maps,option);
for j = 1:iteration
W = reshape(IST(W(:),T,mat(:),lambda3),size(opt.mask));
disp(['iteration  : ',num2str(j)]);
end
pca_results = reshape(reshape(W,[row*col,echo])*conj(V'),[row,col,echo]);
save(fullfile('results',key),'pca_results')

%% epgfit for pca_shrink   
% Z = abs(reshape(final_results,[row*col,echo]));
% Z(abs(Z(:,1))<=opt.th,:) = 0;
% T2 = zeros(row*col,1);
% B1 = zeros(row*col,1);
% I = zeros(row*col,1);
% parfor j= 1:row*col
%     temp = Z(j,:);
%      if abs(temp(1))==0
%         T2(j)=0; B1(j)=1;I(j)=0;
%     else
% [T2(j),B1(j),I(j)] = IniFit(squeeze(Z(j,:)),opt);
%      end
% end
% Tm = [T2,B1];
% 
% TIZ = [Tm,I*scale,Z*scale];
% save(fullfile('results',key),'TIZ')










