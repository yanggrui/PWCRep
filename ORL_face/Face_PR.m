function y=Face_PR()
%% Experiments in the real-world ORL face dataset
clc
close all
clear all

load('D:\coding\Sparse_representation改进算法\Exp\ORL_face\X.mat')

%% each row of X represents a image, the size of X is 400*10304
%% Construct the graph structure
%% Randomly take Num=300 image to construct graph and use the remaining 100 image to evaluate the PSD

Rand_id=randperm(size(X,1)); %%  randon order for the 
Num=50;   %% creat graph using 50 randomly selected image
mG_id=Rand_id(1:Num);
Rem_id=Rand_id(Num+1:end);


param.k=10;
X_mG=X(mG_id,:)';   %% the feature matrix X in our paper consists of the seleted 50 images 
tic
G=gsp_nn_graph(X_mG,param); %% creat graph
t_gnn=toc







X_pre=X([Rem_id mG_id],:);
figure(1)
for j=1:16
    subplot(4,4,j)
    imshow(uint8(reshape(X_pre(j,:),112,92)))
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_face.png'])


%% the proposed method
tic
param=struct; 
param.bwd=round(G.N/20); %% use the default parameters
Partition=Multilayer_PWCRep(G,param);
t_PWCRep=toc;
fprintf('\n The total time for partitioning using PWCRep: %0.4f.\n', t_PWCRep);
L_Partition=length(Partition);
fprintf('\n The sampling size for sampling: %d.\n', L_Partition);
c_hat=zeros(L_Partition,400);
%% sample image based on the proposed method
tic
for j=1:400
    X_j=X_pre(j,:)';
    c_hat(:,j)=Sampling_PWCRep(X_j,Partition);
end
t_PCWRep1=toc %% This process takes little time

    
%% reconstruction based on the propsoed method
y_x=zeros(G.N,400);
tic
for j=1:400
    c_hat_j=c_hat(:,j);
    y_x(:,j)=Reconstruction_PWCRep(G,c_hat_j,Partition);
end
t_rec_PCWRep=toc;
fprintf('\n The total time for the reconstruction of 400 images using PWCRep: %0.4f.\n', t_rec_PCWRep);


NMSE_pro=zeros(1,400);
figure(2)
for j=1:400
    xj=reshape(y_x(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_pro(j)=norm(y_x(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_face.png'])




%% Random method
tic
K=round(G.N/20);
M=L_Partition; %% keep the same sampling size
[iset_rand,opt_weight]=SamplingMethod_Random(G,K,M);
t_samp_rand=toc;
fprintf('\n The total time for sampling using the Random method: %0.4f.\n', t_samp_rand);

xs=zeros(L_Partition,400);
%% sample image based on the Random method
tic
for j=1:400
    X_j=X_pre(j,:)';
    xs(:,j)=X_j(iset_rand);
end
t_random1=toc  %% This process takes little time

    
%% reconstruction based on the random method
xr=zeros(G.N,400);
gamma=1;
tic
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=Reconstruction_Random(G,xs_j,iset_rand,opt_weight,gamma);
end
t_rec_rand=toc;
fprintf('\n The total time for the reconstruction of 400 images using the Random method: %0.4f.\n', t_rec_rand);



NMSE_rand=zeros(1,400);
figure(3)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_rand(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_face.png'])


   

%% Ed_free method
tic
M=length(Partition);
K=round(G.N/20);
order=15;
[T,iset_Ed_free]=Sampling_Ed_free(G,M,K,order);
t_samp_Ed_free=toc;
fprintf('\n The total time for sampling using the Ed-free method: %0.4f.\n', t_samp_Ed_free);

xs=zeros(L_Partition,400);
%% sample image based on the Ed-free method
tic
for j=1:400
    X_j=X_pre(j,:)';
    xs(:,j)=X_j(iset_Ed_free);
end
t_Ed_free=toc  %% This process takes little time

    
%% reconstruction based on the Ed-free method
xr=zeros(G.N,400);
In=speye(G.N);
C=In(iset_Ed_free,:);
mu=0.01;
tic
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=(C'*C+mu*G.L)\(C'*xs_j);
end
t_rec_Ed_free=toc;
fprintf('\n The total time for the reconstruction of 400 images using the Ed-free method: %0.4f.\n', t_rec_Ed_free);

NMSE_Ed_free=zeros(1,400);
figure(4)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_Ed_free(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_face.png'])



%% EN method
tic
M=length(Partition); %% keep the same sampling size
param=struct;
param.bwd=round(G.N/20);
param.order=100;
Uest=EstimationEigenspace_EN(G,param);
iset_EN=Sampling_EN(Uest,M);
t_samp_EN=toc;

fprintf('\n The total time for sampling using the EN method: %0.4f.\n', t_samp_EN);

xs=zeros(L_Partition,400);
%% sample image based on the EN method
tic
for j=1:400
    X_j=X_pre(j,:)';
    xs(:,j)=X_j(iset_EN);
end
t_EN=toc  %% This process takes little time

    
%% reconstruction based on the EN method
xr=zeros(G.N,400);
In=speye(G.N);
C=In(iset_EN,:);
mu=0.01;
tic
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=(C'*C+mu*G.L)\(C'*xs_j);
end
t_rec_EN=toc;
fprintf('\n The total time for the reconstruction of 400 images using the EN method: %0.4f.\n', t_rec_EN);



NMSE_EN=zeros(1,400);
figure(5)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_EN(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end

saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_face.png'])




figure(8)
plot(1:16,NMSE_rand(1:16),'-o',1:16,NMSE_Ed_free(1:16),'-*', ...
    1:16,NMSE_EN(1:16),'-x',1:16,NMSE_pro(1:16),'-p',...
    'linewidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Photo index','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_face.png'])



%% image with noise
SNR=20;
X_pre_noise=X_pre;
for j=1:size(X_pre,1)
     Xj=X_pre(j,:);
     X_pre_noise(j,:)=awgn(Xj,SNR,'measured');
end
figure(11)
for j=1:16
    subplot(4,4,j)
    imshow(uint8(reshape(X_pre_noise(j,:),112,92)))
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ori_noise_face.png'])



%% the proposed method
c_hat=zeros(L_Partition,400);
%% sample image based on the proposed method
tic
for j=1:400
    X_j=X_pre_noise(j,:)';
    c_hat(:,j)=Sampling_PWCRep(X_j,Partition);
end


    
%% reconstruction based on the propsoed method
y_x=zeros(G.N,400);
for j=1:400
    c_hat_j=c_hat(:,j);
    y_x(:,j)=Reconstruction_PWCRep(G,c_hat_j,Partition);
end



NMSE_pro_noise=zeros(1,400);
figure(12)
for j=1:400
    xj=reshape(y_x(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_pro_noise(j)=norm(y_x(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\pro_rec_noise_face.png'])



%% Random method
xs=zeros(L_Partition,400);
%% sample image based on the Random method
for j=1:400
    X_j=X_pre_noise(j,:)';
    xs(:,j)=X_j(iset_rand);
end

%% reconstruction based on the Random method
xr=zeros(G.N,400);
gamma=1;
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=Reconstruction_Random(G,xs_j,iset_rand,opt_weight,gamma);
end

NMSE_rand_noise=zeros(1,400);
figure(13)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_rand_noise(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Random_rec_noise_face.png'])





%% Ed-free method
xs=zeros(L_Partition,400);
%% sample image based on the Ed-free method
for j=1:400
    X_j=X_pre_noise(j,:)';
    xs(:,j)=X_j(iset_Ed_free);
end

%% reconstruction based on the Ed-free method
xr=zeros(G.N,400);
In=speye(G.N);
C=In(iset_Ed_free,:);
mu=0.01;
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=(C'*C+mu*G.L)\(C'*xs_j);
end

NMSE_Ed_free_noise=zeros(1,400);
figure(14)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_Ed_free_noise(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Ed_free_rec_noise_face.png'])




%% EN method
xs=zeros(L_Partition,400);
%% sample image based on the EN method
for j=1:400
    X_j=X_pre_noise(j,:)';
    xs(:,j)=X_j(iset_EN);
end
    
%% reconstruction based on the EN method
xr=zeros(G.N,400);
In=speye(G.N);
C=In(iset_EN,:);
mu=0.01;
tic
for j=1:400
   xs_j=xs(:,j);
   xr(:,j)=(C'*C+mu*G.L)\(C'*xs_j);
end

NMSE_EN_noise=zeros(1,400);
figure(15)
for j=1:400
    xj=reshape(xr(:,j), 112, 92);
    if j<=16
        subplot(4,4,j)
        imshow(uint8(xj))
    end
    NMSE_EN_noise(j)=norm(xr(:,j)-X_pre(j,:)')/norm(X_pre(j,:));
end
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\EN_rec_noise_face.png'])



figure(18)
plot(1:16,NMSE_rand_noise(1:16),'-o',1:16,NMSE_Ed_free_noise(1:16),'-*', ...
    1:16,NMSE_EN_noise(1:16),'-x',1:16,NMSE_pro_noise(1:16),'-p',...
    'linewidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Photo index','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_noise_face.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_noise_face.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\ORL_face\results' ...
    '\Com_noise_face.png'])


fprintf('\n The average NMSE for different sampling methods:\n');

m_rand=mean(NMSE_rand)
m_Ed_free=mean(NMSE_Ed_free)
m_EN=mean(NMSE_EN)
m_pro=mean(NMSE_pro)

m_rand_noise=mean(NMSE_rand_noise)
m_Ed_free_noise=mean(NMSE_Ed_free_noise)
m_EN_noise=mean(NMSE_EN_noise)
m_pro_noise=mean(NMSE_pro_noise)

y=0;
end