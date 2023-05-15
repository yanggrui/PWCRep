function y=Com_Single_VS_Multilater_PCWRep()
%% the proposed single-layer algorithm VS. the multi-layer algorithm
clc
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   sensor graph   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=1000;
G=gsp_sensor(N);

tic
param=struct
[~,epsilon,All_Partition,~,param]=Multilayer_PWCRep(G,param);
t_PK=toc
fprintf('* Time to perform the Multilayer algorithm: %d sec (Graph size=%d,Bandwidth=%d).\n',t_PK,G.N,param.bwd)
epsilon

if ~isfield(G,'e')
    G=gsp_compute_fourier_basis(G);
end
K=param.bwd;
%% Initialize original signal
x_F1=G.U(:,1:K)*rand(K,1);
x_F2=awgn(x_F1,20,'measured');

U_l1=greedy(G.W);
x_F3=U_l1(:,1:K)*rand(K,1);
x_F4=awgn(x_F3,20,'measured');


L_Jset=length(All_Partition) %% the number of node partition constructed in Algorithm 2
Jset_multi=zeros(1,L_Jset);
Jset_single=zeros(1,L_Jset);

err_multi=zeros(4,L_Jset);
err_single=zeros(4,L_Jset);

for m=1:length(epsilon)
    epsilon_m=epsilon(m);
    param=struct;
    param.epsilon=epsilon_m;
    %% Partition based on Algorithm 1 using epsilon_m
    tic
    Partition_single=Singlelayer_PWCRep(G,param);
    t_single=toc;
    fprintf('* Time to perform the Single layer algorithm: %d sec (Graph size=%d,Bandwidth=%d).\n',t_single,G.N,K)

    Jset_single(m)=length(Partition_single);


    %% representation error of single-layer algorithm
    c_F1=Sampling_PWCRep(x_F1,Partition_single);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition_single);
    err_single(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition_single);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition_single);
    err_single(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    c_F3=Sampling_PWCRep(x_F3,Partition_single);
    y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition_single);
    err_single(3,m)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

    c_F4=Sampling_PWCRep(x_F4,Partition_single);
    y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition_single);
    err_single(4,m)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;


    %% the m-th node partition constructed in Algorithm 2
    Partition_multi=All_Partition{m};
    Jset_multi(m)=length(Partition_multi);

    %% representation error of multi-layer algorithm
    c_F1=Sampling_PWCRep(x_F1,Partition_multi);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition_multi);
    err_multi(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition_multi);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition_multi);
    err_multi(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    c_F3=Sampling_PWCRep(x_F3,Partition_multi);
    y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition_multi);
    err_multi(3,m)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

    c_F4=Sampling_PWCRep(x_F4,Partition_multi);
    y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition_multi);
    err_multi(4,m)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;

end
Jset_single_sensor=Jset_single
err_single_sensor=err_single

Jset_multi_sensor=Jset_multi
err_multi_sensor=err_multi

figure(1)
plot(epsilon,err_single(1,:),'-o',epsilon,err_multi(1,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_sensor.png'])

figure(2)
plot(epsilon,err_single(2,:),'-o',epsilon,err_multi(2,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_sensor.png'])

figure(3)
plot(epsilon,err_single(3,:),'-o',epsilon,err_multi(3,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_sensor.png'])
figure(4)
plot(epsilon,err_single(4,:),'-o',epsilon,err_multi(4,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_sensor.png'])

figure(11)
plot(epsilon,Jset_single,'-o',epsilon,Jset_multi,'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('\epsilon','Fontsize',12)
ylabel('Piece number','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_sensor.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   community graph   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

N=1000;
G=gsp_community(N);

tic
param=struct
[~,epsilon,All_Partition,~,param]=Multilayer_PWCRep(G,param);
t_PK=toc
fprintf('* Time to perform the Multilayer algorithm: %d sec (Graph size=%d,Bandwidth=%d).\n',t_PK,G.N,param.bwd)
epsilon

if ~isfield(G,'e')
    G=gsp_compute_fourier_basis(G);
end
K=param.bwd;
%% Initialize original signal
x_F1=G.U(:,1:K)*rand(K,1);
x_F2=awgn(x_F1,20,'measured');

U_l1=greedy(G.W);
x_F3=U_l1(:,1:K)*rand(K,1);
x_F4=awgn(x_F3,20,'measured');


L_Jset=length(All_Partition) %% the number of node partition constructed in Algorithm 2
Jset_multi=zeros(1,L_Jset);
Jset_single=zeros(1,L_Jset);

err_multi=zeros(4,L_Jset);
err_single=zeros(4,L_Jset);

for m=1:length(epsilon)
    epsilon_m=epsilon(m);
    param=struct;
    param.epsilon=epsilon_m;
    %% Partition based on Algorithm 1 using epsilon_m
    tic
    Partition_single=Singlelayer_PWCRep(G,param);
    t_single=toc;
    fprintf('* Time to perform the Single layer algorithm: %d sec (Graph size=%d,Bandwidth=%d).\n',t_single,G.N,K)

    Jset_single(m)=length(Partition_single);


    %% representation error of single-layer algorithm
    c_F1=Sampling_PWCRep(x_F1,Partition_single);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition_single);
    err_single(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition_single);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition_single);
    err_single(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    c_F3=Sampling_PWCRep(x_F3,Partition_single);
    y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition_single);
    err_single(3,m)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

    c_F4=Sampling_PWCRep(x_F4,Partition_single);
    y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition_single);
    err_single(4,m)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;


    %% the m-th node partition constructed in Algorithm 2
    Partition_multi=All_Partition{m};
    Jset_multi(m)=length(Partition_multi);

    %% representation error of multi-layer algorithm
    c_F1=Sampling_PWCRep(x_F1,Partition_multi);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition_multi);
    err_multi(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition_multi);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition_multi);
    err_multi(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    c_F3=Sampling_PWCRep(x_F3,Partition_multi);
    y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition_multi);
    err_multi(3,m)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

    c_F4=Sampling_PWCRep(x_F4,Partition_multi);
    y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition_multi);
    err_multi(4,m)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;

end
Jset_single_community=Jset_single
err_single_community=err_single

Jset_multi_community=Jset_multi
err_multi_community=err_multi

figure(5)
plot(epsilon,err_single(1,:),'-o',epsilon,err_multi(1,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F1_community.png'])

figure(6)
plot(epsilon,err_single(2,:),'-o',epsilon,err_multi(2,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F2_community.png'])

figure(7)
plot(epsilon,err_single(3,:),'-o',epsilon,err_multi(3,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F3_community.png'])
figure(8)
plot(epsilon,err_single(4,:),'-o',epsilon,err_multi(4,:),'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northwest')
xlabel('\epsilon','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\F4_community.png'])

figure(22)
plot(epsilon,Jset_single,'-o',epsilon,Jset_multi,'-*','linewidth',1.5)
l=legend('\bf{SLPR}','\bf{MLPR}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('\epsilon','Fontsize',12)
ylabel('Piece number','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\MLPRvsSLPR\results' ...
    '\J_community.png'])

y=0;
end





