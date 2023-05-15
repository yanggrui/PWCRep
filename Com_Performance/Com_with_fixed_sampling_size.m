function y=Com_with_fixed_sampling_size()
%% Compare the reconstruction performance of different sampling methods 
%% with fiexed sampling size

clc
close all
clear all

%% sensor graph

Nset=1000:200:3000
L_Nset=length(Nset);

%% Initialization
time_rand=zeros(1,L_Nset);

time_Ed_free=zeros(1,L_Nset);

time_EN=zeros(1,L_Nset);

time_pro=zeros(1,L_Nset);

err_rand=zeros(L_Nset,4);
err_Ed_free=zeros(L_Nset,4);
err_EN=zeros(L_Nset,4);
err_pro=zeros(L_Nset,4);


K_vs_M=zeros(L_Nset,1);


for m=1:L_Nset
    N=Nset(m)
    G=gsp_sensor(N);

    tic
    param=struct;
    [Partition,~,~,~,param]=Multilayer_PWCRep(G,param);
    time_pro(m)=toc;

    K=param.bwd
    M=length(Partition)
    
    K_vs_M(m)=M/K;

    tic
    [iset_rand,opt_weight]=Sampling_Random(G,K,M);
    time_rand(m)=toc;

    tic
    [T,iset_Ed_free]=Sampling_Ed_free(G,M,K,param.order);
    time_Ed_free(m)=toc;

    tic
    param=struct;
    param.bwd=K;
    param.order=100;
    Uest=EstimationEigenspace_EN(G,param);
    iset_EN=Sampling_EN(Uest,M);
    time_EN(m)=toc;
    
    %% Initialize original signal
    G=gsp_compute_fourier_basis(G);
    x_F1=G.U(:,1:K)*rand(K,1);
    x_F2=awgn(x_F1,20,'measured');

    U_l1=greedy(G.W);
    x_F3=U_l1(:,1:K)*rand(K,1);
    x_F4=awgn(x_F3,20,'measured');

    %% sampling and reconstruction
    %% proposed method
    c_F1=Sampling_PWCRep(x_F1,Partition);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition);
    err_pro(m,1)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition);
    err_pro(m,2)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    c_F3=Sampling_PWCRep(x_F3,Partition);
    y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition);
    err_pro(m,3)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

    c_F4=Sampling_PWCRep(x_F4,Partition);
    y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition);
    err_pro(m,4)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;

    %% graph node sampling
    mu=0.01;
    %% random method
    xs_F1=x_F1(iset_rand);
    xr_F1=Reconstruction_Random(G,xs_F1,iset_rand,opt_weight,mu);
    err_rand(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

    xs_F2=x_F2(iset_rand);
    xr_F2=Reconstruction_Random(G,xs_F2,iset_rand,opt_weight,mu);
    err_rand(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

    xs_F3=x_F3(iset_rand);
    xr_F3=Reconstruction_Random(G,xs_F3,iset_rand,opt_weight,mu);
    err_rand(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

    xs_F4=x_F4(iset_rand);
    xr_F4=Reconstruction_Random(G,xs_F4,iset_rand,opt_weight,mu);
    err_rand(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;

    %% Ed-free method
    xs_F1=x_F1(iset_Ed_free);
    In=speye(G.N);
    C=In(iset_Ed_free,:);
    %     T_k=T^6;

    xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
    %     xr_F1=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F1);  %  instability
    err_Ed_free(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

    xs_F2=x_F2(iset_Ed_free);
    xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
    %     xr_F2=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F2);  %  instability
    err_Ed_free(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

    xs_F3=x_F3(iset_Ed_free);
    xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
    %     xr_F3=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F3);  %  instability
    err_Ed_free(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

    xs_F4=x_F4(iset_Ed_free);
    xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
    %     xr_F4=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F4);  %  instability
    err_Ed_free(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;


    %% EN method
    In=speye(G.N);
    C=In(iset_EN,:);

    xs_F1=x_F1(iset_EN);
    xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
    err_EN(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

    xs_F2=x_F2(iset_EN);
    xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
    err_EN(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

    xs_F3=x_F3(iset_EN);
    xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
    err_EN(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

    xs_F4=x_F4(iset_EN);
    xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
    err_EN(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;
end

figure(1)
plot(Nset,err_rand(:,1),'-o',Nset,err_Ed_free(:,1),'-*',Nset,err_EN(:,1),'-x', ...
    Nset,err_pro(:,1),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F1_sensor_fixed_size.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F1_sensor_fixed_size.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F1_sensor_fixed_size.png'])

figure(2)
plot(Nset,err_rand(:,2),'-o',Nset,err_Ed_free(:,2),'-*',Nset,err_EN(:,2),'-x', ...
    Nset,err_pro(:,2),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Graph size','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F2_sensor_fixed_size.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F2_sensor_fixed_size.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F2_sensor_fixed_size.png'])


figure(3)
plot(Nset,err_rand(:,3),'-o',Nset,err_Ed_free(:,3),'-*',Nset,err_EN(:,3),'-x', ...
    Nset,err_pro(:,3),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F3_sensor_fixed_size.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F3_sensor_fixed_size.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F3_sensor_fixed_size.png'])


figure(4)
plot(Nset,err_rand(:,4),'-o',Nset,err_Ed_free(:,4),'-*',Nset,err_EN(:,4),'-x', ...
    Nset,err_pro(:,4),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Graph size； N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F4_sensor_fixed_size.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F4_sensor_fixed_size.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Com_F4_sensor_fixed_size.png'])


figure(33)
plot(Nset,K_vs_M,'LineWidth',1.5)
xlabel('Graph size','Fontsize',12)
ylabel('Factor','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Factor_sensor_fixed_size.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Factor_sensor_fixed_size.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_with_fixed_sampling_size\results\Factor_sensor_fixed_size.png'])




y=0;
end