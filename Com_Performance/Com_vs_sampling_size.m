function y=Com_vs_sampling_size()
%% Compare the reconstruction performance of different sampling methods 
%% under different sampling size


%% sensor graph
close all
clear all


%% sensor graph
N=1000;
G=gsp_sensor(N);
G=gsp_compute_fourier_basis(G);
K=round(G.N/20);

%% Initialize original signal
x_F1=G.U(:,1:K)*rand(K,1);
x_F2=awgn(x_F1,20,'measured');

U_l1=greedy(G.W);
x_F3=U_l1(:,1:K)*rand(K,1);
x_F4=awgn(x_F3,20,'measured');


tic
param=struct;
param.bwd=K;
[~,~,All_Partition,param]=Multilayer_PWCRep_J(G,param);
time_pro=toc;

L_all_Partition=length(All_Partition);




m=0;
for j=1:length(All_Partition)

    Partition=All_Partition{j};
    M=length(Partition)
    if M<5*K   %% we focus on the case where sampling size is less than 5*K;
        m=m+1;
        Jset(m)=M;


        tic
        [iset_rand,opt_weight]=Sampling_Random(G,K,M);
        time_rand(m)=toc;

        tic
        M_Ed_free=M;
        [T,iset_Ed_free]=Sampling_Ed_free(G,M_Ed_free,K,param.order);
        time_Ed_free(m)=toc;

        tic
        param=struct;
        param.bwd=K;
        param.order=100;
        Uest=EstimationEigenspace_EN(G,param);
        iset_EN=Sampling_EN(Uest,M);
        time_EN(m)=toc;



        %% sampling and reconstruction
        %% proposed method
        c_F1=Sampling_PWCRep(x_F1,Partition);
        y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition);
        NMSE_pro(m,1)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

        c_F2=Sampling_PWCRep(x_F2,Partition);
        y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition);
        NMSE_pro(m,2)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

        c_F3=Sampling_PWCRep(x_F3,Partition);
        y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition);
        NMSE_pro(m,3)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

        c_F4=Sampling_PWCRep(x_F4,Partition);
        y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition);
        NMSE_pro(m,4)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;

        %% graph node sampling
        mu=0.01;
        %% random method
        xs_F1=x_F1(iset_rand);
        xr_F1=Reconstruction_Random(G,xs_F1,iset_rand,opt_weight,mu);
        NMSE_rand(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_rand);
        xr_F2=Reconstruction_Random(G,xs_F2,iset_rand,opt_weight,mu);
        NMSE_rand(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_rand);
        xr_F3=Reconstruction_Random(G,xs_F3,iset_rand,opt_weight,mu);
        NMSE_rand(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_rand);
        xr_F4=Reconstruction_Random(G,xs_F4,iset_rand,opt_weight,mu);
        NMSE_rand(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;

        %% Ed-free method
        xs_F1=x_F1(iset_Ed_free);
        In=speye(G.N);
        C=In(iset_Ed_free,:);
        %     T_k=T^6;

        xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
        %     xr_F1=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F1);  % instability
        NMSE_Ed_free(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_Ed_free);
        xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
        %     xr_F2=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F2);  % instability
        NMSE_Ed_free(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_Ed_free);
        xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
        %     xr_F3=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F3);  % instability
        NMSE_Ed_free(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_Ed_free);
        xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
        %     xr_F4=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F4);   % instability
        NMSE_Ed_free(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;


        %% EN method
        In=speye(G.N);
        C=In(iset_EN,:);

        xs_F1=x_F1(iset_EN);
        xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
        NMSE_EN(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_EN);
        xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
        NMSE_EN(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_EN);
        xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
        NMSE_EN(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_EN);
        xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
        NMSE_EN(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;
    end
end

figure(1)
plot(Jset,NMSE_rand(:,1),'-o',Jset,NMSE_Ed_free(:,1),'-*',Jset,NMSE_EN(:,1),'-x', ...
    Jset,NMSE_pro(:,1),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_sensor.png'])

figure(2)
plot(Jset,NMSE_rand(:,2),'-o',Jset,NMSE_Ed_free(:,2),'-*',Jset,NMSE_EN(:,2),'-x', ...
    Jset,NMSE_pro(:,2),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_sensor.png'])


figure(3)
plot(Jset,NMSE_rand(:,3),'-o',Jset,NMSE_Ed_free(:,3),'-*',Jset,NMSE_EN(:,3),'-x', ...
    Jset,NMSE_pro(:,3),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_sensor.png'])

figure(4)
plot(Jset,NMSE_rand(:,4),'-o',Jset,NMSE_Ed_free(:,4),'-*',Jset,NMSE_EN(:,4),'-x', ...
    Jset,NMSE_pro(:,4),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_sensor.png'])





%% community graph
clear all


N=1000;
G=gsp_community(N);
G=gsp_compute_fourier_basis(G);
K=round(G.N/20);
x_F1=G.U(:,1:K)*rand(K,1);
x_F2=awgn(x_F1,20,'measured');

U_l1=greedy(G.W);
x_F3=U_l1(:,1:K)*rand(K,1);
x_F4=awgn(x_F3,20,'measured');


tic
param=struct;
param.bwd=K;
[~,~,All_Partition,param]=Multilayer_PWCRep_J(G,param);
time_pro=toc;




m=0;
for j=1:length(All_Partition)

    Partition=All_Partition{j};
    M=length(Partition)
    if M<5*K   %% we focus on the case where sampling size is less than 5*K;
        m=m+1;
        Jset(m)=M;


        tic
        [iset_rand,opt_weight]=Sampling_Random(G,K,M);
        time_rand(m)=toc;

        tic
        M_Ed_free=M;
        [T,iset_Ed_free]=Sampling_Ed_free(G,M_Ed_free,K,param.order);
        time_Ed_free(m)=toc;

        tic
        param=struct;
        param.bwd=K;
        param.order=100;
        Uest=EstimationEigenspace_EN(G,param);
        iset_EN=Sampling_EN(Uest,M);
        time_EN(m)=toc;



        %% sampling and reconstruction
        %% proposed method
        c_F1=Sampling_PWCRep(x_F1,Partition);
        y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition);
        NMSE_pro(m,1)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

        c_F2=Sampling_PWCRep(x_F2,Partition);
        y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition);
        NMSE_pro(m,2)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

        c_F3=Sampling_PWCRep(x_F3,Partition);
        y_x_F3=Reconstruction_PWCRep(G,c_F3,Partition);
        NMSE_pro(m,3)=norm(y_x_F3-x_F3)^2/norm(x_F3)^2;

        c_F4=Sampling_PWCRep(x_F4,Partition);
        y_x_F4=Reconstruction_PWCRep(G,c_F4,Partition);
        NMSE_pro(m,4)=norm(y_x_F4-x_F3)^2/norm(x_F3)^2;

        %% graph node sampling
        mu=0.01;
        %% random method
        xs_F1=x_F1(iset_rand);
        xr_F1=Reconstruction_Random(G,xs_F1,iset_rand,opt_weight,mu);
        NMSE_rand(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_rand);
        xr_F2=Reconstruction_Random(G,xs_F2,iset_rand,opt_weight,mu);
        NMSE_rand(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_rand);
        xr_F3=Reconstruction_Random(G,xs_F3,iset_rand,opt_weight,mu);
        NMSE_rand(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_rand);
        xr_F4=Reconstruction_Random(G,xs_F4,iset_rand,opt_weight,mu);
        NMSE_rand(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;

        %% Ed-free method
        In=speye(G.N);
        C=In(iset_Ed_free,:);
        %     T_k=T^6;

        xs_F1=x_F1(iset_Ed_free);
        xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
        %     xr_F1=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F1);  % instability
        NMSE_Ed_free(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_Ed_free);
        xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
        %     xr_F2=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F2);  % instability
        NMSE_Ed_free(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_Ed_free);
        xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
        %     xr_F3=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F3);  % instability
        NMSE_Ed_free(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_Ed_free);
        xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
        %     xr_F4=T_k(:,iset_Ed_free)*(T_k(iset_Ed_free,iset_Ed_free)\xs_F4);   % instability
        NMSE_Ed_free(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;


        %% EN method
        In=speye(G.N);
        C=In(iset_EN,:);

        xs_F1=x_F1(iset_EN);
        xr_F1=(C'*C+mu*G.L)\(C'*xs_F1);
        NMSE_EN(m,1)=norm(xr_F1-x_F1)^2/norm(x_F1)^2;

        xs_F2=x_F2(iset_EN);
        xr_F2=(C'*C+mu*G.L)\(C'*xs_F2);
        NMSE_EN(m,2)=norm(xr_F2-x_F1)^2/norm(x_F1)^2;

        xs_F3=x_F3(iset_EN);
        xr_F3=(C'*C+mu*G.L)\(C'*xs_F3);
        NMSE_EN(m,3)=norm(xr_F3-x_F3)^2/norm(x_F3)^2;

        xs_F4=x_F4(iset_EN);
        xr_F4=(C'*C+mu*G.L)\(C'*xs_F4);
        NMSE_EN(m,4)=norm(xr_F4-x_F3)^2/norm(x_F3)^2;
    end
end

figure(5)
plot(Jset,NMSE_rand(:,1),'-o',Jset,NMSE_Ed_free(:,1),'-*',Jset,NMSE_EN(:,1),'-x', ...
    Jset,NMSE_pro(:,1),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F1_community.png'])

figure(6)
plot(Jset,NMSE_rand(:,2),'-o',Jset,NMSE_Ed_free(:,2),'-*',Jset,NMSE_EN(:,2),'-x', ...
    Jset,NMSE_pro(:,2),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F2_community.png'])


figure(7)
plot(Jset,NMSE_rand(:,3),'-o',Jset,NMSE_Ed_free(:,3),'-*',Jset,NMSE_EN(:,3),'-x', ...
    Jset,NMSE_pro(:,3),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F3_community.png'])

figure(8)
plot(Jset,NMSE_rand(:,4),'-o',Jset,NMSE_Ed_free(:,4),'-*',Jset,NMSE_EN(:,4),'-x', ...
    Jset,NMSE_pro(:,4),'-p','LineWidth',1.5)
l=legend('\it{Random}','\it{Ed-free}','\it{EN}','\it{Prop.}');
set(l,'Fontsize',12,'Location','best')
xlabel('Sampling size: J','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_performance' ...
    '\Com_vs_sampling_size\results\Com_F4_community.png'])
y=0;
end