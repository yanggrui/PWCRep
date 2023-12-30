function y=Com_diff_Rep()
%% This function is to compare different piecewise representation methods.


close all
clear all


%%  sensor graph

Nset=1000:200:3000
L_Nset=length(Nset);

%% Initialize the representation error
err_prop=zeros(2,L_Nset); err_Rand_assign=zeros(2,L_Nset);
err_Rand_Partition=zeros(2,L_Nset);err_l1=zeros(2,L_Nset);
%% 


for m=1:L_Nset
    N=Nset(m);
    G=gsp_sensor(N);


    param=struct;
    [Partition,~,~,~,~]=Multilayer_PWCRep(G,param);

    








    %% Initialize original signal
    G=gsp_compute_fourier_basis(G); K=round(G.N/20);
    x_F1=G.U(:,1:K)*rand(K,1);
    x_F2=awgn(x_F1,20,'measured');

    %% proposed 
    c_F1=Sampling_PWCRep(x_F1,Partition);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition);
    err_prop(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition);
    err_prop(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    %% rand assign
    y_F1=Rep_rand_assign(x_F1,Partition);
    err_Rand_assign(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

    y_F2=Rep_rand_assign(x_F2,Partition);
    err_Rand_assign(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;

   %% Rand Partition
   y_F1=Rep_rand_Partition(x_F1,Partition);
   err_Rand_Partition(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

   y_F2=Rep_rand_Partition(x_F2,Partition);
   err_Rand_Partition(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;

    
   %% greedy
   U_l1=greedy(G.W);
   y_F1=U_l1(:,1:K)*U_l1(:,1:K)'*x_F1;
   err_l1(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

   y_F2=U_l1(:,1:K)*U_l1(:,1:K)'*x_F2;
   err_l1(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;
end

figure(1)
semilogy(Nset,err_Rand_assign(1,:),'-o',Nset,err_Rand_Partition(1,:),'-*', ...
    Nset,err_l1(1,:),'-x', Nset,err_prop(1,:),'-p','LineWidth',1.5)
l=legend('\it{Ra}','\it{Rp}','\it{greedy}','\it{Prop.}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_sensor.png'])

figure(2)
semilogy(Nset,err_Rand_assign(2,:),'-o',Nset,err_Rand_Partition(2,:),'-*', ...
    Nset,err_l1(2,:),'-x', Nset,err_prop(2,:),'-p','LineWidth',1.5)
l=legend('\it{Ra}','\it{Rp}','\it{greedy}','\it{Prop.}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_sensor.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_sensor.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_sensor.png'])





%%  community graph

Nset=1000:200:3000
L_Nset=length(Nset);

%% Initialize the representation error
err_prop=zeros(2,L_Nset); err_Rand_assign=zeros(2,L_Nset);
err_Rand_Partition=zeros(2,L_Nset);err_l1=zeros(2,L_Nset);
%% 


for m=1:L_Nset
    N=Nset(m);
    G=gsp_community(N);


    param=struct;
    [Partition,~,~,~,~]=Multilayer_PWCRep(G,param);

    








    %% Initialize original signal
    G=gsp_compute_fourier_basis(G); K=round(G.N/20);
    x_F1=G.U(:,1:K)*rand(K,1);
    x_F2=awgn(x_F1,20,'measured');

    %% proposed 
    c_F1=Sampling_PWCRep(x_F1,Partition);
    y_x_F1=Reconstruction_PWCRep(G,c_F1,Partition);
    err_prop(1,m)=norm(y_x_F1-x_F1)^2/norm(x_F1)^2;

    c_F2=Sampling_PWCRep(x_F2,Partition);
    y_x_F2=Reconstruction_PWCRep(G,c_F2,Partition);
    err_prop(2,m)=norm(y_x_F2-x_F1)^2/norm(x_F1)^2;

    %% rand assign
    y_F1=Rep_rand_assign(x_F1,Partition);
    err_Rand_assign(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

    y_F2=Rep_rand_assign(x_F2,Partition);
    err_Rand_assign(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;

   %% Rand Partition
   y_F1=Rep_rand_Partition(x_F1,Partition);
   err_Rand_Partition(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

   y_F2=Rep_rand_Partition(x_F2,Partition);
   err_Rand_Partition(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;

    
   %% greedy
   U_l1=greedy(G.W);
   y_F1=U_l1(:,1:K)*U_l1(:,1:K)'*x_F1;
   err_l1(1,m)=norm(y_F1-x_F1)^2/norm(x_F1)^2;

   y_F2=U_l1(:,1:K)*U_l1(:,1:K)'*x_F2;
   err_l1(2,m)=norm(y_F2-x_F1)^2/norm(x_F1)^2;
end

figure(3)
semilogy(Nset,err_Rand_assign(1,:),'-o',Nset,err_Rand_Partition(1,:),'-*', ...
    Nset,err_l1(1,:),'-x', Nset,err_prop(1,:),'-p','LineWidth',1.5)
l=legend('\it{Ra}','\it{Rp}','\it{greedy}','\it{Prop.}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_Community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_Community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F1_Community.png'])

figure(4)
semilogy(Nset,err_Rand_assign(2,:),'-o',Nset,err_Rand_Partition(2,:),'-*', ...
    Nset,err_l1(2,:),'-x', Nset,err_prop(2,:),'-p','LineWidth',1.5)
l=legend('\it{Ra}','\it{Rp}','\it{greedy}','\it{Prop.}');
set(l,'Fontsize',12,'Location','northeast')
xlabel('Graph size: N','Fontsize',12)
ylabel('\bf{NMSE}','Fontsize',12)
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_Community.fig'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_Community.jpg'])
saveas(gcf,['D:\coding\Sparse_representation改进算法\Exp\Com_diff_Rep' ...
    '\results\Com_diff_Rep_F2_Community.png'])




y=0;
end



function y=Rep_rand_assign(x,Partition)
y=zeros(size(x));
L_part=length(Partition);
for j=1:L_part
    id_assign=randi(length(Partition{j}));
    y(Partition{j})=x(Partition{j}(id_assign));
end
end

function y=Rep_rand_Partition(x,Partition)
y=zeros(size(x));
N=length(x);
L_part=length(Partition);
Rem_nodes=1:N;
for j=1:L_part
    L_Rem=length(Rem_nodes);
    size_part_j=length(Partition{j});
    id_part_j=datasample(1:L_Rem, size_part_j, ...
            'Replace', false, 'Weights', ones(1,L_Rem)/L_Rem);
    part_j=Rem_nodes(id_part_j);
    Rem_nodes(id_part_j)=[];
    y(part_j)=mean(x(part_j));
end

end

