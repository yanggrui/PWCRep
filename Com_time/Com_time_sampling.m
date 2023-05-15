function y=Com_time_sampling()
%% This function is to compare the computation time of different sampling methods
%% Graph size of N=10000,20000
clc
close all
clear all


%% sensor graph


Nset=[10000,20000];


L_Nset=length(Nset);

%% Sampling
t_sam_rand=zeros(L_Nset,1);
t_sam_Ed_free=zeros(L_Nset,1);
t_sam_EN=zeros(L_Nset,1);
t_sam_pro=zeros(L_Nset,1);

%% Reconstruction
t_rec_rand=zeros(L_Nset,1);
t_rec_Ed_free=zeros(L_Nset,1);
t_rec_EN=zeros(L_Nset,1);
t_rec_pro=zeros(L_Nset,1);




Kset=zeros(L_Nset,1);
Mset=zeros(L_Nset,1);


for m=1:L_Nset
    N=Nset(m);
    G=gsp_sensor(N);
    x=randn(G.N,1);
    %% There is no restriction on the original signal x
    %% Since we are only interested in sampling time and reconstruction time


    %% sampling
    tic
    param=struct;
    param.order=15;
    [Partition,~,~,~,param]=Multilayer_PWCRep(G,param);
    c=Sampling_PWCRep(x,Partition);
    t_sam_pro(m)=toc;

    K=param.bwd
    M=length(Partition)
    Kset(m)=K;
    Mset(m)=M;


    tic
    [iset_rand,opt_weight]=Sampling_Random(G,K,M);
    xs_rand=x(iset_rand);
    t_sam_rand(m)=toc;

    tic
    [T,iset_Ed_free]=Sampling_Ed_free(G,M,K,param.order);
    xs_Ed_free=x(iset_Ed_free);
    t_sam_Ed_free(m)=toc;

    tic
    param=struct;
    param.bwd=K;
    param.order=100;
    Uest=EstimationEigenspace_EN(G,param);
    iset_EN=Sampling_EN(Uest,M);
    xs_EN=x(iset_EN);
    t_sam_EN(m)=toc;

    %% reconstruction

    %% compute the computation time of reconstructing 1000 signals
    signal_num=10000;
    tic
    for i=1:signal_num
        y_x=Reconstruction_PWCRep(G,c,Partition);
    end
    t_rec_pro(m)=toc;



    tic
    for i=1:signal_num
        gamma=1;
        xr_rand=Reconstruction_Random(G,xs_rand,iset_rand,opt_weight,gamma);
    end
    t_rec_rand(m)=toc;

    tic
    In=speye(G.N);
    C=In(iset_Ed_free,:);
    mu=0.01;
    for i=1:signal_num
        xr_Ed_free=(C'*C+mu*G.L)\(C'*xs_Ed_free);
    end
    t_rec_Ed_free(m)=toc;

    tic
    In=speye(G.N);
    C=In(iset_EN,:);
    mu=0.01;
    for i=1:signal_num
        xr_EN=(C'*C+mu*G.L)\(C'*xs_EN);
    end
    t_rec_EN(m)=toc;

end
fprintf('\n The results on sensor graph of size N=%d,%d:\n',Nset(1),Nset(2));
Nset=Nset
Kset=Kset
Sampling_size_set=Mset

t_sam_rand=t_sam_rand
t_sam_Ed_free=t_sam_Ed_free
t_sam_EN=t_sam_EN
t_sam_pro=t_sam_pro

fprintf('\n The computation time of reconstructing %d signals:\n',signal_num);
t_rec_rand=t_rec_rand
t_rec_Ed_free=t_rec_Ed_free
t_rec_EN=t_rec_EN
t_rec_pro=t_rec_pro


y=0;
end