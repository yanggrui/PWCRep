function xr=Reconstruction_Random(G,x_sampling,iset_sampling,opt_weight,gamma)
%% Reconstruction for random sampling method M4: without replacement

P=sparse(1:G.N, 1:G.N, 1./sqrt(opt_weight), G.N, G.N);  %%注意此处对应文献中的P^-1/2;

% Sampling matrix (without replacement)
Samplingsize=length(iset_sampling);
M = sparse(1:Samplingsize, iset_sampling, 1, Samplingsize, ...
    G.N);

% Prepare matrices for reconstruction

MP = M*P; MtM = MP'*MP;

% Reconstruction
% regularisation parameters

xr = (MtM + gamma*G.L)\(P*P*M'*x_sampling);





end