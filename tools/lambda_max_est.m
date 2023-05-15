function y=lambda_max_est(G,Partition,param)


%% maximum iteration
ite_max=15;


if ~isfield(param,'T')
    if ~isfield(G,'lmax')
        G=gsp_estimate_lmax(G);
    end
    %%  the order of the polynomial approximation
    if ~isfield(param,'order')
        param.order=30;
    end
    param.lk= gsp_fast_estimate_lk(G, param.bwd, param);
    [~, param.jch_co] = jackson_cheby_poly_coefficients(0, param.lk, [0, G.lmax], param.order);
end
x=randn(G.N,1);
x=x/norm(x);
for j=1:ite_max
    x=One_iteration(G,Partition,x,param);
end
y=norm(x).^(1/ite_max);
end

function y=One_iteration(G,Partition,x,param)
T_x=gsp_cheby_op(G, param.jch_co, x);
HT_x=zeros(size(x));
for j=1:length(Partition)
    HT_x(Partition{j})=T_x(Partition{j})-mean(T_x(Partition{j}));
end
y=gsp_cheby_op(G, param.jch_co, HT_x);
end