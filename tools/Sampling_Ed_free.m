function [T,selected_nodes]=Sampling_Ed_free(G,M,K,p)
%% This function perform the node selection using the localization operator.
%% M indicates the size of sampling set.



if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
end

param.order=p;

lk= gsp_fast_estimate_lk(G, K, param);
[~,jch_co] = jackson_cheby_poly_coefficients(0, lk, [0, G.lmax], param.order);
T= gsp_cheby_op(G,jch_co, speye(G.N));




abs_T=abs(T);
sensor = selection([],abs_T);
selected_nodes = sensor;

for i=1:M-1
    sensor = selection(selected_nodes,abs_T);
    selected_nodes = [selected_nodes sensor];
end

end
function sensor = selection(selected,T_g_tmp)

if ~isempty(selected)
    T=sum(T_g_tmp(:,selected),2);
    T2=mean(T)-T;
    T2(logical(T2<0))=0;
    T_g=(T_g_tmp)*(T2);   
else
     T_g=sum(T_g_tmp);
end
    
T_g(selected) = 0;
[~,sensor]  =max(T_g);

end