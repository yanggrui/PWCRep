function [Partition,varargout]=Singlelayer_PWCRep(G,param)


%% Initialize epsilon
if ~isfield(param,'epsilon')
    param.epsilon=1*1e-1;
end



if ~isfield(param,'dis_node')
    if ~isfield(param,'T')
        if ~isfield(G,'lmax')
            G=gsp_estimate_lmax(G);
        end
        %% bandwidth
        if ~isfield(param,'bwd')
            param.bwd=round(G.N/20);
        end
        %%  the order of the polynomial approximation
        if ~isfield(param,'order')
            param.order=15;
        end
        param.lk= gsp_fast_estimate_lk(G, param.bwd, param);
        [~, param.jch_co] = jackson_cheby_poly_coefficients(0, param.lk, [0, G.lmax], param.order);
        param.T= gsp_cheby_op(G, param.jch_co, speye(G.N));
    end
    param.dis_node=zeros(G.N);
    diag_T=diag(param.T);
    for j=1:G.N
        param.dis_node(:,j)=sqrt(diag_T+ones(G.N,1)*diag_T(j)-2*param.T(:,j));
    end
end
dis_node=param.dis_node;


Rem_node=1:G.N;
J=0;
while ~isempty(Rem_node)
    L_rem=length(Rem_node);
    id_i=randi(L_rem);
    m=Rem_node(id_i);
    J=J+1;
    id=find(dis_node(m,Rem_node)<=param.epsilon/2);
    Partition{J}=Rem_node(id);
    Rem_node(id)=[];
end

if nargout==2
    varargout{1}=param;
end
end