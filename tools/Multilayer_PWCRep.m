function [Partition,varargout]=Multilayer_PWCRep(G,param)


%% Initialize presentation error level
if ~isfield(param,'error_tol')
    param.error_tol=0.1;
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









%% Initialize the initial epsilon
if ~isfield(param,'epsilon')
    param.epsilon=5*1e-2;
end
if ~isfield(param,'step')
    param.step=5*1e-2;
end


%% first step
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
Max_l=lambda_max_est(G,Partition,param);


dis_node_temp=dis_node;
dis_node=zeros(length(Partition));
for i=1:length(Partition)-1
    for j=i+1:length(Partition)
        dis_partition_ij=dis_node_temp(Partition{i},Partition{j});
        dis_node(i,j)=max(max(dis_partition_ij));
        dis_node(j,i)=dis_node(i,j);
    end
end




clear Rem_node id J dis_node_temp

if nargout>=4
    varargout{1}=param.epsilon;
    varargout{2}={Partition};
    varargout{3}=Max_l;
    varargout{4}=param;
end
while Max_l<=param.error_tol
    param.epsilon=param.epsilon+param.step;
    [Partition_temp,dis_node_temp]=Sub_Partition(Partition,dis_node,param.epsilon);
    Max_l=lambda_max_est(G,Partition_temp,param);
    if Max_l>param.error_tol
        param.epsilon=param.epsilon-param.step;
        break;
    else
        Partition=Partition_temp;
        dis_node=dis_node_temp;
        if nargout>=4
            varargout{1}=[varargout{1} param.epsilon];
            varargout{2}=[varargout{2} {Partition}];
            varargout{3}=[varargout{3} Max_l];
            varargout{4}=param;
        end
    end
end
end

function [New_Partition,New_dis_node]=Sub_Partition(Partition,dis_node,epsilon)
M=size(dis_node,1);
Rem_node=1:M;


J=0;
while ~isempty(Rem_node)
    m=Rem_node(1);
    J=J+1;
    id=find(dis_node(m,Rem_node)<=epsilon/2);
    SubPartition{J}=Rem_node(id);
    Rem_node(id)=[];
end

J=length(SubPartition);
for j=1:J
    New_Partition{j}=cell2mat(Partition(SubPartition{j}));
end


New_dis_node=zeros(J);
for i=1:J-1
    for j=i+1:J
        dis_partition_ij=dis_node(SubPartition{i},SubPartition{j});
        New_dis_node(i,j)=max(max(dis_partition_ij));
        New_dis_node(j,i)=New_dis_node(i,j);
    end
end


end