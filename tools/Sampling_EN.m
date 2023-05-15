function selected_iset=Sampling_EN(UK,M)
%% This function is to performance the sampling method of EN

selected_iset=[];
while M>0
iset_m=OneSelection_EN(UK,M,selected_iset);
selected_iset=[selected_iset iset_m];
M=M-length(iset_m);
end

end

function iset=OneSelection_EN(UK,M,ex_iset)
[N,K]=size(UK);
d=min(M,K);
 
v=zeros(N,d);
 
iset=zeros(1,d);
tau_output=zeros(1,d);

Rem_node=1:N;
Rem_node(ex_iset)=[];
L_rem=length(Rem_node);
tau=zeros(N,1);
for i=1:L_rem
j=Rem_node(i);
tau(j)=UK(j,:)*UK(j,:)';
end

for m=1:d
    [~,im]=max(tau);
    tau_output(m)=sqrt(tau(im)); 
    iset(m)=im;
    vm1=UK*(UK(im,:)'/tau_output(m));
    if m>1 && m<d
        vm1=vm1-v(:,1:m-1)*(v(im,1:m-1)'/tau_output(m));
    end
    if m==d
        break;
    end
    v(:,m)=vm1;     
    tau(Rem_node)=tau(Rem_node)-abs(vm1(Rem_node)).^2;
end
end

