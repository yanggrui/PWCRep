function y=Reconstruction_PWCRep(G,c,Partition)
%% This function is the reconstruction function of PWC representation

J=length(Partition);
y=zeros(G.N,1);
for j=1:J
    y(Partition{j})=c(j);
end

end