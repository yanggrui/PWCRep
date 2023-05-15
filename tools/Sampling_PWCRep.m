function c=Sampling_PWCRep(x,Partition)
%% This function is the sampling function of PWC representation



%% Sampling
J=length(Partition);
c=zeros(J,1);
for j=1:J
    c(j)=mean(x(Partition{j}));
end
end