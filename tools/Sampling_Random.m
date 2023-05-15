function [iset_sampling,opt_weight]=Sampling_Random(G,K,Samplingsize)
%% The sampling method using random sampling distribution
%% Without replacement

  [~, cum_coh] = gsp_estimate_Samplingdistribution(G,K); %% �����������Ϊ1e-7,����ʽ�ƽ�ʹ��Jackson-Chebyshev����ʽ�������Ϊ30��
  % Sampling size
%   fprintf('* Ratio of the size of the sampling set to bandwidth K: %d\n',Multi_size)
  opt_weight = cum_coh/sum(cum_coh(:));
  iset_sampling = datasample(1:G.N, Samplingsize, ...
            'Replace', false, 'Weights', opt_weight);    %% sampling set (without replacement)      
  
end