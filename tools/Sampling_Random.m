function [iset_sampling,opt_weight]=Sampling_Random(G,K,Samplingsize)
%% The sampling method using random sampling distribution
%% Without replacement

  [~, cum_coh] = gsp_estimate_Samplingdistribution(G,K); %% 迭代容忍误差为1e-7,多项式逼近使用Jackson-Chebyshev多项式，最高项为30；
  % Sampling size
%   fprintf('* Ratio of the size of the sampling set to bandwidth K: %d\n',Multi_size)
  opt_weight = cum_coh/sum(cum_coh(:));
  iset_sampling = datasample(1:G.N, Samplingsize, ...
            'Replace', false, 'Weights', opt_weight);    %% sampling set (without replacement)      
  
end