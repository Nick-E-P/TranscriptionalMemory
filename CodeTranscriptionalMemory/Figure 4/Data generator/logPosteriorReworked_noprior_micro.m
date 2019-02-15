function [logpdf,gradlogpdf] = logPosteriorReworked_noprior_micro(X,x,yTOT,gamma,R0R0coeff)
%LOGPOSTERIORREWORKED Summary of this function goes here
%   Detailed explanation goes here
[HyperpriorTOT, gradlogHyperprior] = HyperPriorLevel(X);

% [logprior, gradlogprior] = PriorLevel(X,yTOT);

[loglik, gradloglik] = LikelihoodLevel3_micro(X,x,yTOT,gamma,R0R0coeff);

logpdf = HyperpriorTOT +  loglik;

gradlogpdf = gradlogHyperprior + gradloglik;

end

