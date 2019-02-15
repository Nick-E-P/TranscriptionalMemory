function [logpdf,gradlogpdf] = logPosteriorReworked(X,x,yTOT,gamma,R0R0coeff)
%LOGPOSTERIORREWORKED Caluculates log posterior of the model
%%
[HyperpriorTOT, gradlogHyperprior] = HyperPriorLevel(X);

[logprior, gradlogprior] = PriorLevel(X,yTOT);

[loglik, gradloglik] = LikelihoodLevel3(X,x,yTOT,gamma,R0R0coeff);

logpdf = HyperpriorTOT + logprior + loglik;

gradlogpdf = gradlogHyperprior + gradlogprior+gradloglik;

end

