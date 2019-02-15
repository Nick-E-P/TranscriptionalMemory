function [HyperpriorTOT, gradlogHyperprior] = HyperPriorLevel(X)
%HYPERPRIORLEVEL Calculates the priors on \tau, \rho and \lambda
logHyperPrior = zeros(length(X),1);
gradlogHyperprior = zeros(length(X),1);

Xvar = X(1);
Mu = 1; 
Sigma = 1;
Z = (Xvar - Mu)./Sigma;
logHyperPrior(1) = sum(-log(Sigma) - .5*log(2*pi) - .5*(Z.^2));
gradlogHyperprior(1) = -Z./Sigma;

Xvar = X(2);
logHyperPrior(2) = log(exp(Xvar)/(exp(Xvar) + 1)^2);
gradlogHyperprior(2) = -exp(-Xvar)*(exp(Xvar) + 1)^2*((2*exp(2*Xvar))/(exp(Xvar) + 1)^3 - exp(Xvar)/(exp(Xvar) + 1)^2);

y_4 = X(9);
n = 5;
logHyperPrior(9) = log((1 - tanh(y_4)^2)^(n - 1));
  
gradlogHyperprior(9) = -tanh(y_4)*(2*n - 2);

HyperpriorTOT = sum(logHyperPrior);
end

