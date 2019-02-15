function [logprior, gradlogprior] = PriorLevel(X,yTOT)
%PRIORLEVEL Caluculates the prior i.e. the "population level" as described
%in the text
logpriorCURR = zeros(size(yTOT,2),1);
gradlogprior = zeros(length(X),1);

%%
for w = 1:2:size(yTOT,2)
    %%
    
    pairnum = round(w/2);
    Xcurr = X([3:9,6*pairnum+4,6*pairnum+5,6*pairnum+6,6*pairnum+7,6*pairnum+8,6*pairnum+9]);
    
    mu_mu = Xcurr(1);
    mu_D1 = Xcurr(2);
    mu_D2 = Xcurr(3);
    sigma_mu = exp(Xcurr(4));
    sigma_D1 = exp(Xcurr(5));
    sigma_D2 = exp(Xcurr(6));
    y_4 = Xcurr(7);
    mu_1 = Xcurr(8);
    D1_1 = log(sqrt(exp(Xcurr(9))));
    D2_1 = log(sqrt(exp(Xcurr(10))));
    mu_2 = Xcurr(11);
    D1_2 = log(sqrt(exp(Xcurr(12))));
    D2_2 = log(sqrt(exp(Xcurr(13))));

    ydet = [mu_1;D1_1;D2_1;mu_2;D1_2;D2_2]-[mu_mu;mu_D1;mu_D2;mu_mu;mu_D1;mu_D2];
    
    SIGMA = ... 
[           sigma_mu^2,          0,          0, sigma_mu^2*tanh(y_4),          0,          0;
                    0, sigma_D1^2,          0,                    0,          0,          0;
                    0,          0, sigma_D2^2,                    0,          0,          0;
 sigma_mu^2*tanh(y_4),          0,          0,           sigma_mu^2,          0,          0;
                    0,          0,          0,                    0, sigma_D1^2,          0;
                    0,          0,          0,                    0,          0, sigma_D2^2];

    
    n = length(ydet);
    L = chol(SIGMA);
    alpha = L\(L'\ydet);

    logpriorCURR(round(w/2)) = -0.5*ydet'*alpha - sum(diag(log(L))) - 0.5*n*log(2*pi);
    
    dSIGMA_1_dsigma_mu = sigma_mu*...
[           2*sigma_mu, 0, 0, 2*sigma_mu*tanh(y_4), 0, 0;
                    0, 0, 0,                    0, 0, 0;
                    0, 0, 0,                    0, 0, 0;
 2*sigma_mu*tanh(y_4), 0, 0,           2*sigma_mu, 0, 0;
                    0, 0, 0,                    0, 0, 0;
                    0, 0, 0,                    0, 0, 0];


    dSIGMA_1_dsigma_D1 = sigma_D1*... 
    [ 0,          0, 0, 0,          0, 0;
     0, 2*sigma_D1, 0, 0,          0, 0;
     0,          0, 0, 0,          0, 0;
     0,          0, 0, 0,          0, 0;
     0,          0, 0, 0, 2*sigma_D1, 0;
     0,          0, 0, 0,          0, 0];


    dSIGMA_1_dsigma_D2 = sigma_D2*...
    [ 0, 0,          0, 0, 0,          0;
     0, 0,          0, 0, 0,          0;
     0, 0, 2*sigma_D2, 0, 0,          0;
     0, 0,          0, 0, 0,          0;
     0, 0,          0, 0, 0,          0;
     0, 0,          0, 0, 0, 2*sigma_D2];


    dSIGMA_1_dsigma_y_4 = ... 
    [                             0, 0, 0, -sigma_mu^2*(tanh(y_4)^2 - 1), 0, 0;
                                 0, 0, 0,                             0, 0, 0;
                                 0, 0, 0,                             0, 0, 0;
     -sigma_mu^2*(tanh(y_4)^2 - 1), 0, 0,                             0, 0, 0;
                                 0, 0, 0,                             0, 0, 0;
                                 0, 0, 0,                             0, 0, 0];

    dmean_du1 = [1;0;0;1;0;0];
    dmean_du2 = [0;1;0;0;1;0];
    dmean_du3 = [0;0;1;0;0;1];

    dmean_dy1 = [-1;0;0;0;0;0];
    dmean_dy2 = 0.5*[0;-1;0;0;0;0]; % 0.5 as original param varaince not std
    dmean_dy3 = 0.5*[0;0;-1;0;0;0];
    dmean_dy4 = [0;0;0;-1;0;0];
    dmean_dy5 = 0.5*[0;0;0;0;-1;0];
    dmean_dy6 = 0.5*[0;0;0;0;0;-1];

    % for parameters affecting variance
    dL_dsigma_mu = 0.5*trace(alpha*alpha'*dSIGMA_1_dsigma_mu - L\(L'\dSIGMA_1_dsigma_mu));
    dL_dsigma_D1 = 0.5*trace(alpha*alpha'*dSIGMA_1_dsigma_D1 - L\(L'\dSIGMA_1_dsigma_D1));
    dL_dsigma_D2 = 0.5*trace(alpha*alpha'*dSIGMA_1_dsigma_D2 - L\(L'\dSIGMA_1_dsigma_D2));
    dL_dsigma_y_4 = 0.5*trace(alpha*alpha'*dSIGMA_1_dsigma_y_4 - L\(L'\dSIGMA_1_dsigma_y_4));

    % for parameters affecting mean
    dL_dmu_mu = dmean_du1'*alpha;
    dL_dmu_D1 = dmean_du2'*alpha;
    dL_dmu_D2 = dmean_du3'*alpha;
    
    dL_dmu_1 = dmean_dy1'*alpha;
    dL_dD1_1 = dmean_dy2'*alpha;
    dL_dD2_1 = dmean_dy3'*alpha;
    
    dL_dmu_2 = dmean_dy4'*alpha;
    dL_dD1_2 = dmean_dy5'*alpha;
    dL_dD2_2 = dmean_dy6'*alpha;
    
    gradlogpriorTEMP = zeros(length(X),1);    
    
    gradlogpriorTEMP(3:9) = [dL_dmu_mu;dL_dmu_D1;dL_dmu_D2;dL_dsigma_mu;dL_dsigma_D1;dL_dsigma_D2;dL_dsigma_y_4];
    gradlogpriorTEMP(6*pairnum+4:6*pairnum+9) = [dL_dmu_1;dL_dD1_1;dL_dD2_1;dL_dmu_2;dL_dD1_2;dL_dD2_2];
    
    gradlogprior = gradlogprior+gradlogpriorTEMP;
    
end
logprior = sum(logpriorCURR);
end

