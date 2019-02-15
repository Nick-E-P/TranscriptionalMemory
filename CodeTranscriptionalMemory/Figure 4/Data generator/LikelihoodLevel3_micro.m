function [loglik, gradloglik] = LikelihoodLevel3(Xref,x,yTOT,gamma,R0R0coeff)
%LIKELIHOODLEVEL Summary of this function goes here
%   first get general matrices for A,..,H
%%

Likcurr = zeros(size(yTOT,2)/2,1);
gradloglik = zeros(length(Xref),1);

k = 1/exp(Xref(1));
y = log(2)/gamma;

rho = 2/(1+exp(-Xref(2)))-1;

if abs(k-y)<0.01
    k = k+0.02;
end

t = zeros(length(x),length(x));
tau = zeros(length(x),length(x));

for i = 1:length(x)
    for j = 1:length(x)
          t(i,j) = min(x(j),x(i));
          tau(i,j) = abs(x(i)-x(j));
    end
end

A1 = exp(- 2*t*y - tau*y);
B1 = exp(- 2*k*t - k*tau);
C1 =exp(- k*t - k*tau - t*y);
D1 =exp(- k*t - t*y - tau*y);
E1 =exp(-tau*y);
F1 =exp(-k*tau);

dK_dTAU_dD1_1 = k^2.*((tau.*y^2)/(k^2 - 2.*k.*y + y^2) - (2.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (y^2.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2 + (2.*k.*tau.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (2.*k.*y^2.*(- 3.*k^2 + 2.*k.*y + y^2))/(- k^3 + k^2.*y + k.*y^2 - y^3)^2).*F1 + k^2.*((4.*y)/(4.*k + 4.*y)^2 + y/(k^2 - 2.*k.*y + y^2) + (8.*k.*y^2)/(2.*k^2 - 2.*y^2)^2 + (y^2.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2 - (k.*y.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2).*A1 + (-k^2.*((y^2.*(t + tau))/(k^2 - 2.*k.*y + y^2) - (2.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (y^2.*(t + tau))/(2.*k^2 - 2.*y^2) + (4.*k.*y^2)/(2.*k^2 - 2.*y^2)^2 + (y^2.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2 + (2.*k.*y^2.*(t + tau))/(- k^3 + k^2.*y + k.*y^2 - y^3) + (2.*k.*y^2.*(- 3.*k^2 + 2.*k.*y + y^2))/(- k^3 + k^2.*y + k.*y^2 - y^3)^2)).*C1 + (-k^2.*((t.*y^2)/(k^2 - 2.*k.*y + y^2) - (2.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (4.*k.*y^2)/(2.*k^2 - 2.*y^2)^2 + (t.*y^2)/(2.*k^2 - 2.*y^2) + (y^2.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2 + (2.*k.*t.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (2.*k.*y^2.*(- 3.*k^2 + 2.*k.*y + y^2))/(- k^3 + k^2.*y + k.*y^2 - y^3)^2)).*D1 + (-k^2.*(y/(k^2 - 2.*k.*y + y^2) + (2.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) - (2.*k.*y^2.*(- 3.*k^2 + 2.*k.*y + y^2))/(- k^3 + k^2.*y + k.*y^2 - y^3)^2 - (k.*y.*(2.*k - 2.*y))/(k^2 - 2.*k.*y + y^2)^2)).*E1;
dK_dTAU_dD2_1 = k^2*((4*k)/(4*k + 4*y)^2 + (4*y)/(4*k + 4*y)^2 - 1/(4*k + 4*y))*A1;
dK_dD1_1 =(y^2/(k^2 - 2*k*y + y^2) + (2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*F1 + (y/(4*k + 4*y) + (2*y^2)/(2*k^2 - 2*y^2) + y^2/(k^2 - 2*k*y + y^2) - (k*y)/(k^2 - 2*k*y + y^2))*A1 + ((2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3) + (k*y)/(k^2 - 2*k*y + y^2))*E1 + (- y^2/(2*k^2 - 2*y^2) - y^2/(k^2 - 2*k*y + y^2) - (2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*C1 + (- y^2/(2*k^2 - 2*y^2) - y^2/(k^2 - 2*k*y + y^2) - (2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*D1;
dK_dD2_1 =E1 + (k/(4*k + 4*y) + y/(4*k + 4*y) - 1)*A1;
dK3_dTAU_dD1_12 =((2.*k^2.*rho.*y^2)/(k^3 - 3.*k^2.*y + 3.*k.*y^2 - y^3) - (4.*k^3.*rho.*y^2)/(k^4 - 2.*k^3.*y + 2.*k.*y^3 - y^4) - (2.*k^3.*rho.*y^2)/(k^4 - 2.*k^2.*y^2 + y^4) - (2.*k^2.*rho.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (2.*k^3.*rho.*tau.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) + (k^2.*rho.*tau.*y^2)/(k^2 - 2.*k.*y + y^2)).*F1 + ((4.*k^3.*rho.*y^2)/(k^4 - 2.*k^3.*y + 2.*k.*y^3 - y^4) + (2.*k^3.*rho.*y^2)/(k^4 - 2.*k^2.*y^2 + y^4) + (2.*k^2.*rho.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) - (2.*k^3.*rho.*t.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3)).*D1 + ((4.*k^3.*rho.*y^2)/(k^4 - 2.*k^3.*y + 2.*k.*y^3 - y^4) + (2.*k^3.*rho.*y^2)/(k^4 - 2.*k^2.*y^2 + y^4) + (2.*k^2.*rho.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) - (2.*k^3.*rho.*t.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3) - (2.*k^3.*rho.*tau.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3)).*C1 + ((2.*k^3.*rho.*y)/(k^3 - 3.*k^2.*y + 3.*k.*y^2 - y^3) - (4.*k^3.*rho.*y^2)/(k^4 - 2.*k^3.*y + 2.*k.*y^3 - y^4) - (k^2.*rho.*y)/(k^2 - 2.*k.*y + y^2) - (2.*k^3.*rho.*y^2)/(k^4 - 2.*k^2.*y^2 + y^4) - (2.*k^2.*rho.*y^2)/(- k^3 + k^2.*y + k.*y^2 - y^3)).*E1 + (- (2.*k^2.*rho.*y^2)/(k^3 - 3.*k^2.*y + 3.*k.*y^2 - y^3) - (2.*k^2.*rho.*t.*y^2)/(k^2 - 2.*k.*y + y^2) - (k^2.*rho.*tau.*y^2)/(k^2 - 2.*k.*y + y^2)).*B1 + ((k^2.*rho.*y)/(k^2 - 2.*k.*y + y^2) - (2.*k^3.*rho.*y)/(k^3 - 3.*k^2.*y + 3.*k.*y^2 - y^3)).*A1;
dK3_drho_dD1_12 =(y^2/(k^2 - 2*k*y + y^2) + (2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*F1 + ((2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3) + (k*y)/(k^2 - 2*k*y + y^2))*E1 + (-y^2/(k^2 - 2*k*y + y^2))*B1 + (-(k*y)/(k^2 - 2*k*y + y^2))*A1 + (-(2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*C1 + (-(2*k*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*D1;
dK3_dD1_12 = ((rho*y^2)/(k^2 - 2*k*y + y^2) + (2*k*rho*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*F1 + ((k*rho*y)/(k^2 - 2*k*y + y^2) + (2*k*rho*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*E1 + (-(2*k*rho*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*C1 + (-(2*k*rho*y^2)/(- k^3 + k^2*y + k*y^2 - y^3))*D1 + (-(rho*y^2)/(k^2 - 2*k*y + y^2))*B1 + (-(k*rho*y)/(k^2 - 2*k*y + y^2))*A1;

%%

for w = 1:2:size(yTOT,2)
    %%
    
    pairnum = round(w/2);
    X = Xref([1:2,6*pairnum+4,6*pairnum+5,6*pairnum+6,6*pairnum+7,6*pairnum+8,6*pairnum+9]);
    
    X([1,3:end]) = exp(X([1,3:end]));    
    
    mu_1 = X(3);
    mu_2 = X(6);
    
    D1_1 = X(4);
    D1_2 = X(7);
    D1_12 = sqrt(D1_1)*sqrt(D1_2);
    
    dD1_12_dD1_1 = D1_2^(1/2)/(2*D1_1^(1/2));
    dD1_12_dD1_2 = D1_1^(1/2)/(2*D1_2^(1/2));
    
    D2_1 = X(5);
    D2_2 = X(8);
    
    R0_1 = mu_1/2;
    R0_2 = mu_2/2;

    y1 = yTOT(:,w);
    y2 = yTOT(:,w+1);

    fin1 = find(y1,1,'last');
    fin2 = find(y2,1,'last');
    
    y1 = y1(1:fin1);
    y2 = y2(1:fin2);
    
    start1 = find(y1,1,'first');
    start2 = find(y2,1,'first');
    
    y1 = y1(start2:fin1); % y2 must always start second!!
    y2 = y2(start2:fin2);

%     meanfn1 = mu_1+ (R0_1-mu_1)*exp(-y*x(1:fin1));
%     meanfn2 = mu_2+ (R0_2-mu_2)*exp(-y*x(1:fin2));
    
    meanfn1 = mu_1+ (R0_1-mu_1)*exp(-y*x(start2:fin1));
    meanfn2 = mu_2+ (R0_2-mu_2)*exp(-y*(x(start2:fin2)-x(start2)));
    
%     dmean_dmu1x = mu_1-0.5*mu_1*exp(-y*x(1:fin1));
%     dmean_dmu2x = mu_2-0.5*mu_2*exp(-y*x(1:fin2));
    
    dmean_dmu1x = mu_1-0.5*mu_1*exp(-y*x(start2:fin1));
    dmean_dmu2x = mu_2-0.5*mu_2*exp(-y*(x(start2:fin2)-x(start2)));
    
    dmean_dmu1 = [dmean_dmu1x;zeros(size(dmean_dmu2x))];
    dmean_dmu2 = [zeros(size(dmean_dmu1x));dmean_dmu2x];
    
    y1 = y1-meanfn1;
    y2 = y2-meanfn2;

    ydetrend = [y1;y2];
    
    
%     K1 = D1_1*dK_dD1_1(1:fin1,1:fin1)+D2_1*dK_dD2_1(1:fin1,1:fin1);
%     K2 = D1_2*dK_dD1_1(1:fin2,1:fin2)+D2_2*dK_dD2_1(1:fin2,1:fin2);
%     K3 = D1_12*dK3_dD1_12(1:fin1,1:fin2);
    
    K1 = D1_1*dK_dD1_1(start2:fin1,start2:fin1)+D2_1*dK_dD2_1(start2:fin1,start2:fin1);
    K2 = D1_2*dK_dD1_1(1:length(y2),1:length(y2))+D2_2*dK_dD2_1(1:length(y2),1:length(y2));
    K3 = D1_12*dK3_dD1_12(start2:fin1,1:length(y2));
    
    
    SIGMA = [[K1,K3];[K3',K2]];
    SIGMA = SIGMA+diag((1^2).*ones(1,(length(y1)+length(y2))));
    
   
    n = length(ydetrend);
    
    L = chol(SIGMA);
    
    alpha = L\(L'\ydetrend);
    Likcurr(round(w/2)) = -0.5*ydetrend'*alpha - sum(diag(log(L))) - 0.5*n*log(2*pi);
 
    
%     dK1_dTAU = D1_1*dK_dTAU_dD1_1(1:fin1,1:fin1)+D2_1*dK_dTAU_dD2_1(1:fin1,1:fin1);
%     dK1_drho = zeros(fin1,fin1);
%     dK1_dD1_1 = dK_dD1_1(1:fin1,1:fin1);
%     dK1_dD1_2 = zeros(fin1,fin1);
%     dK1_dD2_1 = dK_dD2_1(1:fin1,1:fin1);
%     dK1_dD2_2 = zeros(fin1,fin1);
%     dK2_dTAU = D1_2*dK_dTAU_dD1_1(1:fin2,1:fin2)+D2_2*dK_dTAU_dD2_1(1:fin2,1:fin2);
%     dK2_drho = zeros(fin2,fin2);
%     dK2_dD1_1 = zeros(fin2,fin2);
%     dK2_dD1_2 = dK_dD1_1(1:fin2,1:fin2);
%     dK2_dD2_1 = zeros(fin2,fin2);
%     dK2_dD2_2 = dK_dD2_1(1:fin2,1:fin2);
%     dK3_dTAU = D1_12*dK3_dTAU_dD1_12(1:fin1,1:fin2);
%     dK3_drho = D1_12*dK3_drho_dD1_12(1:fin1,1:fin2);
%     dK3_dD1_1 = dD1_12_dD1_1*dK3_dD1_12(1:fin1,1:fin2);
%     dK3_dD1_2 = dD1_12_dD1_2*dK3_dD1_12(1:fin1,1:fin2);
%     dK3_dD2_1 = zeros(fin1,fin2);
%     dK3_dD2_2 = zeros(fin1,fin2);

    dK1_dTAU = D1_1*dK_dTAU_dD1_1(start2:fin1,start2:fin1)+D2_1*dK_dTAU_dD2_1(start2:fin1,start2:fin1);
    dK1_drho = zeros(length(y1),length(y1));
    dK1_dD1_1 = dK_dD1_1(start2:fin1,start2:fin1);
    dK1_dD1_2 = zeros(length(y1),length(y1));
    dK1_dD2_1 = dK_dD2_1(start2:fin1,start2:fin1);
    dK1_dD2_2 = zeros(length(y1),length(y1));
    dK2_dTAU = D1_2*dK_dTAU_dD1_1(1:length(y2),1:length(y2))+D2_2*dK_dTAU_dD2_1(1:length(y2),1:length(y2));
    dK2_drho = zeros(length(y2),length(y2));
    dK2_dD1_1 = zeros(length(y2),length(y2));
    dK2_dD1_2 = dK_dD1_1(1:length(y2),1:length(y2));
    dK2_dD2_1 = zeros(length(y2),length(y2));
    dK2_dD2_2 = dK_dD2_1(1:length(y2),1:length(y2));
    dK3_dTAU = D1_12*dK3_dTAU_dD1_12(start2:fin1,1:length(y2));
    dK3_drho = D1_12*dK3_drho_dD1_12(start2:fin1,1:length(y2));
    dK3_dD1_1 = dD1_12_dD1_1*dK3_dD1_12(start2:fin1,1:length(y2));
    dK3_dD1_2 = dD1_12_dD1_2*dK3_dD1_12(start2:fin1,1:length(y2));
    dK3_dD2_1 = zeros(length(y1),length(y2));
    dK3_dD2_2 = zeros(length(y1),length(y2));
    
    
    
    K_inv = L\(L'\eye(length(y1)+length(y2),length(y1)+length(y2)));
    pre_term = alpha*alpha'- K_inv;
%     sum(A.*B',2)
    
    % calculate gradients - tau (parameter tau)
    SIGMAdifftau = X(1)*[[dK1_dTAU,dK3_dTAU];[dK3_dTAU',dK2_dTAU]];
%     tic
%     dL_dtau = 0.5*trace(alpha*alpha'*SIGMAdifftau - L\(L'\SIGMAdifftau));
%     toc
%     tic
%     dL_dtau = 0.5*trace((alpha*alpha'- K_inv)*SIGMAdifftau);
%     toc
%     tic
    dL_dtau = 0.5*trace(diag(sum(pre_term.*SIGMAdifftau',2)));
%     toc
    % gradient for rho
    
    SIGMAdiffrho = (2*exp(-X(2)))/(exp(-X(2)) + 1)^2*[[dK1_drho,dK3_drho];[dK3_drho',dK2_drho]];
%     dL_drho = 0.5*trace(alpha*alpha'*SIGMAdiffrho - L\(L'\SIGMAdiffrho));
    dL_drho = 0.5*trace(diag(sum(pre_term.*SIGMAdiffrho',2)));

    % gradient for mu1
    
    dL_dmu1 = (dmean_dmu1)'*alpha;
    
    % gradient for D1_1
    SIGMAdiffD1_1 = D1_1*[[dK1_dD1_1,dK3_dD1_1];[dK3_dD1_1',dK2_dD1_1]];
    
%     dL_dD1_1 = 0.5*trace(alpha*alpha'*SIGMAdiffD1_1 - L\(L'\SIGMAdiffD1_1));
    dL_dD1_1 = 0.5*trace(diag(sum(pre_term.*SIGMAdiffD1_1',2)));

  
    % gradient for D2_1
    
    SIGMAdiffD2_1 = D2_1*[[dK1_dD2_1,dK3_dD2_1];[dK3_dD2_1',dK2_dD2_1]];
%     dL_dD2_1 = 0.5*trace(alpha*alpha'*SIGMAdiffD2_1 - L\(L'\SIGMAdiffD2_1)); 
    dL_dD2_1 = 0.5*trace(diag(sum(pre_term.*SIGMAdiffD2_1',2))); 

    
    % gradient for mu2
    
    dL_dmu2 = (dmean_dmu2)'*alpha;
    
    % gradient for D1_2
    SIGMAdiffD1_2 = D1_2*[[dK1_dD1_2,dK3_dD1_2];[dK3_dD1_2',dK2_dD1_2]];
%     dL_dD1_2 = 0.5*trace(alpha*alpha'*SIGMAdiffD1_2 - L\(L'\SIGMAdiffD1_2));
    dL_dD1_2 = 0.5*trace(diag(sum(pre_term.*SIGMAdiffD1_2',2))); 

    % gradient for D2_2
    
    SIGMAdiffD2_2 = D2_2*[[dK1_dD2_2,dK3_dD2_2];[dK3_dD2_2',dK2_dD2_2]];
%     dL_dD2_2 = 0.5*trace(alpha*alpha'*SIGMAdiffD2_2 - L\(L'\SIGMAdiffD2_2)); 
    dL_dD2_2 = 0.5*trace(diag(sum(pre_term.*SIGMAdiffD2_2',2))); 

    
    gradloglikTEMP = zeros(length(Xref),1);
    gradloglikTEMP(1:2) = [dL_dtau;dL_drho];
    gradloglikTEMP(6*pairnum+4:6*pairnum+9) = [dL_dmu1;dL_dD1_1;dL_dD2_1;dL_dmu2;dL_dD1_2;dL_dD2_2];
    
    gradloglik = gradloglik+gradloglikTEMP;
    
%     gradloglik = [dL_dtau;dL_drhot;dL_dmu;dL_dD1;dL_dD2];
end


loglik = sum(Likcurr);

end

