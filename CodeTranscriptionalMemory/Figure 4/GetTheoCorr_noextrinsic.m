function [theo_corr] = GetTheoCorr_fullmodel(X,x,yTOT,gamma,R0R0coeff)
%GETTHEOCORR_FULLMODEL Summary of this function goes here
%   Detailed explanation goes here


%%
Xref = X;
k = 1/exp(Xref(1));
y = log(2)/gamma;

K1_store = zeros(length(x),size(yTOT,2)/2);
K2_store = zeros(length(x),size(yTOT,2)/2);
K3_store = zeros(length(x),size(yTOT,2)/2);
E_X = zeros(length(x),size(yTOT,2)/2);
E_Y = zeros(length(x),size(yTOT,2)/2);


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


A = exp(- t*y - tau*y);
B = exp(- t*y);
C = exp(- t*k);
D = exp(-tau*y);
E = exp(-2*t*y);
F = exp(- k*t - k*tau);
G = exp(-k*tau);
H = exp(-2*k*t);

%%

for w = 1:2:size(yTOT,2)
    w
    pairnum = round(w/2);
    X = Xref([1:2,6*pairnum+4,6*pairnum+5,6*pairnum+6,6*pairnum+7,6*pairnum+8,6*pairnum+9]);
    
    X([1,3:end]) = exp(X([1,3:end]));
    
    rho = 2/(1+exp(-X(2)))-1;
    
    mu_1 = X(3);
    mu_2 = X(6);
    
    D1_1 = X(4);
    D1_2 = X(7);
    
    D2_1 = X(5);
    D2_2 = X(8);
    
    R0_1 = mu_1/2;
    R0_2 = mu_2/2;


    meanfn1 = mu_1+ (R0_1-mu_1)*exp(-y*x(1:end));
    meanfn2 = mu_2+ (R0_2-mu_2)*exp(-y*x(1:end));
    
    A1 = A(1:end,1:end);
    B1 = B(1:end,1:end);
    C1 = C(1:end,1:end);
    D1 = D(1:end,1:end);
    E1 = E(1:end,1:end);
    F1 = F(1:end,1:end);
    G1 = G(1:end,1:end);
    H1 = H(1:end,1:end);
    
    A2 = A(1:end,1:end);
    B2 = B(1:end,1:end);
    C2 = C(1:end,1:end);
    D2 = D(1:end,1:end);
    E2 = E(1:end,1:end);
    F2 = F(1:end,1:end);
    G2 = G(1:end,1:end);
    H2 = H(1:end,1:end);
 
    % goes to the right of K1 when constructing total SIGMA
    A3 = A(1:end,1:end);
    B3 = B(1:end,1:end);
    C3 = C(1:end,1:end);
    D3 = D(1:end,1:end);
    E3 = E(1:end,1:end);
    F3 = F(1:end,1:end);
    G3 = G(1:end,1:end);
    H3 = H(1:end,1:end);
    
    % define K matrices
    
%     K1 = A1*((B1*(D2_1*k + D2_2*k + D1_1*y + D1_2*y + D2_1*y + D2_2*y))/(8*(k + y)) + (y^2*(B1 - C1)*(D1_1 + D1_2))/(4*(k + y)*(k - y))) - D1*D2_1*(E1 - 1) - (D1_1*G1*y^2*(H1 - 1))/(k - y)^2 - (D1*D1_1*k*y*(E1 - 1))/(k - y)^2 + (y^2*(A1 - F1)*(D1_1 + D1_2)*(3*B1*k - 2*C1*k + B1*y - 2*C1*y))/(4*(k^2 - y^2)*(k - y)) + (2*D1*D1_1*k*y^2*(B1*C1 - 1))/((k + y)*(k - y)^2) + (2*D1_1*G1*k*y^2*(B1*C1 - 1))/((k + y)*(k - y)^2)
%     K2 = A2*((B2*(D2_1*k + D2_2*k + D1_1*y + D1_2*y + D2_1*y + D2_2*y))/(8*(k + y)) + (y^2*(B2 - C2)*(D1_1 + D1_2))/(4*(k + y)*(k - y))) - D2*D2_2*(E2 - 1) - (D1_2*G2*y^2*(H2 - 1))/(k - y)^2 - (D2*D1_2*k*y*(E2 - 1))/(k - y)^2 + (y^2*(A2 - F2)*(D1_1 + D1_2)*(3*B2*k - 2*C2*k + B2*y - 2*C2*y))/(4*(k^2 - y^2)*(k - y)) + (2*D2*D1_2*k*y^2*(B2*C2 - 1))/((k + y)*(k - y)^2) + (2*D1_2*G2*k*y^2*(B2*C2 - 1))/((k + y)*(k - y)^2)
%     K3 = A3*((B3*R0R0coeff*(D2_1*k + D2_2*k + D1_1*y + D1_2*y + D2_1*y + D2_2*y))/(8*(k + y)) + (R0R0coeff*y^2*(B3 - C3)*(D1_1 + D1_2))/(4*(k + y)*(k - y))) + (R0R0coeff*y^2*(A3 - F3)*(D1_1 + D1_2)*(3*B3*k - 2*C3*k + B3*y - 2*C3*y))/(4*(k^2 - y^2)*(k - y)) - (D1_1^(1/2)*D1_2^(1/2)*G3*rho*y^2*(H3 - 1))/(k - y)^2 - (D3*D1_1^(1/2)*D1_2^(1/2)*k*rho*y*(E3 - 1))/(k - y)^2 + (2*D3*D1_1^(1/2)*D1_2^(1/2)*k*rho*y^2*(B3*C3 - 1))/((k + y)*(k - y)^2) + (2*D1_1^(1/2)*D1_2^(1/2)*G3*k*rho*y^2*(B3*C3 - 1))/((k + y)*(k - y)^2)

    K1 = A1.*((B1.*(D2_1.*k + D2_2.*k + D1_1.*y + D1_2.*y + D2_1.*y + D2_2.*y))/(8.*(k + y)) + (y^2.*(B1 - C1).*(D1_1 + D1_2))/(4.*(k + y).*(k - y))) - D1.*D2_1.*(E1 - 1) - (D1_1.*G1.*y^2.*(H1 - 1))/(k - y)^2 - (D1.*D1_1.*k.*y.*(E1 - 1))/(k - y)^2 + (y^2.*(A1 - F1).*(D1_1 + D1_2).*(3.*B1.*k - 2.*C1.*k + B1.*y - 2.*C1.*y))/(4.*(k^2 - y^2).*(k - y)) + (2.*D1.*D1_1.*k.*y^2.*(B1.*C1 - 1))/((k + y).*(k - y)^2) + (2.*D1_1.*G1.*k.*y^2.*(B1.*C1 - 1))/((k + y).*(k - y)^2);
    K2 = A2.*((B2.*(D2_1.*k + D2_2.*k + D1_1.*y + D1_2.*y + D2_1.*y + D2_2.*y))/(8.*(k + y)) + (y^2.*(B2 - C2).*(D1_1 + D1_2))/(4.*(k + y).*(k - y))) - D2.*D2_2.*(E2 - 1) - (D1_2.*G2.*y^2.*(H2 - 1))/(k - y)^2 - (D2.*D1_2.*k.*y.*(E2 - 1))/(k - y)^2 + (y^2.*(A2 - F2).*(D1_1 + D1_2).*(3.*B2.*k - 2.*C2.*k + B2.*y - 2.*C2.*y))/(4.*(k^2 - y^2).*(k - y)) + (2.*D2.*D1_2.*k.*y^2.*(B2.*C2 - 1))/((k + y).*(k - y)^2) + (2.*D1_2.*G2.*k.*y^2.*(B2.*C2 - 1))/((k + y).*(k - y)^2);
    K3 = A3.*((B3.*R0R0coeff.*(D2_1.*k + D2_2.*k + D1_1.*y + D1_2.*y + D2_1.*y + D2_2.*y))/(8.*(k + y)) + (R0R0coeff.*y^2.*(B3 - C3).*(D1_1 + D1_2))/(4.*(k + y).*(k - y))) + (R0R0coeff.*y^2.*(A3 - F3).*(D1_1 + D1_2).*(3.*B3.*k - 2.*C3.*k + B3.*y - 2.*C3.*y))/(4.*(k^2 - y^2).*(k - y)) - (D1_1^(1/2).*D1_2^(1/2).*G3.*rho.*y^2.*(H3 - 1))/(k - y)^2 - (D3.*D1_1^(1/2).*D1_2^(1/2).*k.*rho.*y.*(E3 - 1))/(k - y)^2 + (2.*D3.*D1_1^(1/2).*D1_2^(1/2).*k.*rho.*y^2.*(B3.*C3 - 1))/((k + y).*(k - y)^2) + (2.*D1_1^(1/2).*D1_2^(1/2).*G3.*k.*rho.*y^2.*(B3.*C3 - 1))/((k + y).*(k - y)^2);


%     SIGMA = [[K1,K3];[K3',K2]];
%     SIGMA = SIGMA+diag((1^2).*ones(1,(end+end)));

    
    K1_store(1:end,round(w/2)) = diag(K1);
    K2_store(1:end,round(w/2)) = diag(K2);
    K3_store(1:end,round(w/2)) = diag(K3);
    
    E_X(1:end,round(w/2)) = meanfn1;
    E_Y(1:end,round(w/2)) = meanfn2;
end


for i = 1:length(x)
    E_cov = mean(K3_store(i,:));
    cov_E = 0;%cov(E_X(i,:),E_Y(i,:));
    E_var = mean([K1_store(i,:),K2_store(i,:)]);
    var_E = 0;%var([E_X(i,:),E_Y(i,:)]);
    theo_corr(i) = (E_cov + cov_E)/(E_var+var_E);
end
end

