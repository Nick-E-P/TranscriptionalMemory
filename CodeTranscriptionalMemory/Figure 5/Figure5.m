ClonekeyN = {'Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jam2_S1_STRSa1'};

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3)),'_MD_HMCprelim'])
    ClonekeyN = {'Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jam2_S1_STRSa1'};
    thvec = (vertcat(chains{:}))';
    mu_mu = thvec(3,:);
    sigma_mu = exp(thvec(6,:));
    variance = (exp(sigma_mu.^2)-1).*exp(2*mu_mu+sigma_mu.^2);
    meanlog = exp(mu_mu+sigma_mu.^2/2);
    CV(z3,1:size(thvec,2)) = (variance).^0.5./meanlog;
    siscorr(z3,1:size(thvec,2)) = tanh(thvec(9,:));
    
    tau_s(z3,1:size(thvec,2)) = exp(thvec(1,:));
    mu_sigma_s = exp(thvec(4,:));
    sigma_s_normed(z3,1:size(thvec,2)) = mu_sigma_s./exp(mu_mu);
    
    rho_est(z3,1:size(thvec,2)) = 2./(1+exp(-thvec(2,:)))-1;
    
end

%%
h = colormap(hsv(size(ClonekeyN,2)));
h = h*0.6;
Clonekey = {'Pgk','Rbpj','Jam2'};

%%
subplot(2,2,3)
Parmcurr = siscorr';%squeeze(Parms(:,4,:));
q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s = [q5;q25;q75;q95];
g = boxplot(Parmcurr,char(Clonekey),'colors',h,'OutlierSize',7);
set(g(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(g(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(g(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(g(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(g(7,:),'Visible','off');
xtickangle(45)
ylim([0 1])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'b'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('\lambda')

subplot(2,2,4)
Parmcurr = rho_est';
q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s = [q5;q25;q75;q95];
g = boxplot(Parmcurr,char(Clonekey),'colors',h,'OutlierSize',7);
set(g(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(g(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(g(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(g(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(g(7,:),'Visible','off');
xtickangle(45)
ylim([0 1])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'c'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('\rho')

%%

subplot(2,2,1:2)

load('Pgk_S1_STRSa1_MD_HMCprelim.mat')
Xref = mean(vertcat(chains{:}));

y = log(2)/gamma;

w = 15;

y1 = yTOT(:,w);
y2 = yTOT(:,w+1);
pairnum = round(w/2);
X = Xref([1:2,6*pairnum+4,6*pairnum+5,6*pairnum+6,6*pairnum+7,6*pairnum+8,6*pairnum+9]);

X([1,3:end]) = exp(X([1,3:end]));    

mu_1 = X(3);
mu_2 = X(6);

R0_1 = mu_1/2;
R0_2 = mu_2/2;

fin1 = find(y1,1,'last');
fin2 = find(y2,1,'last');

meanfn1 = mu_1+ (R0_1-mu_1)*exp(-y*x(1:fin1));
meanfn2 = mu_2+ (R0_2-mu_2)*exp(-y*x(1:fin2));

plot(x(1:fin1)/x(fin1),y1(1:fin1),'r')
hold on
plot(x(1:fin1)/x(fin1),meanfn1,'--k')
plot(1+x(1:fin2)/x(fin2),y2(1:fin2),'r')
plot(1+x(1:fin2)/x(fin2),meanfn2,'--k')

w = 25;

y1 = yTOT(:,w);
y2 = yTOT(:,w+1);
pairnum = round(w/2);
X = Xref([1:2,6*pairnum+4,6*pairnum+5,6*pairnum+6,6*pairnum+7,6*pairnum+8,6*pairnum+9]);

X([1,3:end]) = exp(X([1,3:end]));    

mu_1 = X(3);
mu_2 = X(6);

R0_1 = mu_1/2;
R0_2 = mu_2/2;

fin1 = find(y1,1,'last');
fin2 = find(y2,1,'last');

meanfn1 = mu_1+ (R0_1-mu_1)*exp(-y*x(1:fin1));
meanfn2 = mu_2+ (R0_2-mu_2)*exp(-y*x(1:fin2));

plot(x(1:fin1)/x(fin1),y1(1:fin1),'b')
plot(x(1:fin1)/x(fin1),meanfn1,'--k')
plot(1+x(1:fin2)/x(fin2),y2(1:fin2),'b')
plot(1+x(1:fin2)/x(fin2),meanfn2,'--k')
hold off

title('Mother cell    Daughter cell','fontweight','normal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'a'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('Bioluminescence [A.U.]')
xlabel('Time (generations)')

plot(X(10:3:end))
