ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3)),'_HMCprelim'])
    ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
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

ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3)),'_rand_HMCprelim'])
    ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
    thvec = (vertcat(chains{:}))';    
    rho_est_rand(z3,1:size(thvec,2)) = 2./(1+exp(-thvec(2,:)))-1;    
end

%%

h = colormap(hsv(size(ClonekeyN,2)));
h = h*0.6;

Clonekey = {'Rasal2','Pgk','Rbpj','Jarid2','Jam2','Nono','Tfdp2','Dstn','Spry4'};

subplot(3,3,1)
Parmcurr = CV';
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
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'a'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('CV \mu_i')


subplot(3,3,2)
Parmcurr = siscorr';
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

subplot(3,3,4)
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
ylim([-0.1 1])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'d'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('\rho')
%%
subplot(3,3,5)
Parmcurr = (rho_est_rand(:,1:200))';
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
ylim([-0.1 1])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'e'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('\rho')

%%
% for microenvironment

ClonekeyN = {'Rbpj_S1_STRSa1','Jam2_S1_STRSa1'};

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3)),'_HMC_micro_related'])
    ClonekeyN = {'Rbpj_S1_STRSa1','Jam2_S1_STRSa1'};
    thvec = (vertcat(chains{:}))';    
    rho_est_related(z3,1:size(thvec,2)) = 2./(1+exp(-thvec(2,:)))-1;  
end

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3)),'_HMC_micro_unrelated'])
    ClonekeyN = {'Rbpj_S1_STRSa1','Jam2_S1_STRSa1'};
    thvec = (vertcat(chains{:}))';    
    rho_est_unrelated(z3,1:size(thvec,2)) = 2./(1+exp(-thvec(2,:)))-1;    
end

h = colormap(hsv(size(ClonekeyN,2)));
h = h*0.6;
Clonekey = {'Rbpj sis','Jam2 sis','Rbpj non-sis','Jam2 non-sis'};

subplot(3,3,6)


Parmcurr = ([rho_est_related;rho_est_unrelated])';
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
ylim([-0.1 1])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'f'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('\rho')
%%

subplot(3,3,3)
scatter(mean(CV,2),mean(siscorr,2),'x')
xlabel('CV \mu_i')
ylabel('\lambda')

a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'c'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

title(['R = ',num2str(corr(mean(CV,2),mean(siscorr,2)),2)],'fontweight','normal')

%%

Clonekey9 = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
Clonekeyref9 = {'Rasal2','Pgk','Rbpj','Jarid2','Jam2','Nono','Tfdp2','Dstn','Spry4'};

 z9 = [2];
    subplot(3,3,7)
    currclonename = Clonekeyref9(z9);
    load([char(Clonekey9(z9)),'_HMCprelim'])
    concatenatedSamples = vertcat(chains{:});
    X = mean(concatenatedSamples,1);
    x = linspace(0,13.5,100);
    R0R0coeff = GetR0R0coeff(X,x,yTOT,gamma);%0.9;R0R0coeff = 0.9;
    [theo_corr] = GetTheoCorr_fullmodel(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_norho] = GetTheoCorr_norho(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_noextrinsic] = GetTheoCorr_noextrinsic(X,x,yTOT,gamma,R0R0coeff);

    normalised_y = zeros(100,size(num1,2));

    t = (0:5:(5*size(num1,1)-1))';

    for i = 1:size(num1,2)
         y = num1(:,i);
         ty = t;
         ty(y==0) = [];
         y(y==0) = [];
         xq = linspace(0,ty(end),100);
         vq2 = interp1(ty,y,xq,'spline'); 
         normalised_y(:,i) = vq2;
    end
   
    corrvec = zeros(1,size(normalised_y,1));
    firstcells = normalised_y(:,1:2:size(normalised_y,1));
    secondcells = normalised_y(:,2:2:size(normalised_y,1));

    for i = 1:size(normalised_y,1)
        corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
    end
   
    xnorm = 100*x/x(end);
    p1 = plot(xnorm,(corrvec),'color',[0 0.8 0]);
    hold on
    plot(xnorm,theo_corr)
    plot(xnorm,theo_corr_norho)
    plot(xnorm,theo_corr_noextrinsic)
    hold off
    xlabel('% cell cycle')
    ylabel('Corr (sisters)')
    ylim([0 1])

       legend('Data','Full model','Model - \rho=0','Model - \sigma=0','location','southwest') 

    title(currclonename,'fontweight','normal')
  a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'g'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


set(gca, 'FontName', 'Arial')

%%

Clonekey9 = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
Clonekeyref9 = {'Rasal2','Pgk','Rbpj','Jarid2','Jam2','Nono','Tfdp2','Dstn','Spry4'};


 z9 = [3];
    subplot(3,3,8)

    currclonename = Clonekeyref9(z9);
    load([char(Clonekey9(z9)),'_HMCprelim'])
    concatenatedSamples = vertcat(chains{:});
    X = mean(concatenatedSamples,1);
    x = linspace(0,13.5,100);
    R0R0coeff = GetR0R0coeff(X,x,yTOT,gamma);%0.9;
    [theo_corr] = GetTheoCorr_fullmodel(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_norho] = GetTheoCorr_norho(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_noextrinsic] = GetTheoCorr_noextrinsic(X,x,yTOT,gamma,R0R0coeff);

    normalised_y = zeros(100,size(num1,2));

    t = (0:5:(5*size(num1,1)-1))';

    for i = 1:size(num1,2)
         y = num1(:,i);
         ty = t;
         ty(y==0) = [];
         y(y==0) = [];
         xq = linspace(0,ty(end),100);
         vq2 = interp1(ty,y,xq,'spline'); 
         normalised_y(:,i) = vq2;
    end
    
    corrvec = zeros(1,size(normalised_y,1));
    firstcells = normalised_y(:,1:2:size(normalised_y,1));
    secondcells = normalised_y(:,2:2:size(normalised_y,1));

    for i = 1:size(normalised_y,1)
        corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
    end
    
    xnorm = 100*x/x(end);
    p1 = plot(xnorm,(corrvec),'color',[0 0.8 0]);
    hold on
    plot(xnorm,theo_corr)
    plot(xnorm,theo_corr_norho)
    plot(xnorm,theo_corr_noextrinsic)
    hold off
    xlabel('% cell cycle')
    ylabel('Corr (sisters)')
    ylim([0 1])
    title(currclonename,'fontweight','normal')
  a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'h'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
set(gca, 'FontName', 'Arial')

%%

Clonekey9 = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1_1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
Clonekeyref9 = {'Rasal2','Pgk','Rbpj','Jarid2','Jam2','Nono','Tfdp2','Dstn','Spry4'};


 z9 = [8];
    subplot(3,3,9)
    currclonename = Clonekeyref9(z9);
    load([char(Clonekey9(z9)),'_HMCprelim'])
    concatenatedSamples = vertcat(chains{:});
    X = mean(concatenatedSamples,1);
    x = linspace(0,13.5,100);
    R0R0coeff = GetR0R0coeff(X,x,yTOT,gamma);
    [theo_corr] = GetTheoCorr_fullmodel(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_norho] = GetTheoCorr_norho(X,x,yTOT,gamma,R0R0coeff);
    [theo_corr_noextrinsic] = GetTheoCorr_noextrinsic(X,x,yTOT,gamma,R0R0coeff);
    normalised_y = zeros(100,size(num1,2));

    t = (0:5:(5*size(num1,1)-1))';

    for i = 1:size(num1,2)
         y = num1(:,i);
         ty = t;
         ty(y==0) = [];
         y(y==0) = [];
         xq = linspace(0,ty(end),100);
         vq2 = interp1(ty,y,xq,'spline'); 
         normalised_y(:,i) = vq2;
    end
  
    corrvec = zeros(1,size(normalised_y,1));
    firstcells = normalised_y(:,1:2:size(normalised_y,1));
    secondcells = normalised_y(:,2:2:size(normalised_y,1));

    for i = 1:size(normalised_y,1)
        corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
    end
    
    xnorm = 100*x/x(end);
    p1 = plot(xnorm,(corrvec),'color',[0 0.8 0]);
    hold on
    plot(xnorm,theo_corr)
    plot(xnorm,theo_corr_norho)
    plot(xnorm,theo_corr_noextrinsic)
    hold off
    xlabel('% cell cycle')
    ylabel('Corr (sisters)')
    ylim([0 1])
    title(currclonename,'fontweight','normal')
  a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'i'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


set(gca, 'FontName', 'Arial')

