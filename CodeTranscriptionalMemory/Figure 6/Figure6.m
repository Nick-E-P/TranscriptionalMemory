[num]=xlsread('Families preliminary.xlsx',2); % RBPJ 

FamIdnt1 = num(6,:);
Level1 = num(5,:);

[num]=xlsread('Families preliminary.xlsx',3); % RBPJ 

FamIdnt2 = num(6,:);
Level2 = num(5,:);

FamIdnt = [FamIdnt1,FamIdnt2];
Level = [Level1,Level2];

% need to re-order
famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistJam2(i) = stdlist(i)/meanlist(i);
end

[I,A] = sort(meanlist)
[Ny,Nya] = sort(A)

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   Newlist(famidntstemp) = Nya(i)*ones(size(cellscurr))
end

subplot(3,3,2)
scatter(Newlist,Level/mean(Level),'bx')
hold on
scatter(1:max(FamIdnt),I/mean(Level),'k')
hold off
xlabel('Family number')
ylabel('Relative expression')
xlim([0 max(FamIdnt)])
ylim([0 3])
title('Jam2','fontweight','normal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'b'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
%%

clear
[num]=xlsread('Dstn families modded.xlsx',3); % Dstn 

FamIdnt = num(9,:);
Level = mean(num(5:7,:)); %using final intensity values

famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistDstn(i) = stdlist(i)/meanlist(i); 
end

[I,A] = sort(meanlist)
[Ny,Nya] = sort(A)

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   Newlist(famidntstemp) = Nya(i)*ones(size(cellscurr))
end

subplot(3,3,1)
scatter(Newlist,Level/mean(Level),'bx')
hold on
scatter(1:max(FamIdnt),I/mean(Level),'k')
hold off
xlabel('Family number')
ylabel('Relative expression')
legend('Individual cells','Family mean','location','northwest')
xlim([0 max(FamIdnt)])
title('Dstn','fontweight','normal')
ylim([0 3])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'a'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

clear
[num]=xlsread('PGK families.xlsx',2); % Pgk 

FamIdnt = num(5,:);
Level = mean(num(2:4,:)); %using final intensity values

famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistPGK(i) = stdlist(i)/meanlist(i);
end

[I,A] = sort(meanlist)
[Ny,Nya] = sort(A)

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   Newlist(famidntstemp) = Nya(i)*ones(size(cellscurr))
end

subplot(3,3,3)
scatter(Newlist,Level/mean(Level),'bx')
hold on
scatter(1:max(FamIdnt),I/mean(Level),'k')
hold off
xlabel('Family number')
ylabel('Relative expression')
xlim([0 max(FamIdnt)])
ylim([0 3])
title('Pgk','fontweight','normal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'c'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

%%
% create boxplot of mean values

[num]=xlsread('Dstn families modded.xlsx',3); % Dstn 

FamIdnt = num(9,:);
Level = mean(num(5:7,:)); %using final intensity values
famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlistDstn(i) = mean(cellscurr)/mean(Level);
   stdlist(i) = std(cellscurr)
end

meanlist = [];
stdlist = [];

[num]=xlsread('Families preliminary.xlsx',2); % RBPJ 

FamIdnt1 = num(6,:);
Level1 = num(5,:);

[num]=xlsread('Families preliminary.xlsx',3); % RBPJ 

FamIdnt2 = num(6,:);
Level2 = num(5,:);

FamIdnt = [FamIdnt1,FamIdnt2];
Level = [Level1,Level2];

famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlistJam2(i) = mean(cellscurr)/mean(Level);
   stdlist(i) = std(cellscurr)
end

meanlist = [];
stdlist = [];

[num]=xlsread('PGK families.xlsx',2); % Pgk 

FamIdnt = num(5,:);
Level = mean(num(2:4,:)); %using final intensity values
Level_N = zeros(size(Level));
famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlistPGK(i) = mean(cellscurr)/mean(Level);
   stdlist(i) = std(cellscurr)
end
meanlist = [];
stdlist = [];

Parmcurr = meanlistDstn;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s1 = [q5;q25;q75;q95];

Parmcurr = meanlistJam2;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s2 = [q5;q25;q75;q95];

Parmcurr = meanlistPGK;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s3 = [q5;q25;q75;q95];

s = [s1,s2,s3];

Parcurr = [meanlistDstn;meanlistJam2(1:24);meanlistPGK(1:24)];
 
Clonekey = {'Dstn','Jam2','Pgk'};

subplot(3,3,4)
h = colormap(hsv(size(Clonekey,2)));
h = h*0.6;
g = boxplot(Parcurr',char(Clonekey),'colors',h,'OutlierSize',7);

set(g(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(g(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(g(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(g(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(g(7,:),'Visible','off');
xtickangle(45)
ylim([0 2.5])
title('Family means','fontweight','normal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'d'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('Relative expression')

%%

[num]=xlsread('Dstn families modded.xlsx',3); % Dstn 

FamIdnt = num(9,:);
Level = mean(num(5:7,:)); %using final intensity values
famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistDstn(i) = stdlist(i)/meanlist(i); 
end

meanlist = [];
stdlist = [];

[num]=xlsread('Families preliminary.xlsx',2); % RBPJ 

FamIdnt1 = num(6,:);
Level1 = num(5,:);

[num]=xlsread('Families preliminary.xlsx',3); % RBPJ 

FamIdnt2 = num(6,:);
Level2 = num(5,:);

FamIdnt = [FamIdnt1,FamIdnt2];
Level = [Level1,Level2];

famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistJam2(i) = stdlist(i)/meanlist(i);
end

meanlist = [];
stdlist = [];

[num]=xlsread('PGK families.xlsx',2); % Pgk 

FamIdnt = num(5,:);
Level = mean(num(2:4,:)); %using final intensity values
Level_N = zeros(size(Level));
famidntstemp = zeros(size(FamIdnt));
Newlist = zeros(size(famidntstemp));

for i = 1:max(FamIdnt)
   famidntstemp =   (FamIdnt==i);
   cellscurr = Level(famidntstemp)
   meanlist(i) = mean(cellscurr)
   stdlist(i) = std(cellscurr)
   CVlistPGK(i) = stdlist(i)/meanlist(i);
   Level_N(famidntstemp) = Level(famidntstemp)-meanlist(i);
end
meanlist = [];
stdlist = [];


Parmcurr = CVlistDstn;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s1 = [q5;q25;q75;q95];

Parmcurr = CVlistJam2;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s2 = [q5;q25;q75;q95];

Parmcurr = CVlistPGK;

q5=quantile(Parmcurr,0.05);
q95=quantile(Parmcurr,0.95);
q75=quantile(Parmcurr,0.75);
q25=quantile(Parmcurr,0.25);
s3 = [q5;q25;q75;q95];

s = [s1,s2,s3];

Parcurr = [CVlistDstn;CVlistJam2(1:24);CVlistPGK(1:24)];
 
Clonekey = {'Dstn','Jam2','Pgk'};

subplot(3,3,5)
h = colormap(hsv(size(Clonekey,2)));
h = h*0.6;
g = boxplot(Parcurr',char(Clonekey),'colors',h,'OutlierSize',7);
set(g(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(g(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(g(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(g(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(g(7,:),'Visible','off');
xtickangle(45)
ylim([-0.1 1])
title('Within family CV','fontweight','normal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'e'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
ylabel('CV')

%%
% shows inference results

% ClonekeyN = {'Dstn_S1_STRSa1','Jam2_S1_STRSa1','Pgk_S1_STRSa1'};
% 
% for z3 = 1:length(ClonekeyN)
%     load([char(ClonekeyN(z3)),'_HMCprelim'])
%     ClonekeyN = {'Dstn_S1_STRSa1','Jam2_S1_STRSa1','Pgk_S1_STRSa1'};
%     thvec = (vertcat(chains{:}))';
%     mu_sigma = thvec(4,:);
%     mu_mu = thvec(3,:);
%     sigma_mu = exp(thvec(6,:));
%     sigma_sigma = exp(thvec(7,:));    
%     variance = (exp(sigma_mu.^2)-1).*exp(2*mu_mu+sigma_mu.^2);
%     meanlog_mu = exp(mu_mu+sigma_mu.^2/2);
%     meanlog_sigma = exp(mu_sigma+sigma_sigma.^2/2);    
%     
%     CV(z3,1:size(thvec,2)) = meanlog_sigma./meanlog_mu;
% 
%     mu_mu_list = exp(thvec(10:3:end,:));
%     mu_sigma_list = sqrt(exp(thvec(11:3:end,:)));
%     CV_list = mu_sigma_list./mu_mu_list;
%      
% end
% 
% subplot(3,3,6)
% h = colormap(hsv(size(ClonekeyN,2)));
% h = h*0.6;
% Clonekey = {'Dstn','Jam2','Pgk'};
% Parmcurr = CV';%squeeze(Parms(:,4,:));
% q5=quantile(Parmcurr,0.05);
% q95=quantile(Parmcurr,0.95);
% q75=quantile(Parmcurr,0.75);
% q25=quantile(Parmcurr,0.25);
% s = [q5;q25;q75;q95];
% g = boxplot(Parmcurr,char(Clonekey),'colors',h,'OutlierSize',7);
% set(g(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
% set(g(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
% set(g(3,:),{'Ydata'},num2cell(s([end end],:),1)')
% set(g(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
% set(g(7,:),'Visible','off');
% xtickangle(45)
% ylim([0 1])
% title('Transcriptional noise strength','fontweight','normal')
% a = xlim();
% b = ylim();
% text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'f'},...
%     'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
% ylabel('E[\sigma_{s,i}]/E[\mu_i]')

%%
subplot(3,3,7)
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'g'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

[num]=xlsread('Clones.xlsx',1); 

Clonekey = {'Dstn','Jam2','Pgk'};

idnt = ones(6,3);
idnt(:,2) = 2*idnt(:,2);
idnt(:,3) = 3*idnt(:,3);

idnt = idnt+0.2*rand(6,3);

idnt = reshape(idnt,1,18);
num = reshape(num,1,18);
%%
subplot(3,3,9)
scatter(idnt,num,'bx')

ylabel('Relative expression')
xlim([0.5 3.5])
ylim([0 2.7])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'h'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
xtickangle(45)
xticks([1 2 3])
xticklabels({'Dstn','Jam2','Pgk'})
