
clear 
% load in data
% [num]=csvread('4.12 with mitosis_50p.csv'); %Dstn
% clonename =  'Dstn';
% num=csvread('4.2 with mitosis_50pairs.csv'); % Nono
% clonename =  'Nono';
% [num]=csvread('4.11 with mitosis.csv'); % Jarid2
% clonename =  'Jarid2';
% [num]=csvread('4.10 with mitosis_50.csv'); % Rbpj
% clonename =  'Rbpj';
[num]=csvread('pgk21_all.csv'); % pgk21 
% clonename =  'pgk21';
% [num]=csvread('8.2 with mitosis_50p.csv'); % Sprouty4
% % clonename =  'Sprouty';
% [num]=xlsread('4.9 50 pairs.xlsx'); % Jam2
% clonename =  'Jam2';
% [num]=xlsread('NC7.xlsx'); % Rasa12
% clonename =  'Rasa12';
% [num]=xlsread('PGK 19 all_new.xlsx')
num(isnan(num))=0; %replace nans by zeros
num(num<0) = 0;
num1 = zeros(100,size(num,2)-18);
for w = 1:size(num,2)
    y1 = num(:,w);
    fin1 = find(y1,1,'last')-15; % approx last 15 time points - defocussing of image during mitosis
    y1 = y1(4:fin1);
    num1(1:length(y1),w) = y1;
end

% normalise signal over whole cell cycle

t = (0:5:(5*size(num1,1)-1))'; % define time points (in mins)

normalised_y = zeros(100,size(num1,2));

for i = 1:size(num1,2)
     y = num1(:,i);
     ty = t;
     ty(y==0) = [];
     y(y==0) = [];
     xq = linspace(0,ty(end),100);
     vq2 = interp1(ty,y,xq,'spline'); 
     normalised_y(:,i) = vq2;
end

% 
%%
% rank cells in population according to initial value
rowstocol = normalised_y';
[a,b] = sort(rowstocol,1);
[w,l] = sort(b);
[q,t] = sort(l(:,1));
trial = l(t,:);

figure()

[q1,t1] = sort(rowstocol(:,1));
trial1 = rowstocol(t1,:);

subplot(2,4,1)

h = colormap(jet(size(trial,1)));

v = randperm(size(trial1,1));
for j = 1:size(trial1,1)
    hold on
    i= v(j);
    plot(1:size(normalised_y,1),trial1(i,:),'color',h(i,:))
end
plot(1:size(normalised_y,1),mean(trial1),'k--','linewidth',2)
hold off
xlabel('% cell cycle')
ylabel('Bioluminescence [A.U]')
ylim([0 25000])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'a'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
title('pGK','fontweight','normal')
yticks([0 10000 20000])
yticklabels({'0','10000','20000'})

%%
subplot(2,4,2)

h = colormap(jet(size(trial,1)));
    
hold on
for p = 1:3
    plot(1:size(normalised_y,1),trial1(p,:),'color',h(p,:))
end

for p = 0:2
    plot(1:size(normalised_y,1),trial1(end-p,:),'color',h(end-p,:))
end

hold off
xlabel('% cell cycle')
ylabel('Bioluminescence [A.U]')
ylim([0 25000])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'b'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
yticks([0 10000 20000])
yticklabels({'0','10000','20000'})

%%

t = (0:5:(5*size(num,1)-1))'/60; % define time points (in mins)

yTOT = num;
cell = 1;
subplot(2,4,3)
hold on

y1 = yTOT(:,cell);
y2 = yTOT(:,cell+1);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1a = y1(1:fin);
y2a = y2(1:fin);
x = t(1:fin);
p1 = plot(x,[y1a],'r')
p1a = plot(x,[y2a],'r')
y1 = yTOT(:,cell+6);
y2 = yTOT(:,cell+7);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1b = y1(1:fin);
y2b = y2(1:fin);
x = t(1:fin);
p2 = plot(x,[y1b],'g')
p2a = plot(x,[y2b],'g')
y1 = yTOT(:,cell+4);
y2 = yTOT(:,cell+5);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1c = y1(1:fin);
y2c = y2(1:fin);
x = t(1:fin);
p3 = plot(x,[y1c],'b')
p3a = plot(x,[y2c],'b')
hold off
xlabel('Time [hours]')
ylabel('Bioluminescence [A.U]')

xlim([0 12])
yticks([0 10000 20000])
yticklabels({'0','10000','20000'})
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'c'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

%%
% calculate evolution of correlation over the cell cycle
corrvec = zeros(1,size(normalised_y,1));
firstcells = normalised_y(:,1:2:size(normalised_y,1));
secondcells = normalised_y(:,2:2:size(normalised_y,1));

for i = 1:size(normalised_y,1)
    corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
end

%%
% do bootstrap error on sister sister correlations
% choose pair randomly
% choose ordering randomly
bootsM = zeros(size(normalised_y));
pairnumsTOT = round(size(normalised_y,2)/2);
repeats = 100;

for i = 1:repeats
    r = randi([1 pairnumsTOT],pairnumsTOT,1);
    for j = 1:length(r)
        cell1 = normalised_y(:,2*r(j)-1);
        cell2 = normalised_y(:,2*r(j));
        p = rand;
        if p<0.5
            bootsM(:,2*j-1) = cell1;
            bootsM(:,2*j) = cell2;
        else
            bootsM(:,2*j-1) = cell2;
            bootsM(:,2*j) = cell1;
        end
    end
    
    firstcells = bootsM(:,1:2:size(bootsM,1));
    secondcells = bootsM(:,2:2:size(bootsM,1));
    
    for k = 1:size(normalised_y,1)
        corrvecB(i,k) = corr(firstcells(k,:)',secondcells(k,:)');    
    end
    
end

subplot(2,4,4)
xtrial = 1:100;
hold on
p1 = plot(xtrial,(corrvec),'g')
plot(xtrial,(corrvec)+std(corrvecB,[],1),'g--')
plot(xtrial,(corrvec)-std(corrvecB,[],1),'g--')

% find pairs of cells with nearest initial levels of reporter, excluding
% sisters
rowstocol = normalised_y';
[q1,t1] = sort(rowstocol(:,1));
sortedM = zeros(size(rowstocol));
for i = 1:2:size(rowstocol,1)
    sortedM(i,:) = rowstocol(i,:);
    I = find(t1==i);
    if t1(I+1) ~=  i+1
        sortedM(i+1,:) = rowstocol(t1(I+1),:);
    else
        sortedM(i+1,:) = rowstocol(t1(I-1),:);
        disp(i)
    end

end

[a,b] = sort(sortedM,1);
b1 = b';
[w,l] = sort(b);
l3 = l;
[q,t] = sort(l(:,1));
trial = l(t,:);

input = l;

normalised_y = input';
corrvec = zeros(1,size(normalised_y,1));
firstcells = normalised_y(:,1:2:size(normalised_y,1));
secondcells = normalised_y(:,2:2:size(normalised_y,1));

for i = 1:size(normalised_y,1)
    corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
end

% get bootstrap error
bootsM = zeros(size(normalised_y));
pairnumsTOT = round(size(normalised_y,2)/2);

for i = 1:repeats
    r = randi([1 pairnumsTOT],pairnumsTOT,1);
    for j = 1:length(r)
        cell1 = normalised_y(:,2*r(j)-1);
        cell2 = normalised_y(:,2*r(j));
        p = rand;
        if p<0.5
            bootsM(:,2*j-1) = cell1;
            bootsM(:,2*j) = cell2;
        else
            bootsM(:,2*j-1) = cell2;
            bootsM(:,2*j) = cell1;
        end
    end
    
    firstcells = bootsM(:,1:2:size(bootsM,1));
    secondcells = bootsM(:,2:2:size(bootsM,1));
    
    for k = 1:size(normalised_y,1)
        corrvecB(i,k) = corr(firstcells(k,:)',secondcells(k,:)');    
    end
    
end


xtrial = 1:100;
p2 = plot(xtrial,(corrvec),'r')
plot(xtrial,(corrvec)+std(corrvecB,[],1),'r--')
plot(xtrial,(corrvec)-std(corrvecB,[],1),'r--')
hold off
ylim([0 1])
xlabel('% cell cycle')
ylabel('Corr')
legend([p1 p2],{'Sisters','Randomised'},'location','southwest')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'d'},...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

%%
clear 
% load in data
[num]=csvread('4.12 with mitosis_50p.csv'); %Dstn
% clonename =  'Dstn';
% num=csvread('4.2 with mitosis_50pairs.csv'); % Nono
% clonename =  'Nono';
% [num]=csvread('4.11 with mitosis.csv'); % Jarid2
% clonename =  'Jarid2';
% [num]=csvread('4.10 with mitosis_50.csv'); % Rbpj
% clonename =  'Rbpj';
% [num]=csvread('pgk21_all.csv'); % pgk21 
% clonename =  'pgk21';
% [num]=csvread('8.2 with mitosis_50p.csv'); % Sprouty4
% % clonename =  'Sprouty';
% [num]=xlsread('4.9 50 pairs.xlsx'); % Jam2
% clonename =  'Jam2';
% [num]=xlsread('NC7.xlsx'); % Rasa12
% clonename =  'Rasa12';
% [num]=xlsread('PGK 19 all_new.xlsx')
num(isnan(num))=0; %replace nans by zeros
num(num<0) = 0;

num1 = zeros(100,size(num,2)-18);
for w = 1:size(num,2)
    y1 = num(:,w);
    fin1 = find(y1,1,'last')-15;
    y1 = y1(4:fin1);
    num1(1:length(y1),w) = y1;
end

% normalise signal over whole cell cycle

t = (0:5:(5*size(num1,1)-1))'; % define time points (in mins)

normalised_y = zeros(100,size(num1,2));

for i = 1:size(num1,2)
     y = num1(:,i);
     ty = t;
     ty(y==0) = [];
     y(y==0) = [];
     xq = linspace(0,ty(end),100);
     vq2 = interp1(ty,y,xq,'spline'); 
     normalised_y(:,i) = vq2;
end

%%

ranked = rank(normalised_y);
rowstocol = normalised_y';
[a,b] = sort(rowstocol,1);
b1 = b';
[w,l] = sort(b);
[q,t] = sort(l(:,1));
trial = l(t,:);

R = trial;
maxlags=size(R,2)-1;
maxtime=size(R,2);
R2 = trial;
n = size(R2,2);
numcells = size(R2,1);

[q1,t1] = sort(rowstocol(:,1));
trial1 = rowstocol(t1,:);

subplot(2,4,5)
h = colormap(jet(size(trial,1)));
v = randperm(size(trial1,1));
for j = 1:size(trial1,1)
    hold on
    i= v(j);
    plot(1:size(normalised_y,1),trial1(i,:),'color',h(i,:))
end
plot(1:size(normalised_y,1),mean(trial1),'k--','linewidth',2)
hold off
xlabel('% cell cycle')
ylabel('Bioluminescence [A.U]')
a = xlim();
b = ylim();
title('Dstn','fontweight','normal')
yticks([0 4000 8000])
yticklabels({'0','4000','8000'})

%%

subplot(2,4,6)
h = colormap(jet(size(trial,1)));
    
hold on
for p = 1:3
    plot(1:size(normalised_y,1),trial1(p,:),'color',h(p,:))
end

for p = 0:2
    plot(1:size(normalised_y,1),trial1(end-p,:),'color',h(end-p,:))
end
hold off
xlabel('% cell cycle')
ylabel('Bioluminescence [A.U]')
a = xlim();
b = ylim();
yticks([0 4000 8000])
yticklabels({'0','4000','8000'})
%%

t = (0:5:(5*size(num,1)-1))'/60; % define time points (in hours)

yTOT = num;%normalised_y;
cell = 1;
subplot(2,4,7)
hold on

y1 = yTOT(:,cell);
y2 = yTOT(:,cell+1);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1a = y1(1:fin);
y2a = y2(1:fin);
x = t(1:fin);
p1 = plot(x,[y1a],'b')
p1a = plot(x,[y2a],'b')
y1 = yTOT(:,cell+2);
y2 = yTOT(:,cell+3);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1b = y1(1:fin);
y2b = y2(1:fin);
x = t(1:fin);
p2 = plot(x,[y1b],'g')
p2a = plot(x,[y2b],'g')
y1 = yTOT(:,cell+4);
y2 = yTOT(:,cell+5);
fin = min(find(y1,1,'last'),find(y2,1,'last'));
y1c = y1(1:fin);
y2c = y2(1:fin);
x = t(1:fin);
p3 = plot(x,[y1c],'r')
p3a = plot(x,[y2c],'r')
hold off
xlabel('Time [hours]')
ylabel('Bioluminescence [A.U]')
xlim([0 12])
a = xlim();
b = ylim();
yticks([0 3000 6000])
yticklabels({'0','3000','6000'})

%%

corrvec = zeros(1,size(normalised_y,1));
firstcells = normalised_y(:,1:2:size(normalised_y,1));
secondcells = normalised_y(:,2:2:size(normalised_y,1));

for i = 1:size(normalised_y,1)
    corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
end

%%
% do bootstrap error on sister sister correlations
% choose pair randomly
% choose ordering randomly
bootsM = zeros(size(normalised_y));
pairnumsTOT = round(size(normalised_y,2)/2);
repeats = 100;

for i = 1:repeats
%     i
    r = randi([1 pairnumsTOT],pairnumsTOT,1);
    for j = 1:length(r)
        cell1 = normalised_y(:,2*r(j)-1);
        cell2 = normalised_y(:,2*r(j));
        p = rand;
        if p<0.5
            bootsM(:,2*j-1) = cell1;
            bootsM(:,2*j) = cell2;
        else
            bootsM(:,2*j-1) = cell2;
            bootsM(:,2*j) = cell1;
        end
    end
    
    firstcells = bootsM(:,1:2:size(bootsM,1));
    secondcells = bootsM(:,2:2:size(bootsM,1));
    
    for k = 1:size(normalised_y,1)
        corrvecB(i,k) = corr(firstcells(k,:)',secondcells(k,:)');    
    end
    
end

subplot(2,4,8)
xtrial = 1:100;
hold on
p1 = plot(xtrial,(corrvec),'g')
plot(xtrial,(corrvec)+std(corrvecB,[],1),'g--')
plot(xtrial,(corrvec)-std(corrvecB,[],1),'g--')

rowstocol = normalised_y';
[q1,t1] = sort(rowstocol(:,1));
trial1 = rowstocol(t1,:);
sortedM = zeros(size(rowstocol));
for i = 1:2:size(rowstocol,1)
    sortedM(i,:) = rowstocol(i,:);
    I = find(t1==i);
    if (I+1)<size(rowstocol,1) && t1(I+1) ~=  i+1 
        sortedM(i+1,:) = rowstocol(t1(I+1),:);
    else
        sortedM(i+1,:) = rowstocol(t1(I-1),:);
        disp(i)
    end

end

[a,b] = sort(sortedM,1);
b1 = b';
[w,l] = sort(b);
l3 = l;
[q,t] = sort(l(:,1));
trial = l(t,:);
input = l;

normalised_y = input';
corrvec = zeros(1,size(normalised_y,1));
firstcells = normalised_y(:,1:2:size(normalised_y,1));
secondcells = normalised_y(:,2:2:size(normalised_y,1));

for i = 1:size(normalised_y,1)
    corrvec(i) = corr(firstcells(i,:)',secondcells(i,:)');    
end

bootsM = zeros(size(normalised_y));
pairnumsTOT = round(size(normalised_y,2)/2);

for i = 1:repeats
    r = randi([1 pairnumsTOT],pairnumsTOT,1);
    for j = 1:length(r)
        cell1 = normalised_y(:,2*r(j)-1);
        cell2 = normalised_y(:,2*r(j));
        p = rand;
        if p<0.5
            bootsM(:,2*j-1) = cell1;
            bootsM(:,2*j) = cell2;
        else
            bootsM(:,2*j-1) = cell2;
            bootsM(:,2*j) = cell1;
        end
    end
    
    firstcells = bootsM(:,1:2:size(bootsM,1));
    secondcells = bootsM(:,2:2:size(bootsM,1));
    
    for k = 1:size(normalised_y,1)
        corrvecB(i,k) = corr(firstcells(k,:)',secondcells(k,:)');    
    end
    
end


xtrial = 1:100;
p2 = plot(xtrial,(corrvec),'r')
plot(xtrial,(corrvec)+std(corrvecB,[],1),'r--')
plot(xtrial,(corrvec)-std(corrvecB,[],1),'r--')
hold off

xlabel('% cell cycle')
ylabel('Corr')
ylim([-0.2 1])
a = xlim();
b = ylim();

