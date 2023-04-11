
% Resample for time warping functions

clear all; close all;
addpath('supplement\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: Random 500 warping functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsize = 16; % Label fontsize
nsize = 18; % Axis fontsize
d = 100;
t = linspace(0,1,d);
time_diff = mean(diff(t));
N = 500; % total number of simulated warping functions


rng(6677);
a = exprnd(1/3,[N,1]);
b = chi2rnd(3,[N,1]);
c = gamrnd(0.5,2,[N,1]);
d2 = exprnd(1/3,[N,1]);

for i =1:N
    r1(i,:) = (exp(a(i)*t)-1)/(exp(a(i))-1);
    r2(i,:) = (exp(b(i)*t)-1)/(exp(b(i))-1); 
    temp(i,:) = (exp(-c(i)*t)-1)/(exp(-c(i))-1);
    r3(i,:) = (exp(d2(i)*temp(i,:))-1)/(exp(d2(i))-1); 
    e = rand(2,1); e1 = e(1); e2 = max(e(2)-e(1),0); e3 = 1-e1-e2;
    r(i,:) = e1*r1(i,:)+e2*r2(i,:)+e3*r3(i,:);
    repo_e(i,:) = [e1,e2,e3];
end
figure(1); clf; % plot the simulated warping functions
h = plot(t,r);
axis equal;
ylim([0,1]);
xlim([0,1]);
set(gca, 'Fontsize', nsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
xticks([0 0.2 0.4 0.6 0.8 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Project warping function into H(0,1) space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1). transform to PDF space
for i = 1: N
     q(i,:) = repo_e(i,1)*(exp(a(i)*t)*a(i))/(exp(a(i))-1)+repo_e(i,2)*(exp(b(i)*t)*b(i))/(exp(b(i))-1)+repo_e(i,3)*(exp(d2(i)*temp(i,:))*d2(i))/(exp(d2(i))-1).*(exp(-c(i)*t)*(-c(i)))/(exp(-c(i))-1);  
end

% 2). transform to H(0,1) space
for i = 1:N
    v(i,:) = log(q(i,:))-trapz(t, log(q(i,:)));
end

% figure(2); clf; % plot the clr transformed warping functions
% plot(t,v);
% pbaspect([1 1 1]);
% set(gca, 'Fontsize', nsize,'linewidth', 1.5)
% set(gcf,'paperpositionmode','auto');
% set(gcf,'windowstyle','normal');
% set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
% set(gca,'fontweight','normal')
% opts.Colors     = get(groot,'defaultAxesColorOrder');
% opts.saveFolder = 'img/';
% opts.width      = 12;
% opts.height     = 10;
% opts.fontType   = 'Times';
% xticks([0 0.2 0.4 0.6 0.8 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3: fPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = mean(v);
C = cov(v);
[U, S, V] = svd(C);
S = S*time_diff;

dn =5;
for j = 1:dn
    U1(j,:) = U(:,j);
    U1(j,:) = U1(j,:)/sqrt(trapz(t, U1(j,:).^2)); 
end

U1_phi = exp(U1(1,:))./(trapz(t,exp(U1(1,:))));
U1_w = cumsum(U1_phi,2)./sum(U1_phi,2);
U1_w = normalize(U1_w,'range');
U2_phi = exp(U1(2,:))./(trapz(t,exp(U1(2,:))));
U2_w= cumsum(U2_phi,2)./sum(U2_phi,2);
U2_w = normalize(U2_w,'range');

figure(1); %plot the first and second pc
hold on;
FIG(1)=plot(t,U1_w, 'b', 'linewidth', 4);
FIG(2)=plot(t,U2_w, 'm', 'linewidth', 4);
legend(FIG([1 2]),{'First','Second'},'location','Southeast','FontSize',12,'Box','off');

%use pca for reconstruction get the coeff hist 
coeff1 = trapz(t, (v-mu).*U1(1,:),2);
coeff2 = trapz(t, (v-mu).*U1(2,:),2);

figure (4); clf;
h = histogram(coeff1,20);
h(1).FaceColor = [0 .5 .5];
pbaspect([1 1 1]);
set(gca, 'Fontsize', nsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

figure (5); clf;
h = histogram(coeff2,20);
h(1).FaceColor = [0 .5 .5];
pbaspect([1 1 1]);
set(gca, 'Fontsize', nsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

prolambda = S(1,1)/S(2,2);
pd1 = fitdist(coeff1,'Kernel');
pd2 = fitdist(coeff2,'Kernel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% Step 4: resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 500;
x_new2 = random(pd1,[n,1]).*U1(1,:)+random(pd2,[n,1]).*U1(2,:);
x_new2 = x_new2 +mu;
% figure (6);clf;
% plot(t, x_new2);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% set(gca, 'Fontsize', nsize,'linewidth', 1.5)
% set(gcf,'paperpositionmode','auto');
% set(gcf,'windowstyle','normal');
% set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
% set(gca,'fontweight','normal')
% opts.Colors     = get(groot,'defaultAxesColorOrder');
% opts.saveFolder = 'img/';
% opts.width      = 12;
% opts.height     = 10;
% opts.fontType   = 'Times';
% pbaspect([1 1 1]);

theta_B3 = exp(x_new2)./(trapz(t,exp(x_new2),2));
xnew_theta2= cumsum(theta_B3,2)./sum(theta_B3,2);
xnew_theta2=normalize(xnew_theta2','range');
figure (7);clf;
plot(t, xnew_theta2);
axis equal;
ylim([0,1]);
xlim([0,1]);
xticks([0 0.2 0.4 0.6 0.8 1]);
set(gca, 'Fontsize', nsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';

figure (8);clf;
for i =1:100
    eigenvalue(i) = S(i,i);
end
plot(cumsum(eigenvalue(1:10))/sum(eigenvalue),'linewidth', 1.5);
set(gca, 'Fontsize', nsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';
xlim([1,10]);
ylim([0 1]);
xticks([1 2 3 4 5 6 7 8 9 10]);
xlabel('Number of Principal Components');
ylabel('Fraction of explained variation')

axes('Position',[.3 .6 .3 .3])
box on
plot(cumsum(eigenvalue(1:5))/sum(eigenvalue),'linewidth', 1.5);
xlim([1,5]);
xticks([1 2 3 4 5]);
set(gca, 'Fontsize', lsize,'linewidth', 1.5)
set(gcf,'paperpositionmode','auto');
set(gcf,'windowstyle','normal');
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
set(gca,'fontweight','normal')
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 12;
opts.height     = 10;
opts.fontType   = 'Times';