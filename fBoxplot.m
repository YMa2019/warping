% boxplot for warping function

clear all; close all;
lsize = 16; % Label fontsize
nsize = 18; % Axis fontsize
addpath('supplement\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: Simulate 500 warping functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default');
d = 100;
t = linspace(0,1,d);
time_diff = mean(diff(t));
N = 500;
for k_b =1:d
    f1(2*k_b-1,:) = sqrt(2)*sin(2*k_b*pi*t) ;
    f1(2*k_b,:) = sqrt(2)*cos(2*k_b*pi*t);
end

x_new2 = zeros(N,d);
%%Generate the stochastic process
x_new2 = zeros(N,d);
dn =2;
%laplace
for k =1: dn
    temp = laprnd(N,1,0,1/(2*k)).*f1(k,:);
    x_new2 = x_new2 +temp;
end

theta_B3 = exp(x_new2)./(trapz(t,exp(x_new2),2));
xnew_theta2= cumsum(theta_B3,2)./sum(theta_B3,2);
figure (1);clf;
plot(t, xnew_theta2,'linewidth', 1, 'Color', [0.5 0.5 0.5]);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: fPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = xnew_theta2;

% 1). transform to PDF space
for i =1:N
    q(i,:) = gradient(r(i,:), t);
end

% 2). transform to L2 space
for i = 1:N
    v(i,:) = log(q(i,:))-trapz(t, log(q(i,:)));
end


% 3). PCA
mu = mean(v);
v_cen = v - mu;
C = cov(v);
[U, S, V] = svd(C);
S = S*time_diff; 

figure (2);clf;
for i =1:10
    eigenvalue(i) = S(i,i);
end
plot(cumsum(eigenvalue(1:10))/sum(eigenvalue),'linewidth', 1.5);
xlim([1,10]);
ylim([0, 1]);
xlabel('Number of Principal Components')
ylabel('Fraction of explained variation')
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
box on;

dn =2;
for j = 1:dn
    U1(j,:) = U(:,j);
    U1(j,:) = U1(j,:)/sqrt(trapz(t, U1(j,:).^2)); 
end

% Calculate coefficients
c1 = trapz(t, v_cen.*U1(1,:),2);
c2 = trapz(t, v_cen.*U1(2,:),2);
P = [c1,c2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3: Calculate the depth for the coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depth = 1./(1 + mahal(P, P)); % for simplicity, we used the Mahalanobis Depth 

% To use the Tukey depth, please read the comment below
% Uncomment the line below to save the coefficients to matlab data file
% which can be utilized to calculate the Tukey depth in R.
% save('pca_d2_Gau_Half.dat','P','-ascii')

% Read the depth file calculated and saved use the R code provided as the
% name "Tukey Depth Calculation"
% depth = csvread('out_d2_Lap_Half.csv',1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 4: Construct the convex hull
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bag, Ind1] =  maxk(depth,N*0.5,1);
%[bag, Ind1] =  maxk(depth(:,dn+1),N*0.5,1); % uncomment for Tukey depth
%calculated by R code

Bag = P(Ind1,:);
rho = sqrt(chi2inv(0.99, dn))/sqrt(chi2inv(0.5, dn));
T = centroid(Bag(unique(convhulln(Bag)),:));
% T = mean(Bag(unique(convhulln(Bag)),:));
temp = (Bag - T)*rho+Bag;
med = 0;
for j = 1:2
    med = mean(T(:,j)).*U1(j,:)+med;
end

med = med + mu;
[k_b,av_b] = convhull(Bag);
intemp = inhull(P,temp,[],1e-10);
indtemp = find(intemp==1);
Loop = P(indtemp,:);
[k_l,av_l] = convhull(Loop);
in = inhull(P,Bag,[],1e-10);
indin = find(in==1);
in2 = inhull(P,Loop,[],1e-10);
ind2in = find(in2==1);
indnot= find(in2==0);

figure(3); clf;
axis square;
hold on
plot(Bag(k_b,1),Bag(k_b,2),'b','linewidth', 1.2);
plot(Loop(k_l,1),Loop(k_l,2),'c','linewidth', 1.2);
scatter(Loop(:,1), Loop(:,2), 10, 'c');
scatter(P(indin,1), P(indin,2),10, 'b');
scatter(T(1), T(2),30, 'o','k','filled')
scatter(P(indnot,1), P(indnot,2),20, 'r','d','filled');
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
box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5: Plot the boxplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mincuv = min(r(indin,:));
maxcuv = max(r(indin,:));
mincuv2 = min(r(indtemp,:));
maxcuv2 = max(r(indtemp,:));

figure(4); clf;
hold on;
patch([t, fliplr(t)], [mincuv2, fliplr(maxcuv2)],[0.8 0.8 0.8], 'EdgeColor',[0.8 0.8 0.8]);
patch([t, fliplr(t)], [mincuv, fliplr(maxcuv)],[0.65 0.65 0.65], 'EdgeColor',[0.65 0.65 0.65]);
plot(t, clr_inv(med,t),'k','LineWidth',3)
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
box on;
if length(indnot)>0
    figure(4); hold on;
    plot(t, r(indnot,:),'r--');
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
end




