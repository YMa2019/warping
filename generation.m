% Generate the time warping functions without prior information using
% Algorithm 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: Create Fourier basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clear all;
addpath('supplement\');
rng(20230410); % change the seed


d = 100;
t = linspace(0,1,d);
n = 10;
for k =1:d
    f1(2*k-1,:) = sqrt(2)*sin(2*k*pi*t) ;
    f1(2*k,:) = sqrt(2)*cos(2*k*pi*t);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Generate the stochastic process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_new2 = zeros(n,d);
dn =20; % the number of basis function we used

% 1). Gaussian process
for k =1: dn
    temp = normrnd(0,1/(1*k),[n,1]).*f1(k,:);
    x_new2 = x_new2 +temp;
end

% 2). laplace  %uncomment line 33 - 36 to get laplacian distributed
% coefficients
% for k =1: dn
%     temp = laprnd(n,1,0,1/k).*f1(k,:);
%     x_new2 = x_new2 +temp;
% end


% 3). uniform %uncomment line 41 - 45 to get laplacian distributed
% coefficients
% c = sqrt(3);
% for k =1: dn
%     temp = unifrnd(-c/k,c/k,[n,1]).*f1(k,:);
%     x_new2 = x_new2 +temp;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3: Plot the simulated time warping functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsize = 16; % Label fontsize
nsize = 18; % Axis fontsize

figure (5);clf;
plot(t, x_new2,'linewidth', 1); % plot the stochastic process in H(0, 1)
xticks([0 0.2 0.4 0.6 0.8 1]);
ylim([-4,4]);
% xlim([0,1]);
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


% apply the inverse-clr transfromation to get warping function
theta_B3 = exp(x_new2)./(trapz(t,exp(x_new2),2));
xnew_theta2= cumsum(theta_B3,2)./sum(theta_B3,2);
figure (6);clf;
plot(t, xnew_theta2,'linewidth', 1);% plot the simulated warping function
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


