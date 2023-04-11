% penalized registration example (diagonal case)
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: simulate f1 and f2 as the two functions to be aligned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 101;
t = linspace(0,1,N);
f1 = @(x) 6*(0.8).^(20*x).*cos(10*pi.*x-pi/4);
f2 = @(x) 5*(0.8).^(20*x).*sin(10*pi.*x);
time_gap = 1/(N-1);
a = [-0.5, 2];
f1 = f1(t)';
gamma2 = (exp(a(2)*t)-1)/(exp(a(2))-1);
f2 = 1.1*interp1(t,f2(t),gamma2);
f2 = f2';
figure(1); clf;
plot(t, f1, t, f2);
q1 = sign(gradient(f1)/time_gap).*sqrt(abs(gradient(f1)/time_gap));
q2 = sign(gradient(f2)/time_gap).*sqrt(abs(gradient(f2)/time_gap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Gradient Descent for Penalized function alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamb = 10; %% set penalty parameter
ini_gamma = t.^(0.2); %% set initial function
learnrate = 0.0002; %% set learning rate
Maxiter = 100;
% simulate the diagnal inverse covariance f_cov
t = t';
f_function1 = @(t) 0.025*(t+0.1);
f_function2 = @(t) 2.5*t;
x1 = linspace(1,100,101);
f_cov = f_function1(x1(1:60)');
f_cov(62:101,:)= f_function2(x1(61:100));

% calculate the optimal warping new_gamma
new_gamma = penaltyFA(f1, f2, t, f_cov, lamb, ini_gamma, learnrate, Maxiter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3: Align f2 to match f1 without penalty by setting lamb = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamb = 0;
learnrate = 0.0003;
gamma_t = penaltyFA(f1, f2, t, f_cov, lamb, ini_gamma, learnrate, Maxiter);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 4: plot the alignment results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lsize = 16; % Label fontsize
nsize = 18; % Axis fontsize
figure(7);clf;
plot(t, gamma_t,'R','LineWidth', 1.5);
hold on;
plot(t, new_gamma,'K','LineWidth', 1.5);
plot(t, t ,'G','LineWidth', 1.5);
legend({'\gamma*_{\lambda=0}','\gamma*_{\lambda=10}','\gamma_{id}'},'location','best','FontSize',12,'Box','off');
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

%%% plot f1,f2,f2_gamma_t %%%
f2_gamma_t = interp1(t,f2,gamma_t);
figure(5);clf;
plot(t,f1,'b','LineWidth', 1.5);
hold on;
plot(t,f2,'g','LineWidth', 1.5);
plot(t, f2_gamma_t, 'r--','LineWidth', 1.5)
f2_gamma_t_new1 = interp1(t,f2,new_gamma);   
plot(t, f2_gamma_t_new1, 'k--','LineWidth', 1.5)
legend({'f_1(t)','f_2(t)','f_2(\gamma*_{\lambda=0}(t))','f_2(\gamma*_{\lambda=10}(t))'},'location','best','FontSize',12,'Box','off');
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
% pbaspect([1 1 1]);

figure(2);clf;
plot(t(1:61),f_cov(1:61),'b','LineWidth', 1.5);
hold on;
plot(t(62),f_cov(62),'bo','LineWidth', 1.5)
plot(t(62:101),f_cov(62:101),'b','LineWidth', 1.5);
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