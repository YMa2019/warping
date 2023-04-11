% Resample for Berkey Growth Data of male

clear all; close all;
addpath('supplement\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: load the growth data and perform the necessary data clearning
%%% and smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load growth.mat;

rng(528);
lsize = 16; % Label fontsize
nsize = 18; % Axis fontsize

hboy = hgtmmat';
hgirl = hgtfmat';
age_new = linspace(1,18,101);
t = normalize(age_new, 'range');
time_diff = mean(diff(t));
age_nol = normalize(age, 'range');

% data smoothing
shrinking = 0.99;
[nm,mm] = size(hboy);
for i =1:nm
    for j =1:mm-1
        if hboy(i,j+1)<=hboy(i,j)
            hboy(i,j+1) = hboy(i,j) + mean(diff(hboy(i,j-6:j)))*(shrinking^j);
        end
    end
end
for i =1:nm
    hboyc(i,:) = csaps(age,hboy(i,:),.99,age_new);
    m_gam(i,:) = normalize(hboyc(i,:), 'range');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Project warping function into H(0,1) space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1). transform to PDF space
for i = 1: nm
    q_m(i,:) = gradient(m_gam(i,:),t);
end

% 2). transform to L2 space
for i = 1:nm
    v_m(i,:) = log(q_m(i,:))-trapz(t, log(q_m(i,:)));
end

% Calculate the growth Mean and plot it
mu_m = mean(v_m);
mu_m_tw = exp(mu_m)./(trapz(t,exp(mu_m)));
mu_m_tw= cumsum(mu_m_tw)./sum(mu_m_tw);
mu_m_tw = normalize(mu_m_tw,'range');
figure(1); clf;
plot(t,m_gam,'Color', [0.5 0.5 0.5]);
hold on; 
FIG(1)=plot(t, mu_m_tw,'b','linewidth',2);
legend(FIG([1]),{'Mean'},'location','Southeast','FontSize',12,'Box','off');
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
%%% Step 3: fPCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,d] = size(v_m);
C_m = cov(v_m);
[U_m, S_m, V_m] = svd(C_m);
S_m = S_m*time_diff;

dn =20; 
ln = 1;

% ln=dn;
for j = 1:dn
    Um(j,:) = U_m(:,j);
    Um(j,:) = Um(j,:)/sqrt(trapz(t, Um(j,:).^2)); 
end

dim=1; % change the dim here to plot the corresponding effect of the n^th eigenfunction as perturbation from the mean
c = -2:0.5:2;
for j =1:length(c)
    um_sigma(j,:) = c(j)*sqrt(S_m(dim,dim))*Um(dim,:)+ mu_m;
    U1_phi_m(j,:) = exp(um_sigma(j,:))./(trapz(t,exp(um_sigma(j,:)),2));
    U1_w_m(j,:)= cumsum(U1_phi_m(j,:),2)./sum(U1_phi_m(j,:),2);
    U1_w_m(j,:) = normalize(U1_w_m(j,:),'range');
end

figure(2); clf;
plot(t, U1_w_m,'linewidth',1.5);
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

%use pca for reconstruction get the coeff hist 
for i =1:dn
    coeff_m(i,:) = trapz(t, (v_m-mu_m).*Um(i,:),2);
    lb = min(coeff_m(i,:))-1e-5; ub = max(coeff_m(i,:))+1e-5;
    pd_m(i,:) = fitdist(coeff_m(i,:)','Kernel','Kernel','normal','Support',[lb,ub],'Width',0.05);
end

dim = dn;
lb = min(coeff_m(dim,:))-1e-5; ub = max(coeff_m(dim,:))+1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% Step 4: resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = N;
x_new2_m = 0;
for i =ln:dn
    re_coeff(i,:) = random(pd_m(i,:),[n,1]);
    temp = re_coeff(i,:)'.*Um(i,:);
    x_new2_m = x_new2_m+temp;
end
x_new2_m = x_new2_m+mu_m;

figure(3); clf;
subplot(1,2,1);
histogram(coeff_m(dim,:),20);
subplot(1,2,2);
histogram(re_coeff(dim,:),20);




theta_B3_m = exp(x_new2_m)./(trapz(t,exp(x_new2_m),2));
xnew_theta2_m= cumsum(theta_B3_m,2)./sum(theta_B3_m,2);
for i=1:n
    xnew_theta2_m(i,:) = (xnew_theta2_m(i,:)-min(xnew_theta2_m(i,:)))/(max(xnew_theta2_m(i,:))-min(xnew_theta2_m(i,:)));
end
figure (4);clf;
plot(t, xnew_theta2_m);
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

figure (5);clf;
for i =1:d
    eigenvalue_m(i) = S_m(i,i);
end
plot(cumsum(eigenvalue_m(1:25))/sum(eigenvalue_m),'linewidth', 1.5);
pbaspect([1 1 1]);
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
xlim([1,25]);
ylim([0,1]);
xlabel('Number of Principal Components');
ylabel('Fraction of explained variation')


