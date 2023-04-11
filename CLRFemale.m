% Resample for Berkey Growth Data of female

clear all; close all;
addpath('supplement\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1: load the growth data and perform the necessary data clearning
%%% and smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load growth.mat;
rng(184);
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
[nf,mf] = size(hgirl);
for i =1:nf
    for j =1:mf-1
        if hgirl(i,j+1)<=hgirl(i,j)
            hgirl(i,j+1) = hgirl(i,j) + mean(diff(hgirl(i,j-6:j)))*(shrinking^j);
        end
    end
end

for i =1:nf
    hgirl_nol(i,:) = normalize(hgirl(i,:), 'range');
    hgirlc(i,:) = csaps(age,hgirl(i,:),.99,age_new);
    f_gam(i,:) = normalize(hgirlc(i,:), 'range');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2: Project warping function into H(0,1) space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1). transform to PDF space
for i = 1: nf
    q_f(i,:) = gradient(f_gam(i,:),t);
end

% 2). transform to L2 space
for i = 1:nf
    v_f(i,:) = log(q_f(i,:))-trapz(t, log(q_f(i,:)));
end

% Calculate the growth Mean and plot it
mu_f = mean(v_f);
mu_f_tw = exp(mu_f)./(trapz(t,exp(mu_f)));
mu_f_tw= cumsum(mu_f_tw)./sum(mu_f_tw);
mu_f_tw = normalize(mu_f_tw,'range');
figure(1); clf;
plot(t,f_gam,'Color', [0.5 0.5 0.5]);
hold on; 
FIG(1)=plot(t, mu_f_tw,'b','linewidth',2);
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
[N,d] = size(v_f);
C_f = cov(v_f);
[U_f, S_f, V_f] = svd(C_f);
S_f = S_f*time_diff;

dn =20; 
ln = 1;

% ln=dn;
for j = 1:dn
    Uf(j,:) = U_f(:,j);
    Uf(j,:) = Uf(j,:)/sqrt(trapz(t, Uf(j,:).^2)); 
end

dim=1;% change the dim here to plot the corresponding effect of the n^th eigenfunction as perturbation from the mean

c = -2:0.5:2;
for j =1:length(c)
    uf_sigma(j,:) = c(j)*sqrt(S_f(dim,dim))*Uf(dim,:)+ mu_f;
    U1_phi_f(j,:) = exp(uf_sigma(j,:))./(trapz(t,exp(uf_sigma(j,:)),2));
    U1_w_f(j,:)= cumsum(U1_phi_f(j,:),2)./sum(U1_phi_f(j,:),2);
    U1_w_f(j,:) = normalize(U1_w_f(j,:),'range');
end

figure(2); clf;
plot(t, U1_w_f,'linewidth',1.5);
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
    coeff_f(i,:) = trapz(t, (v_f-mu_f).*Uf(i,:),2);
    lb = min(coeff_f(i,:))-1e-5; ub = max(coeff_f(i,:))+1e-5;
    pd_f(i,:) = fitdist(coeff_f(i,:)','Kernel','Kernel','normal','Support',[lb,ub],'Width',0.02);
end

dim = dn;
lb = min(coeff_f(dim,:))-1e-5; ub = max(coeff_f(dim,:))+1e-5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% Step 4: resampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = N;
x_new2_f = 0;
for i =ln:dn
    re_coeff(i,:) = random(pd_f(i,:),[n,1]);
    temp = re_coeff(i,:)'.*Uf(i,:);
    x_new2_f = x_new2_f+temp;
end
x_new2_f = x_new2_f+mu_f;

figure(6); clf;
subplot(1,2,1);
histogram(coeff_f(dim,:),20);
subplot(1,2,2);
histogram(re_coeff(dim,:),20);


theta_B3_f = exp(x_new2_f)./(trapz(t,exp(x_new2_f),2));
xnew_theta2_f= cumsum(theta_B3_f,2)./sum(theta_B3_f,2);
for i=1:n
    xnew_theta2_f(i,:) = (xnew_theta2_f(i,:)-min(xnew_theta2_f(i,:)))/(max(xnew_theta2_f(i,:))-min(xnew_theta2_f(i,:)));
end
figure (4);clf;
plot(t, xnew_theta2_f);
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
    eigenvalue_f(i) = S_f(i,i);
end
plot(cumsum(eigenvalue_f(1:d))/sum(eigenvalue_f),'linewidth', 1.5);
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
