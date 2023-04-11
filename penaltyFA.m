function new_gamma = penaltyFA(f1, f2, t, f_cov, lamb, ini_gamma, learnrate, Maxiter)
    lamb = lamb; %% set penalty parameter
    ini_gamma = ini_gamma'; %% set initial function
    learnrate = learnrate; %% set learning rate
    Maxiter = Maxiter;
    
    N = size(f1,1);
    time_gap = 1/(N-1);
    q1 = sign(gradient(f1)/time_gap).*sqrt(abs(gradient(f1)/time_gap));
    q2 = sign(gradient(f2)/time_gap).*sqrt(abs(gradient(f2)/time_gap));

    phi = log(gradient(ini_gamma)/time_gap); %Calculate corresponding "initial phi=log(dot_gamma)"
    dot_q2 = gradient(q2)/time_gap; %Calculate derivative of q2
    
    cnt = 1;
    phi_repo = []; %save all the updated phi's
    phi_repo(:,1) = phi;
%     check1 = trapz(t, exp(phi)); %check if the int_0^1(exp(phi)(t)dt=1
    
    temp_in = cumtrapz(t,exp(phi));
    temp_in = round(temp_in/temp_in(end)*(N-1))+1;
    loss = trapz(t,-2*q1.*q2(temp_in).*sqrt(exp(phi)))...
           +lamb*trapz(t,f_cov.*phi.^2)+lamb*trapz(t,f_cov)*(trapz(t,phi))^2-2*lamb*(trapz(t,f_cov.*phi))*(trapz(t,phi)); 
    
    temp_t = temp_in;
    
    
    
    for iter = 1:Maxiter
        learn_rate = learnrate*(1-0.001)^(iter-1); %set learning rate Decay
    
        f_grad = trapz(t, q1.*dot_q2(temp_t).*sqrt(exp(phi)))-cumtrapz(t, q1.*dot_q2(temp_t).*sqrt(exp(phi)));    
        f_grad = -2*f_grad.*exp(phi) - q1.*q2(temp_t).*sqrt(exp(phi));
        f_grad = f_grad + 2*lamb*(f_cov.*phi + trapz(t, phi)*trapz(t, f_cov) - trapz(t,phi).*f_cov - trapz(t, f_cov.*phi));
    
        phi = phi-learn_rate.*f_grad;
    
        phi = phi-log(trapz(t,exp(phi)));
    
        cnt = cnt+1;
%         check1(cnt) = trapz(t, exp(phi));
        phi_repo(:,cnt) = phi;
    
        temp_t = cumtrapz(t,exp(phi));
        temp_t = round(temp_t/temp_t(end)*(N-1))+1;
        loss(cnt) = trapz(t,-2*q1.*q2(temp_t).*sqrt(exp(phi)))+lamb*trapz(t,f_cov.*phi.^2)...
                    +lamb*trapz(t,f_cov)*(trapz(t,phi))^2-2*lamb*(trapz(t,f_cov.*phi))*(trapz(t,phi));
    end
    
    new_gamma = cumsum(exp(phi_repo(:,end)))./sum(exp(phi_repo(:,end)));
    new_gamma = (new_gamma-min(new_gamma))/(max(new_gamma)-min(new_gamma));
end 