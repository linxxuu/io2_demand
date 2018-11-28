function f = gmm_obj(theta)

    global sh_g0 sh_g x1_ng x2_ng x3_ng p_ng c1_ng c2_ng c3_ng 
%     theta = param_gmm;
    b0=theta(1);
    alpha= theta(2);
    b1=theta(3);
    b2= theta(4);
    b3= theta(5);
    gamma0= theta(6);
    gamma1= theta(7);
    gamma2= theta(8);
    gamma3=theta(9);
    
    % demand side of equations: cost shifters and x are exgonous to residual of demand equations 
    f0= log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng;
    f1= x1_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng); 
    f2= x2_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng); 
    f3= x3_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng); 

    f4= c1_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng);
    f5= c2_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng);
    f6= c3_ng.*(log(sh_g0)- b0- b1*x1_ng- b2*x2_ng- b3*x3_ng- alpha*p_ng);

    % supply side of equations : cost shifters are uncorrelated with
    % residual term of price equation 
    f7= p_ng - (gamma0 + gamma1*c1_ng + gamma2*c2_ng + gamma3*c3_ng)-1./(alpha*(1-sh_g));
    f8= c1_ng .* (p_ng- (gamma0 + gamma1*c1_ng+ gamma2*c2_ng + gamma3*c3_ng) - 1./(alpha*(1-sh_g)));
    f9= c2_ng .* (p_ng- (gamma0 + gamma1*c1_ng+ gamma2*c2_ng + gamma3*c3_ng) - 1./(alpha*(1-sh_g)));
    f10=c3_ng .* (p_ng- (gamma0 + gamma1*c1_ng+ gamma2*c2_ng + gamma3*c3_ng) - 1./(alpha*(1-sh_g)));
    
    sig = [f0 f1 f2 f3 f4 f5 f6 f7 f8 f9 f10];
    sig = sig'*sig;
       
    g=[mean(f0) mean(f1) mean(f2) mean(f3) mean(f4) mean(f5) mean(f6) mean(f7) mean(f8) mean(f9) mean(f10)]';
    
%     disp(size(sig));
%     disp(size(g));
%     disp(size(sig\g));
    f= 970*g'*(sig\g);
    
    % f= 970*g'*g;
end

