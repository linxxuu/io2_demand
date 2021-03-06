function f = gmm_obj(theta)

    global s s0 x1 x2 x3 p c1 c2 c3 

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
    f0= log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p;
    f1= x1.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 
    f2= x2.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 
    f3= x3.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 

    f4= c1.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);
    f5= c2.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);
    f6= c3.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);

    % supply side of equations : cost shifters are uncorrelated with
    % residual term of price equation 
    f7= p - (gamma0 + gamma1*c1 + gamma2*c2 + gamma3*c3)-1./(alpha*(1-s));
    f8= c1 .* (p- (gamma0 + gamma1*c1+ gamma2*c2 + gamma3*c3) - 1./(alpha*(1-s)));
    f9= c2 .* (p- (gamma0 + gamma1*c1+ gamma2*c2 + gamma3*c3) - 1./(alpha*(1-s)));
    f10=c3 .* (p- (gamma0 + gamma1*c1+ gamma2*c2 + gamma3*c3) - 1./(alpha*(1-s)));
    
    sig = [f0 f1 f2 f3 f4 f5 f6 f7 f8 f9 f10];
    sig = sig'*sig;
    
    g=[mean(f0) mean(f1) mean(f2) mean(f3) mean(f4) mean(f5) mean(f6) mean(f7) mean(f8) mean(f9) mean(f10)]';

    f= 970*g'*(sig\g);

    % f= 970*g'*g;
end

