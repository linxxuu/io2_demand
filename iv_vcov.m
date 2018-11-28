function sigma = iv_vcov(theta)
    % use delta method to find VCOV 
    % 1. construct D : rows = moments, cols = vars
    %   get the mean derivative (=sample avg)
    % 2. get asymptotic variance of the efficient GMM
    
    global s s0 sh_g x1 x2 x3 p c1 c2 c3 
    
    b0=theta(1);
    alpha= theta(2);
    b1=theta(3);
    b2= theta(4);
    b3= theta(5);

    
    %% 1. construct D
    D= zeros(7,5);
    
    % demand side of equations: cost shifters and x are exgonous to residual of demand equations 
    
    f0= log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p;
    D(1,1) = -1;
    D(1,2) = mean(-p);
    D(1,3) = mean(-x1);
    D(1,4) = mean(-x2);
    D(1,5) = mean(-x3);
    f1= x1.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 
    D(2,1) = mean(-x1);
    D(2,2) = mean(-x1.*p);
    D(2,3) = mean(-x1.*x1);
    D(2,4) = mean(-x1.*x2);
    D(2,5) = mean(-x1.*x3);
    f2= x2.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 
    D(3,1) = mean(-x2);
    D(3,2) = mean(-x2.*p);
    D(3,3) = mean(-x2.*x1);
    D(3,4) = mean(-x2.*x2);
    D(3,5) = mean(-x2.*x3);
    f3= x3.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p); 
    D(4,1) = mean(-x3);
    D(4,2) = mean(-x3.*p);
    D(4,3) = mean(-x3.*x1);
    D(4,4) = mean(-x3.*x2);
    D(4,5) = mean(-x3.*x3);
    
    f4= c1.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);
    D(5,1) = mean(-c1);
    D(5,2) = mean(-c1.*p);
    D(5,3) = mean(-c1.*x1);
    D(5,4) = mean(-c1.*x2);
    D(5,5) = mean(-c1.*x3);
    f5= c2.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);
    D(6,1) = mean(-c2);
    D(6,2) = mean(-c2.*p);
    D(6,3) = mean(-c2.*x1);
    D(6,4) = mean(-c2.*x2);
    D(6,5) = mean(-c2.*x3);
    f6= c3.*(log(s)-log(s0)- b0- b1*x1- b2*x2- b3*x3- alpha*p);
    D(7,1) = mean(-c3);
    D(7,2) = mean(-c3.*p);
    D(7,3) = mean(-c3.*x1);
    D(7,4) = mean(-c3.*x2);
    D(7,5) = mean(-c3.*x3);
    
    %% 2. get asymptotic vcov of efficient GMM
    
    eps = [f0 f1 f2 f3 f4 f5 f6];
    S = eps'*eps;
    
    sigma = 970*(D'*(S\D))\eye(size(theta,1));
    
end
