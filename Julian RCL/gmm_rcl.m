function f = gmm_rcl(sigg)

% gmm objective function, random coefficients model (RCL)
    global X1 Z s s0
    %delta_save s s0
    
    %clear delta theta1 res1 mm sig g 

   
    delta=ones(970,1);
    delta=log(s) - log(s0);
    delta_fixpoint= delta_fp(sigg,delta); 
    delta=delta_fixpoint;
    zz= Z'*Z;
    theta1= inv(X1'*Z*inv(zz)*Z'*X1)* X1'*Z*inv(zz)*Z'*delta;
    
    res1= delta-X1*theta1;
    mm= Z'*res1; % this is a 10 by 1 matrix;
    
    sig= mm*mm';
    
    g= sum(mm,2); 
    
    
    

    f= 970*mm'/zz*mm;
 % f=g'*g;

    
    %clear delta res1 mm sig g 

    

disp(['GMM objective:  ' num2str(f)])
disp(num2str(theta1))
%disp(num2str(sigma))
disp('Note: 2nd entry is price coeff.')
disp(num2str(sigg))
end