function f = gmm_obj3(theta2)

% gmm objective function in random coefficients models
    global X1 Z 
    
    clear delta theta1 res1 mm sig g 

   
    
    delta= cverg_delta(theta2); 
    
    zz= Z'*Z;
    theta1= inv(X1'*Z*inv(zz)*Z'*X1)* X1'*Z*inv(zz)*Z'*delta;
    
    res1= delta-X1*theta1;
    mm= Z'*res1; % 10 by 1 matrix;
    
    sig= mm*mm';
    
    g= sum(mm,2); 
    
    
    

         f= 970*g'*inv(zz)*g;

    
    clear delta res1 mm sig g 

    

disp(['GMM objective:  ' num2str(f)])
disp(num2str(theta1))
