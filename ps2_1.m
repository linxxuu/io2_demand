% PS2 Demand

%% Input Data
clc
clear all

%M= csvread('Data.csv');
Data= dlmread('Data.csv');

mid = Data(:, 1);
j = Data(:, 2);
s = Data(:, 3);
x1 = Data(:, 4);
x2 = Data(:, 5);
x3 = Data(:, 6);
p = Data(:, 7);
w_cost1= Data(:, 8);
w_cost2= Data(:, 9);
w_cost3= Data(:, 10);
group_id= Data(:, 11);


% calculate outside options: 
s0= zeros(970,1);
for i=1:970
    
    s0(i)= 1- (mid(:) == mid(i))'*s;
    
end
 
    
%% 1.1 Pooled OLS 

Y= log(s)-log(s0);
X= [ones(970,1) p x1 x2 x3];

bhat = (X'*X)\(X'*Y);
se1ols= sqrt(diag(mean((Y-X*bhat).^2)*((X'*X)\eye(size(X,2)))));
vars = ['Const     ';'Price     ';'X1        ';'X2        ';'X3        '];
str = [bhat(1:5,1) se1ols(1:5,1)];
 disp('**********************************************');
 disp('Pooled OLS Estimation in Logit Model:');
 disp('**********************************************');
 disp([' Vars','       Coeff','    ','Std Err']);
 disp([vars,num2str(str)]);
 disp(' ');
 disp(' ');
 disp(' ');
 
 %% 1.2.1 IV
 
 % first stage: cost shifters and product characteristics 
 %Question: whether including constant term in first stage? 

 Z= [ones(970,1) w_cost1 w_cost2 w_cost3 x1 x2 x3];
 bhatFS = (Z'*Z)\(Z'*p);
 se1FS= sqrt(diag(mean((p-Z*bhatFS).^2)*((Z'*Z)\eye(size(Z,2)))));
 vars = ['const     ';'w1        ';'w2        ';'w3        ';'X1        ';'X2        ';'X3        '];
 str = [bhatFS(1:7,1) se1FS(1:7,1)];
 disp('**********************************************');
 disp('IV: First Stage Estimation in Logit Model:');
 disp('**********************************************');
 disp([' Vars','       Coeff','    ','Std Err']);
 disp([vars,num2str(str)]);
 disp(' ');
 disp(' ');
 disp(' ');
 % need to calculate F and p 
 
 
 phat= Z*bhatFS;
 X= [ones(970,1) p x1 x2 x3];
 Y= log(s);
 Xhat= [ones(970,1) phat x1 x2 x3];
 bhat2SLS = (Xhat'*X)\(Xhat'*Y); 
 se2SLS= sqrt(diag(mean((Y-X*bhat2SLS).^2)*((X'*X)\eye(size(X,2)))));
 vars = ['Const     ';'Price     ';'X1        ';'X2        ';'X3        '];
 str = [bhat2SLS(1:5,1) se2SLS(1:5,1)];
 disp('**********************************************');
 disp('IV Estimation by cost shifters in Logit Model:');
 disp('**********************************************');
 disp([' Vars','       Coeff','    ','Std Err']);
 disp([vars,num2str(str)]);
 disp(' ');
 disp(' ');
 disp(' ');
 
 disp('Endogeous bias is upward bias, suggestting that the unobserved') 
 disp('endogeous which may be postively correlatd with price and also ')
 disp('have postive effect on market share, e.g. advertising will inc-')
 disp('rease price but will also make products more attractive to consumers')
 
%% 1.2.2 Own-price elasticities /cross price elasticities 



epsilon1= bhat2SLS(2,1).*p.*(1-s);
epsilon2= -bhat2SLS(2,1).*p.*s;
array1= [ mid j p s epsilon1 epsilon2];




for i=1:max(j)
    epsilon_jj(i,1)= ((array1(:, 2) == i)'*epsilon1)/ sum(array1(:, 2) == i);
end


epsilon_matrix= zeros(max(j),max(j));
for r= 1:max(j)

    for c= 1:max(j)
        if r==c
            epsilon_matrix(r,c)= epsilon_jj(r,1);
       
        else 
            epsilon_matrix(r,c)= (array1(:, 2) == c)'*epsilon2 / sum(array1(:, 2) == c);
          
        end
    end
end  
 disp('**********************************************');
 disp('Slutsky Substitution Matrix is: ');
 disp('**********************************************');
 epsilon_matrix
 
 
%% 1.3 GMM 

 global alpha gamma0 gamma1 gamma2 gamma3 b_x1 b_x2 b_x3 b0 
 
 
 
 
 