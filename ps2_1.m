% PS2 Demand

%% Input Data
clc
clear all
global s s0 x1 x2 x3 p c1 c2 c3 j mid

%M= csvread('Data.csv');
Data= dlmread('Data.csv');

mid = Data(:, 1);
j = Data(:, 2);
s = Data(:, 3);
x1 = Data(:, 4);
x2 = Data(:, 5);
x3 = Data(:, 6);
p = Data(:, 7);
c1= Data(:, 8);
c2= Data(:, 9);
c3= Data(:, 10);
group_id= Data(:, 11);


% calculate outside options: 
s0= zeros(970,1);
for i=1:970
    
    s0(i)= 1- (mid(:) == mid(i))'*s;
    
end

    
%% 1.1 Pooled OLS 

Y= log(s)-log(s0);
X= [ones(970,1) p x1 x2 x3];
    n = size(X,1);
    n_vars = size(X, 2);

bhat = (X'*X)\(X'*Y);
se1ols= sqrt(diag(n/(n-n_vars)*mean((Y-X*bhat).^2)*((X'*X)\eye(size(X,2)))));

vars = ['const     ';'Price     ';'X1        ';'X2        ';'X3        '];
varstr = [bhat(1:5,1) se1ols(1:5,1)];
disp('**********************************************');
disp('Pooled OLS Estimation in Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(varstr)]);
disp(' ');
disp(' ');
disp(' ');
 
%% 1.2.1 IV

% first stage: cost shifters and product characteristics 
% *Question*: whether including constant term in first stage? 
% Ante: Yes. First stage should generally include _all controls_ used in
% second stage (as we are interested in conditional correlation of the
% instrument and the variable of interest). Constant can be viewed as "one
% control".

Z= [ones(970,1) c1 c2 c3 x1 x2 x3];
    n = size(Z,1);
    n_vars = size(Z, 2);

bhatFS = (Z'*Z)\(Z'*p);
se1FS= sqrt(diag(n/(n-n_vars)*mean((p-Z*bhatFS).^2)*((Z'*Z)\eye(size(Z,2)))));

vars = ['const     ';'w1        ';'w2        ';'w3        ';...
    'X1        ';'X2        ';'X3        '];
str = [bhatFS se1FS];
disp('**********************************************');
disp('IV: First Stage Estimation in Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(str)]);
disp(' ');
disp(' ');
disp(' ');

 % NOTE: need to calculate F and p 
 
 
phat= Z*bhatFS;
X= [ones(970,1) p x1 x2 x3];
Y= log(s)-log(s0);
Xhat= [ones(970,1) phat x1 x2 x3];
    n = size(X,1);
    n_vars = size(X, 2);
     
bhat2SLS = (Xhat'*X)\(Xhat'*Y); 
%NOTE: these 2SLS SE are wrong.
se2SLS= sqrt(diag(mean(n/(n-n_vars)*(Y-X*bhat2SLS).^2)*...
    ((Xhat'*Xhat)\eye(size(X,2)))));
 
vars = ['Const     ';'Price     ';'X1        ';'X2        ';'X3        '];
str = [bhat2SLS se2SLS];
disp('**********************************************');
disp('IV Estimation by cost shifters in Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(str)]);
disp(' ');
disp(' ');
disp(' ');

disp('Endogeous bias is upward bias, suggestting that the unobserved') 
disp('endogeous which may be postively correlatd with price and also ')
disp('have postive effect on market share, e.g. advertising will inc-')
disp('rease price but will also make products more attractive to consumers')
 
%% 1.2.2 Own-price elasticities /cross price elasticities 

%NOTE: this should be done market by market, and only then averaged over
%all the markets to get the answer. I am not sure it gives the same answer.
%Ante

epsilon1= bhat2SLS(2,1).*p.*(1-s);
epsilon2= -bhat2SLS(2,1).*p.*s;
array1= [mid j p s epsilon1 epsilon2];

epsilon_jj = NaN([max(j) 1]);
for i=1:max(j)
    epsilon_jj(i,1)= ((array1(:, 2) == i)'*epsilon1)/ sum(array1(:, 2) == i);
end

epsilon_matrix= zeros(max(j),max(j));
for r= 1:max(j)
    for c= 1:max(j)
        if r==c
            epsilon_matrix(r,c)= epsilon_jj(r,1);
        else 
            epsilon_matrix(r,c)= (array1(:, 2) == c)'*epsilon2/...
                sum(array1(:, 2) == c);          
        end
    end
end  
disp('**********************************************');
disp('Slutsky Substitution Matrix is: ');
disp('**********************************************');
disp(epsilon_matrix);

%% 1.3. 2. GMM, version: Lin

start_theta = [bhat2SLS; bhatFS(1:4)];
options = optimset('Display','iter', 'TolX',1e-3, 'TolFun',1e-8,  'PlotFcns', @optimplotfval);
[param_gmm,fval,exitflag,output,grad1,hess1] = fminunc('gmm_obj', start_theta, options);
sig_param_gmm = gmm_obj_vcov(param_gmm);
%NOTE: SE seem fishy here. check function sig_param_gmm. Delta method used.
se_param_gmm = sqrt(diag(sig_param_gmm));
 
vars = ['b0        ';'alpha     ';'b1        ';'b2        ';'b3        ';...
    'gamma0    ';'gamma1    ';'gamma2    ';'gamma3    '];
disp('**********************************************');
disp('GMM estimation by cost shifters and supply equations:');
disp('**********************************************');

disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(param_gmm), blanks(9)', blanks(9)', blanks(9)', num2str(se_param_gmm)]);
  
  
%% 2 Nested logit
% 2.1 create sg contains group share: 
sh_gp= zeros(970,1);
for i=1:970
    sh_gp(i)= (group_id==group_id(i) & mid==mid(i))'*s; 
end
sh_g = s./sh_gp;

% 2.2. center all vars around product 1 (within group-market pair)
x1_ng = x1;
x2_ng = x2;
x3_ng = x3;
p_ng = p;
c1_ng = c1;
c2_ng = c2;
c3_ng = c3;
sh_g0 = sh_g;

for mid_i_=1:50
    for group_i_=1:3
        which_ = (group_id == group_i_) & (mid == mid_i_);
        zero_ng = find(which_, 1);
        %adjust _ng vars
        sh_g0(which_) = sh_g(which_)./sh_g(zero_ng);
        x1_ng(which_) = x1_ng(which_) - x1_ng(zero_ng);
        x2_ng(which_) = x2_ng(which_) - x2_ng(zero_ng);
        x3_ng(which_) = x3_ng(which_) - x3_ng(zero_ng);
        p_ng(which_) = p_ng(which_) - p_ng(zero_ng);
        c1_ng(which_) = c1_ng(which_) - c1_ng(zero_ng);
        c2_ng(which_) = c2_ng(which_) - c2_ng(zero_ng);
        c3_ng(which_) = c3_ng(which_) - c3_ng(zero_ng);
    end
end

% 2.3. within group-market: NestedLogit IV
    %1st stage
Z_ng= [ones(970,1) c1_ng c2_ng c3_ng x1_ng x2_ng x3_ng];
bhatFS_NL = (Z_ng'*Z_ng)\(Z_ng'*p_ng);

    %2nd stage
phat_NL= Z_ng*bhatFS_NL;
X_ng= [ones(970,1) p_ng x1_ng x2_ng x3_ng];
Y_ng= log(sh_g0);
Xhat_NL= [ones(970,1) phat_NL x1_ng x2_ng x3_ng];
 n = size(X_ng,1);
 n_vars = size(X_ng, 2);

bhat_1sNL = (Xhat_NL'*X_ng)\(Xhat_NL'*Y_ng); 
%NOTE: these 2SLS SE are wrong.
se_1sNL= sqrt(diag(mean(n/(n-n_vars)*(Y-X*bhat_1sNL).^2)*...
 ((Xhat_NL'*Xhat_NL)\eye(size(X_ng,2)))));

vars = ['Const     ';'Price     ';'X1        ';'X2        ';'X3        '];
str = [bhat_1sNL se_1sNL];
disp('**********************************************');
disp('IV Estimation by cost shifters in Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(str)]);
disp(' ');
disp(' ');
disp(' ');

 % 2.4 between group-market: Logit
data_g = zeros(150,5);
for mid_i_=1:50
for group_i_=1:3
    row_i_ = 3*mid_i_+group_i_-3;
    data_g(row_i_,1)= mid_i_;
    data_g(row_i_,2)= group_i_;

    which_ = (group_id == group_i_) & (mid == mid_i_);
    s_ng = (which_)'*s;
    Ign_tilde = log(sum(exp(X_ng(which_,:)*bhat_1sNL)));
    s0_ng= 1- (mid(:) == mid_i_)'*s;

    data_g(row_i_,3)= s_ng;
    data_g(row_i_,4)= Ign_tilde;
    data_g(row_i_,5)= s0_ng;
end
end
 
Y_g= log(data_g(:,3))-log(data_g(:,5));
X_g = [ones(size(data_g(:,4))) data_g(:,4)];
bhat_2sNL = (X_g'*X_g)\(X_g'*Y_g); 
sigma_NL = 1- bhat_2sNL(2);

bhat_NL = bhat_1sNL*(1-sigma_NL);

% 2.5. Bootstrapping S.E.
N_runs = 100;
cofs_bstrap = zeros(length(bhat_NL), N_runs);

for run_i_ = 1:N_runs
    mids_ = randsample(50, 25);
    subsample_i = ismember(mid,mids_);
    subsample_ng = ismember(data_g(:,1),mids_);
    %Within group-market: 1st stage
    n_ss = length(find(subsample_i));
    Z_ng_bs= [ones(n_ss,1) c1_ng(subsample_i) c2_ng(subsample_i) c3_ng(subsample_i)...
        x1_ng(subsample_i) x2_ng(subsample_i) x3_ng(subsample_i)];
    p_ng_bs = p_ng(subsample_i);

    bhatFS_NL_bs = (Z_ng_bs'*Z_ng_bs)\(Z_ng_bs'*p_ng_bs);
    %Within group-market: 2nd stage
    phat_NL_bs= Z_ng_bs*bhatFS_NL_bs;
    X_ng_bs= [ones(n_ss,1) p_ng_bs x1_ng(subsample_i) x2_ng(subsample_i) x3_ng(subsample_i)];
    Y_ng_bs= log(sh_g0(subsample_i));
    Xhat_NL_bs= [ones(n_ss,1) phat_NL_bs x1_ng(subsample_i) x2_ng(subsample_i) x3_ng(subsample_i)];
    
    bhat_1sNL_bs = (Xhat_NL_bs'*X_ng_bs)\(Xhat_NL_bs'*Y_ng_bs);    
    %Between groups- within market
    Y_g_bs= log(data_g(subsample_ng,3))-log(data_g(subsample_ng,5));
    X_g_bs = [ones(size(data_g(subsample_ng,4))) data_g(subsample_ng,4)];
    bhat_2sNL_bs = (X_g_bs'*X_g_bs)\(X_g_bs'*Y_g_bs); 
    sigma_NL_bs = 1- bhat_2sNL_bs(2);

    bhat_NL_bs = bhat_1sNL_bs*(1-sigma_NL_bs);

    cofs_bstrap(:,run_i_) = bhat_NL_bs;
end

se_NL = std(cofs_bstrap,0,2);

vars = ['Const     ';'Price     ';'X1        ';'X2        ';'X3        '];
str = [bhat_NL se_NL];
disp('**********************************************');
disp('IV Estimation by cost shifters in Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(str)]);
disp(' ');
disp(' ');
disp(' ');

 % 3 Random coefficients logit model
 %% 3.1 
 
 % NS: number of simulations (each simulation is similar to each consumer) 
 NS= 20;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 