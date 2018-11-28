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

