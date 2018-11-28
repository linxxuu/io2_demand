for i=1:50
    for k=1:3
        which = (group_id == k) & (mid == i);
        tab = tabulate(j(which));
        if tab(1) ~= 1
            disp(tab);
        end
    end
end


disp([vars, num2str([bhat2SLS; bhatFS(1:4)]), blanks(9)', blanks(9)', blanks(9)',...
    num2str(param_gmm), blanks(9)', blanks(9)', blanks(9)',...
    num2str(param_gmm_ng), blanks(9)', blanks(9)', blanks(9)']);

%%
% % 2.3. within group-market: NestedLogit GMM
% start_theta = param_gmm; %[bhat2SLS; bhatFS(1:4)]; %param_gmm;
% options = optimset('Display','iter', 'TolX',1e-3, 'TolFun',1e-8,  'PlotFcns', @optimplotfval);
% [param_gmm_ng,fval,exitflag,output,grad1,hess1] = fminunc('gmm_ng_obj', start_theta, options);
% 
% vars = ['b0        ';'alpha     ';'b1        ';'b2        ';'b3        ';...
%  'gamma0    ';'gamma1    ';'gamma2    ';'gamma3    '];
% disp('**********************************************');
% disp('(Nested Logit) GMM estimation by cost shifters and supply equations:');
% disp('**********************************************');
% 
% disp([' Vars','       Coeff','    ','Std Err']);
% disp([vars, num2str(param_gmm), blanks(9)', blanks(9)', blanks(9)']);


%% 1.3.1 GMM, version: Ante
Y= log(s)-log(s0);
X= [ones(970,1) p x1 x2 x3];
c = [c1 c2 c3];

%suppy side residual
supp_e0 = @(coef) (p + 1./(coef(2).*(1-s)) ...
    - c*coef(3:5)' - coef(1)*ones(970,1));
%demand side residual
dem_e0 = @(coef) (Y-(X*[coef]'));

res_gmm13 = @(coef) ...
    (sum((c'*supp_e0(coef([6 2 7:9]))).^2) + ...
    sum((c'*dem_e0(coef(1:5))).^2));

%supply side starting values
Y_13 = p + 1./(bhat2SLS(2).*(1-s));
X_13 = [ones(970,1) c1 c2 c3];
cg_start13 = (X_13'*X_13)\(X_13'*Y_13);

start_gmm13 = [bhat2SLS; cg_start13]';
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);
bhat_gmm13 = fminsearch(res_gmm13, start_gmm13, options)';

%find SE of GMM
se_gmm13


vars = ['Const     ';'Price     ';'X1        ';'X2        ';...
    'X3        ';'C1        ';'C2        ';'C3        '];
str = [bhat_gmm13 se_gmm13];
disp('**********************************************');
disp('GMM: Demand+Supply (only cost shifters) Logit Model:');
disp('**********************************************');
disp([' Vars','       Coeff','    ','Std Err']);
disp([vars, num2str(str)]);
disp(' ');
disp(' ');
disp(' ');
 