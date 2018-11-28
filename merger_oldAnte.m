chg_= 10;
p_old = p;
sh_old = s;
i=1;
clear p_hist s_hist
top2_0 = zeros(size(top2_));
while chg_ > 0.1
    Data_temp = Data_sorted;
    Data_temp(:, 7) = p_old;
    Data_temp(:, 3) = sh_old;
    [sh_new, p_new] = update_ps_logit(Data_temp, top2_, param_gmm);
    chg_ = max(abs(sh_new-sh_old)) + max(abs(p_new - p_old));
    p_old = p_new;
    sh_old = sh_new;
    disp(chg_);
    p_hist(:,i) = p_new;
    s_hist(:,i) = sh_new;
    Data_temp(:, 7) = p_new;
    Data_temp(:, 3) = sh_new;
    
    
    if i>1
        p_i = find(p_hist(:,i)-p_hist(:,i-1) == max(p_hist(:,i)-p_hist(:,i-1)));
        s_i = find(s_hist(:,i)-s_hist(:,i-1) == max(s_hist(:,i)-s_hist(:,i-1)));
        disp('p_i= ');
        disp(p_i);
        disp(max(p_hist(:,i)-p_hist(:,i-1)));
        disp('s_i= ');
        disp(s_i);
        disp(max(s_hist(:,i)-s_hist(:,i-1)));
    end
    i=i+1;
end