function p_new = p_merg(sh_1, mc_1, p_2, sh_2, mc_2, alpha)
    
    p_new = max(mc_1 + sh_2*(p_2-mc_2)/(1-sh_1) - 1/(alpha*(1-sh_1)),0);

end