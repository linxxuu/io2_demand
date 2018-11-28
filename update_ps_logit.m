function [sh_new, p_new] = update_ps_logit(Data_temp, top2_, param)
    
    sh_new = Data_temp(:, 3);
    x1 = Data_temp(:, 4);
    x2 = Data_temp(:, 5);
    x3 = Data_temp(:, 6);
    p_new = Data_temp(:, 7);
    c1= Data_temp(:, 8);
    c2= Data_temp(:, 9);
    c3= Data_temp(:, 10);
    
    alpha = param(2);
    gamma = param(6:9);
    u_cof = param(1:5);

    %update prices, given market shares
    for mid_i_=1:50
        market_which_ = (Data_temp(:,1) == mid_i_);
        market_i_ = Data_temp(market_which_,:);
        
        start_market = find(market_which_,1);
        %update prices, given market shares
        for prod_j_ = start_market:(start_market+size(market_i_,1)-1)
            if top2_(prod_j_) == 1
                mc_1 = [1 c1(prod_j_) c2(prod_j_) c3(prod_j_)]*gamma;
                mc_2 = [1 c1(prod_j_+1) c2(prod_j_+1) c3(prod_j_+1)]*gamma;
                p_new(prod_j_) = p_merg(sh_new(prod_j_), mc_1,...
                    p_new(prod_j_+1), sh_new(prod_j_+1), mc_2, alpha);
            elseif top2_(prod_j_) == 2
                mc_1 = [1 c1(prod_j_) c2(prod_j_) c3(prod_j_)]*gamma;
                mc_2 = [1 c1(prod_j_-1) c2(prod_j_-1) c3(prod_j_-1)]*gamma; 
                p_new(prod_j_) = p_merg(sh_new(prod_j_), mc_1,...
                    p_new(prod_j_-1), sh_new(prod_j_-1), mc_2, alpha);      
            else
                mc_1 = [1 c1(prod_j_) c2(prod_j_) c3(prod_j_)]*gamma; 
                p_new(prod_j_) = p_merg(sh_new(prod_j_), mc_1,...
                    0, 0, 0, alpha);
            end
            
            %update marketshares, given new prices
            
            market_data = [ones(size(x1(market_which_))) p_new(market_which_)...
                x1(market_which_) x2(market_which_) x3(market_which_)];
            total_market_ = sum(exp(market_data*u_cof))+1;
            sh_new(start_market:(start_market+size(market_i_,1)-1)) =...
                exp(market_data*u_cof)./total_market_;
            
        end 
    end
    

end