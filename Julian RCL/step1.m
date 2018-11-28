function f=step1(sigmavector,deltavector)
%returns the shares given delta (and sigma); along the way, it uses the
%mufunc function.
global vfull ns X2 mid
%sigma0=sigma(1);
%sigma1=sigma(2);
%sigma2=sigma(3);
%sigma3=sigma(4);


mu_ijt=mufunc(X2,sigmavector);

delta_long=repmat(deltavector,1,ns);
exp_delta_mu=exp(delta_long+mu_ijt);
sum_exp_delta_mu=zeros(970,ns);
for samplenumb=1:ns
for i=1:970
      sum_exp_delta_mu(i,samplenumb)= ((mid(:) == mid(i))'*(exp_delta_mu(:,samplenumb)))+1;
end
end

%use share function
shares_j=sum(exp_delta_mu./sum_exp_delta_mu,2)/ns;
f=shares_j;
end