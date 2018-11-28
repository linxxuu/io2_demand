function f=delta_fp(sigmavector, delta_old)
global s delta_save
diff=1000;
count=0;
while diff>1/1
delta_new=delta_old+log(s)-log(step1(sigmavector,delta_old));  
diff=(delta_new-delta_old)'*(delta_new-delta_old);
delta_old=delta_new;
count=count+1
end
delta_save=delta_new;
f=delta_new;
end