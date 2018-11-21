function f = cverg_delta(theta2w)

global s delta0

old_delta= delta0;
new_delta=zeros(970,1);

while abs(old_delta-new_delta)> 10e-3
    log_s_old= log(mkt_share(theta2w,old_delta));
    new_delta2=old_delta+log(s)-log_s_old;
    old_delta=new_delta;
    new_delta=new_delta2;
end

f=new_delta;
clear new_delta old_delta new_delta2


end