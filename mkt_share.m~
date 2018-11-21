function f = mkt_share(thet,delta)
% calculate market shares if we know parameters (theta1= alpha, beta0~3)
% and (theta2= (sigma0~sigma3))
% delta input is 970*1
global ns v X1 X2 mid nc cid den

clear mu e_delta expu term1 ind_share mrk term_0  term_1 term_2 term_3 term_4


mu= repmat(thet,[nc,1]).*v*(X2');

e_delta= kron(delta', ones(nc,1));

% mu+e_delta (size is 5000*970) is for 5000 consumers (100 consumers in each market), their
% utilities for each products (jn) in each market
expu= exp(mu+e_delta);

term_0= den.*expu;

term_1= sum(term_0,2);

term_2 =1./(1+term_1);

term_3= repmat(term_2,1,970);

term_4= den.*term_3;

ind_share= expu.*term_4;

mrk= sum(ind_share,1)./ns;
% calculate denomenator for the 5000 consumers* 970 products matrix.
% for i=1:970
%  denomentator =  (mid==mid(i))'.* expu() 
%  expu(:,i)'*((mid==mid(i))'.*(1.(1+(mid==mid(i))'*expu(:,i)))');
% end

% calculate for each product jn in market n, market share= exp(ujn)/(1+sum exp(ukn))
% for i=1:970
%     
% mks(i)= expu(:,i)'*((mid==mid(i))'.*(1.(1+(mid==mid(i))'*expu(:,i)))');
% end
% 
% f=mks; 
f= mrk';
clear mu e_delta expu ind_share mrk delta thet term_0  term_1 term_2 term_3 term_4

end