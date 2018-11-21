%% 1. Parameters 
clc
clear all
global s s0 x1 x2 x3 p c1 c2 c3 j mid X1 X2 ns  v Z delta0


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


X1=[p ones(970,1) x1 x2 x3]; 

X2=[p x1 x2 x3];


Z=[ones(970,1) c1 c2 c3 c1.^2 c2.^2 c3.^2 x1.^2 x2.^2 x3.^2]; 

s0= zeros(970,1);
for i=1:970
    
    s0(i)= 1- (mid(:) == mid(i))'*s;
    
end

%% 0. Prepare random coefficients draws 
ns=100;

% so for each market (mid), we need to draw #ns consumers, each contains
% (v0i, v1i,v2i,v3i,v4i) (dimension k=0,1,2,3)
v=[];
cid=[];
for i=1:max(mid)
v1= normrnd(0,1,[ns,4]);
cid1= i*ones(ns,1); 
v=[v;v1];
cid=[cid;cid1];
clear v1 cid1
end

global nc den

nc= ns*max(mid);
%% 1. 


delta0= log(s)-log(s0);

den= zeros(nc,970);
for i=1:970
    
den(:,i)=(cid==mid(i));


end

clear i



 theta2= [0.1 0.2 0.3 0.4];
 theta2 = fminunc('gmm_obj3',theta2);