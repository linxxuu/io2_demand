% rcl file

%% Input Data
clc
clear all
cd 'C:\Users\jschaerer\Documents\emp.IO 2018\ps2\Julian RCL';
%M= csvread('Data.csv');
Data= dlmread('Data.csv');

global s s0 x1 x2 x3 p c1 c2 c3 j mid X1 X2 Z group_id
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
group_id= Data(:, 11);

X2=[p x1 x2 x3]
X1=[ones(970,1) X2];
Z=[ones(970,1) x1.^2 x2.^2 x3.^2 c1 c1.^2 c2 c2.^2 c3 c3.^2]; 

% calculate outside options: 
s0= zeros(970,1);
for i=1:970
    s0(i)= 1- (mid(:) == mid(i))'*s;
end


% parameters
global theta1 delta_save alpha gamma0 gamma1 gamma2 gamma3 b_x1 b_x2 b_x3 b_0 sigma0 sigma1 sigma2 sigma3 sigma

%starting values
alpha=1;
gamma0=1;
gamma1=1;
gamma2=1;
gamma3=1;
b_x1=1;
b_x2=1;
b_x3=1;
b_0=1;

sigma0=0.1;
sigma1=0.2;
sigma2=0.3;
sigma3=0.4;
sigma=[sigma0 sigma1 sigma2 sigma3]';
% draw random numbers
global vfull ns mu delta
ns=100;
%randn('seed',1);
vfull=randn([1 4*ns]);

vfull=repmat(vfull, 970,1);
delta=log(s)-log(s0);
%delta_save=delta;
%delta_save not needed, after all

%%
theta2= [0.586 0.936 0.45 0.28]';
[param,fval,exitflag,output,grad,hess] = fminunc('gmm_rcl', theta2, optimset('Display','iter','PlotFcn','optimplotfval','MaxFunEvals',10000,'TolX',1e-6,'TolFun',1e-8));
 vars = ['b0        ';'alpha     ';'b1        ';'b2        ';'b3        '];
%%REST: can be ignored (old stuff)

%delta=delta_fp(sigma, delta);
%[param,fval,exitflag,output,grad]=fminunc(@gmm_rcl, [-1.17 2.88 1.11 .75]')

%theta2 = fminunc('gmm_rcl',theta2);

% recommendation: run one of the following two. (2nd is nicer because it
% plots the obj. function value)
%theta2= [0.586 0.936 0.45 0.28]';
%[param,fval,exitflag,output,grad,hess]=fminunc('gmm_rcl', theta2)


 % or different plotting option:'PlotFcn','optimplotfval',
 
 
