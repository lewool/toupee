%% set up

x  = prevChoiceValues(2:end)';
y  = allChoiceValues(2:end)';
z = prevContrastValues(2:end)';

X = [ones(length(x),1) x];
Z = [ones(length(z),1) z];
%% regular old OLS

b1 = X\y; 
yy1 = X*b1;
u1 = y - yy1;

b2 = Z\y; 
yy2 = Z*b1;
u2 = y - yy2;

corr(x,u1)


%% TSLS

d = Z\x; 
x_hat = Z*d;
err = x - x_hat;

b = [ones(length(x_hat),1) x_hat]\y; 
y_hat = [ones(length(x_hat),1) x_hat]*b;
noise = y - y_hat;