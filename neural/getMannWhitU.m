function [U1, U2, N] = getMannWhitU(group1, group2)

%Compare the values between group 1 and group 2 to determine how many times
%a value from group1 is greater than one from group2 after all possible
%comparisons. group1 and group2 are column vectors. U1 is the U-statistic
%for group1, U2 for group2, N is the total number of comparisons (length of
%group1 x length of group 2)

ng1 = numel(group1);
ng2 = numel(group2);
N = ng1*ng2;
t = tiedrank([group1(:); group2(:)]); 
R1 = sum(t(1:ng1)); 
R2 = sum(t(ng1+1:end));
U1 = R1 - ng1*(ng1+1)/2;
U2 = N - U1;
