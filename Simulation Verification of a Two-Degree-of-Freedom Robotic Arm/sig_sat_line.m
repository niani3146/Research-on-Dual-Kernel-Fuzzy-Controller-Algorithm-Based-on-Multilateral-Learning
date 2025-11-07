function y = sig_sat_line(x,bound)
[m,n]=size(x);
if min([m,n])>=2
    disp('wrong')
elseif max([m,n])==1
    y=2*bound/(1+exp(-x))-bound;
else
    y=2*bound./(1+exp(-x))-bound;
end