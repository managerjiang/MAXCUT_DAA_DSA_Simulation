function y=f(x)
global A;
global bi;
global yv;
global sigma;
global C;
global n;

sumy=0;
sumqua=0;

for i=1:n
    sumy=sumy+yv(i,1)*(trace(A{i,1}'*(x*x'))-bi);
    sumqua=sumqua+(trace(A{i,1}'*(x*x'))-bi)^2;
end
y=trace(C'*(x*x'))-sumy+sigma/2*sumqua;
end