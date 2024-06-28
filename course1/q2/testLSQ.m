f=@funct;
F=@Func;
C=@Lsqcoef;
P=@Legendre;
n=7;
Coef=C(F,n);%LSQ coefficients  
x=-1:0.01:1;
freal=funct(x);%true value
LSQ=Coef'*P(x,n);
plot(x,freal)
hold on
plot(x,LSQ)
legend('f(x)','LSQ(x)')
xlabel('x')
ylabel('y=f(x)')


