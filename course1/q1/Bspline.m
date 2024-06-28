%Find coefficients 
function C= Bspline(f)
 n=length(f);
 trid=@tridisolve;
 a=zeros(1,n-3);
 b=zeros(1,n-2);
 c=zeros(1,n-3);
 d=zeros(1,n-2);
 a(1:n-3)=1;
 b(1:n-2)=4;
 c(1:n-3)=1;
 
 d(1)=f(2)-f(1)/6;
 d(2:n-3)=f(3:n-2);
 d(n-2)=f(n-1)-f(n)/6;
 
 C=zeros(1,n+2);
 C(3:n)=trid(a,b,c,d);
 C(2)=f(1)/6;
 C(1)=2*C(2)-C(3);%first-derivative conditions
 C(n+1)=f(n)/6;
 C(n+2)=2*C(n+1)-C(n);
 end