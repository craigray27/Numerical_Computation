%Compute  the  LSQ  coefficients
function C=Lsqcoef(Func,n)
G=@Gauss;
C=G(Func,n);%Compute inner product(f(x),Pi(x))
for i=1:n+1
    C(i)=(2*i-1)*C(i)/2; %inner product(Pk(x),Pk(x)=2/(2k+1)
end