%F is vector-valued function that F(i)=f(x)*Pi(x) where Pi(x) is Legendre
function F=Func(x,n)
P=@Legendre;
f=@funct;
dimension=length(x);%length of vector x
F=zeros(n+1,dimension);
z=P(x,n);
for i=1:n+1
F(i,:)=f(x).*z(i,:);%F(i,:)=(Fi(x0),Fi(x1),...,Fi(xd))
end
end