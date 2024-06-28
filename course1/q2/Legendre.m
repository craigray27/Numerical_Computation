%Vector-valued function that P(i)=Pi(x) where Pi(x) is Legendre
function P=Legendre(x,n)
 dimension=length(x);
 P=zeros(n+1,dimension);
 P(1,:)=1;
 P(2,:)=x;
 for i=3:n+1
     %P(i,:)=(Pi(x0),Pi(x1),...,Pi(xd))
     P(i,:)=(2*i-3)*x.*P(i-1,:)/(i-1)-(i-2)*P(i-2,:)/(i-1);
 end
end