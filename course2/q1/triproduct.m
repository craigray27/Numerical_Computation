%Function to compute A*x where A is a tridiagnal matrix, A=(abc)
function y=triproduct(a,b,c,x)
 n=length(x);
 y(1)=b(1)*x(1)+c(1)*x(2);
 y(2:n-1)=b(2:n-1).*x(2:n-1)+c(2:n-1).*x(3:n)+a(1:n-2).*x(1:n-2);
 y(n)=a(n-1)*x(n-1)+b(n)*x(n);
end