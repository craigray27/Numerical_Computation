%Vector-valued function to do Gauss integrate
function F=testguassfunc(x,n)
dimension=length(x);
F=zeros(n,dimension);
if n==1
F(1,:)=exp(-x.^2);
end