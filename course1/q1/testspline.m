F=@func;
B=@splinefunc;
Bs=@Bspline;
e=@evaluate;
b=1;
a=-1;
n=[16,32,64,128];
err=zeros(4,1);
k=1;
for n1=n(1:end)
xnode=zeros(1,n1+1);
xnode(1)=a;
for i=2:n1+1
xnode(i)=xnode(i-1)+(b-a)/n1;
end

fnode=F(xnode);
Coef=Bs(fnode);

[xhat,q]=e(Coef,xnode,a,b);

err(k)=max(abs(F(xhat)-q));%||f-q||
k=k+1;
end

n1=128;
xnode=zeros(1,n1+1);
xnode(1)=a;
for i=2:n1+1
xnode(i)=xnode(i-1)+(b-a)/n1;
end
plot(log(n),log(err));
hold on
scatter(log(n),log(err));
legend('Estimated log(n) and log(error)','n=16k(k=1,2,3,4)')
xlabel('log(n)')
ylabel('log(error)')

fnode=F(xnode);
Coef=Bs(fnode);

[xhat,q]=e(Coef,xnode,a,b);
x=-1:0.01:1;
y=F(x);
plot(x,y)
hold on
plot(xhat,q)
hold on
scatter(xnode,fnode)
legend('f(x)','q(x)','node points')
xlabel('x')
ylabel('y=f(x)')
title('n=128')