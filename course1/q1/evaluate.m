function [xhat,q]= evaluate(Coef,xnode,a,b)
B=@splinefunc;
n=length(xnode);%index of xi is from 0 to n-1
h=xnode(2)-xnode(1);
distance=(b-a)/(20*n-20);%20*(n-1)+1 point with same distances

xhat=zeros(1,20*n-19);
q=zeros(1,20*n-19);
xhat(1)=xnode(1);
xhat(end)=xnode(end);
q(1)=Coef(1)+4*Coef(2)+Coef(3);%condition for q(x) in node point
q(end)=Coef(end-2)+4*Coef(end-1)+Coef(end);

for i=2:20*n-20
    xhat(i)=xhat(i-1)+distance;
    k=floor((xhat(i)-xnode(1))/h)+3;%location of xhat(i)
    if (k<4)
    %B-1(x)=B((x-x0+h)/h)
    q(i)=Coef(1)*B((xhat(i)-xnode(1)+h)/h)+Coef(2)*B((xhat(i)-xnode(1))/h)+Coef(3)*B((xhat(i)-xnode(2))/h)+Coef(4)*B((xhat(i)-xnode(3))/h);
    elseif(k>n)
    %B(n+1)(x)=B((x-xn-h)/h)
    q(i)=Coef(k-2)*B((xhat(i)-xnode(k-3))/h)+Coef(k-1)*B((xhat(i)-xnode(k-2))/h)+Coef(k)*B((xhat(i)-xnode(k-1))/h)+Coef(k+1)*B((xhat(i)-xnode(k-1)-h)/h);
    else
    q(i)=Coef(k-2)*B((xhat(i)-xnode(k-3))/h)+Coef(k-1)*B((xhat(i)-xnode(k-2))/h)+Coef(k)*B((xhat(i)-xnode(k-1))/h)+Coef(k+1)*B((xhat(i)-xnode(k))/h);
    end
end