function J=Jaccobi(Func,x,epsilon)
n=length(x);
J=zeros(n,n);
deltax=zeros(1,n); %delatax=x+epsilon*ej
for i=1:n
    deltax(i)=epsilon;
    J(:,i)=((Func(x+deltax)-Func(x))/epsilon)';
    deltax(i)=0;
end
end