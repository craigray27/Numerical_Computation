function X=Newton(n,Func,x0,epsilon)
 J=@Jaccobi;
 k=length(x0);
 Xs=zeros(n,k);
 Xs(1,:)=x0;
 for i=2:n
     Xs(i,:)=Xs(i-1,:)+(J(Func,Xs(i-1,:),epsilon)\(-Func(Xs(i-1,:))))'; 
 end
 X=Xs(end,:);
end