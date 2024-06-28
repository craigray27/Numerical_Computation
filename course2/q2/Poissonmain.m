fx=@f;
gx=@g;
P=@Poisson;
N=64;
h=1/N;
[u,U,Freal]=P(N,fx,gx);

%contour plot, delete'%' if you wanna see
   %x1=h:h:1-h;
   %y1=h:h:1-h;
   %[X,Y]=meshgrid(x1,y1);
   %contourf(X,Y,U);
   %ylabel('y')
   %xlabel('x')
   %zlabel('u(x,y)')
 Nj=[16,32,64,128,256,512,1024];
 error=zeros(length(Nj),1);
 error2=zeros(length(Nj),1);
 for i=1:length(Nj)
   [u,~,Freal]=P(Nj(i),fx,gx);
   error(i)=sqrt(1/Nj(i)*sum((Freal-u).^2));
   error2(i)=max(abs((Freal-u)));
 end
 
 lge=log(error);
 lgN=log(Nj);
 plot(log(Nj),lge)
 hold on
 scatter(log(Nj),lge)
 xlabel('log(Nj)')
 ylabel('log(L-2 norm)')