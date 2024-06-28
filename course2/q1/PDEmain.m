L=1;M=1000;N=100;
h=1/99;deltat=200/999;
P=2.26;
betav=0.1;betau=0.01*betav;
n=N-2;
k0=0.067;
pde=@PDE;
[u,v]=pde(L,N,M,h,deltat,betau,betav,k0,P);
   
   %final time plot, delete'%' if you wanna see
   %plot(0:h:1,v(end,:))
   %ylabel('v')
   %xlabel('x')
   
   %contour plot, delete'%' if you wanna see
   %x1=0:h:1;
   %t1=0:deltat:200;
   %[X,Y]=meshgrid(x1,t1);
   %contourf(X,Y,v);
   %ylabel('t')
   %xlabel('x')
   
 hj=[h,h/2,h/4,h/8];
e=zeros(3,1);
[u2,~]=pde(L,199,M,hj(2),deltat,betau,betav,k0,P); % h1/2 
[u3,~]=pde(L,397,M,hj(3),deltat,betau,betav,k0,P); % h2/2 
[u4,~]=pde(L,793,M,hj(4),deltat,betau,betav,k0,P); % h3/2 
e(1)=abs(hj(2)*sum(u2(:))-hj(1)*sum(u(:)));
e(2)=abs(hj(3)*sum(u3(:))-hj(2)*sum(u2(:)));
e(3)=abs(hj(4)*sum(u4(:))-hj(3)*sum(u3(:)));

lge=log(e);lgh=log(hj(1:3));
plot(log(hj(1:3)),log(e))
hold on
scatter(log(hj(1:3)),log(e))

   
   
   
 
   
   