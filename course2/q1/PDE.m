%Function to solve u and v
function [u,v]=PDE(L,N,M,h,deltat,betau,betav,k0,P)
   tridpro=@triproduct;
   %N:number of grids of x, M: number of grids of t
   %h,deltat: step size of x, t respectively
   n=N-2;
   u=zeros(M,N);
   v=zeros(M,N);
   u(1,:)=2.0;
   u(1,1:floor(0.9/h))=0.1;
   v(1,:)=1.97;
   
   mu=deltat/(2*h^2);
   
   A1u(1:n-1)=-betau*mu;
   A2u(2:n-1)=1+betau*mu*2;
   A2u(1)=1+betau*mu;
   A2u(n)=1+betau*mu;
   A3u(1:n-1)=-betau*mu;   
   B1u=-A1u;
   B2u(2:n-1)=1-betau*mu*2;
   B2u(1)=1-betau*mu;
   B2u(n)=1-betau*mu;
   B3u=-A3u;
   
   A1v(1:n-1)=-betav*mu;
   A2v(2:n-1)=1+betav*mu*2;
   A2v(1)=1+betav*mu;
   A2v(n)=1+betav*mu;
   A3v(1:n-1)=-betav*mu;   
   B1v=-A1v;
   B2v(2:n-1)=1-betav*mu*2;
   B2v(1)=1-betav*mu;
   B2v(n)=1-betav*mu;
   B3v=-A3v;
   
   A1uu=[A1u(1),A1u];
   A3uu=[A3u(1),A3u];
   Au=spdiags([A1uu' A2u' A3uu'],[-1 0 1],n,n);
   
   A1vv=[A1v(1),A1v];
   A3vv=[A3v(1),A3v];
   Av=spdiags([A1vv' A2v' A3vv'],[-1 0 1],n,n);

  for i=2:M
       fuv=v(i-1,2:end-1).*(k0+((u(i-1,2:end-1).^2)./(1+u(i-1,2:end-1).^2)))-u(i-1,2:end-1);
       du=tridpro(B1u,B2u,B3u,u(i-1,2:end-1))+deltat*fuv;
       u(i,2:end-1)=(Au\du')';
       
       dv=tridpro(B1v,B2v,B3v,v(i-1,2:end-1))-deltat*fuv;
       v(i,2:end-1)=(Av\dv')';
  end
   u(:,1)=u(:,2);
   u(:,end)=u(:,end-1);
   v(:,1)=v(:,2);
   v(:,end)=v(:,end-1);
end