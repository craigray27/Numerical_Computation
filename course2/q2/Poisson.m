%Function to solve Poisoon PDE
function [u,U,Freal]=Poisson(N,fx,gx)
%N+1 is the grids for both x and y, 1/N is step size
%u is the result stored in vector form(u11,u12,u13....u1N,....uNN)
%U is the result stored in Matrix form
%Freal is analytical solution stored in vector form
 h=1/N;
 m=N-1;
 nm=2*((m-1)*(3*m-1)+2*m-1)-m^2;
 row=zeros(nm,1);
 col=zeros(nm,1);
 data=zeros(nm,1);
 F=zeros(m*m,1);
 G=zeros(m*m,1);
 Freal=zeros(m*m,1);
 
 x=h;
 y=h;
 q=1;
 row(1:3)=1;
 col(1)=q;
 col(2)=q+1;
 col(3)=q+m;
 data(1)=4;
 data(2)=-1;
 data(3)=-1;
 
 F(1)=fx(x,y);
 G(1)=gx(x)+gx(y);
 Freal(1)=cos(4*pi*x)*cos(4*pi*y);
 k=2;
 q=2;
 i=4;
 
 while(k<=m*m)
     if mod(k-1,m)~=0
         x=x+h;
         F(k)=fx(x,y); 
         Freal(k)=cos(4*pi*x)*cos(4*pi*y);
     else
         x=h;
         y=y+h;
         F(k)=fx(x,y); 
         Freal(k)=cos(4*pi*x)*cos(4*pi*y);
     end
     
     if (k<=m) && (mod(k,m)~=0)
         row(i:i+3)=k;
         col(i:i+2)=q-1:1:q+1;
         col(i+3)=q+m;
         data(i)=-1;
         data(i+1)=4;
         data(i+2)=-1;
         data(i+3)=-1;
         G(k)=gx(x);
         q=q+1;
         k=k+1;
         i=i+4;
       
     elseif (k<=m) && (mod(k,m)==0)
         row(i:i+2)=k;
         col(i:i+1)=q-1:1:q;
         col(i+2)=q+m;
         data(i)=-1;
         data(i+1)=4;
         data(i+2)=-1;
         G(k)=gx(x)+g(y);
         q=q+1;
         k=k+1;
         i=i+3;
         
      elseif (k>m) && (mod(k-1,m)==0) && (k<=m*m-m)
          row(i:i+3)=k;
          col(i)=q-m;
          col(i+1:i+2)=q:1:q+1;
          col(i+3)=q+m;
          data(i)=-1;
          data(i+1)=4;
          data(i+2)=-1;
          data(i+3)=-1;
          q=q+1;

          i=i+4;
          G(k)=g(y);
          k=k+1;
          
       elseif (k>m) && (mod(k-1,m)~=0) && (mod(k,m)~=0) &&(k<=m*m-m) 
          row(i:i+4)=k;
          col(i)=q-m;
          col(i+1:i+3)=q-1:1:q+1;
          col(i+4)=q+m;
          data(i)=-1;
          data(i+1)=-1;
          data(i+2)=4;
          data(i+3)=-1;
          data(i+4)=-1;
          q=q+1;

          i=i+5;
          G(k)=0;
          k=k+1;
       elseif (k>m) && (mod(k,m)==0) &&(k<=m*m-m)
          row(i:i+3)=k;
          col(i)=q-m;
          col(i+1:i+2)=q-1:1:q;
          col(i+3)=q+m;
          data(i)=-1;
          data(i+1)=-1;
          data(i+2)=4;
          data(i+3)=-1;
          q=q+1;

          i=i+4;
          G(k)=g(y);
          k=k+1;
          
        elseif (mod(k-1,m)==0) &&(k>m*m-m)
            row(i:i+2)=k;
            col(i)=q-m;
            col(i+1:i+2)=q:1:q+1;
            data(i)=-1;
            data(i+1)=4;
            data(i+2)=-1;
            q=q+1;
     
            i=i+3;
            G(k)=g(x)+g(y);
            k=k+1;
            
        elseif (mod(k,m)~=0) &&(k>m*m-m)
            row(i:i+3)=k;
            col(i)=q-m;
            col(i+1:i+3)=q-1:1:q+1;
            data(i)=-1;
            data(i+1)=-1;
            data(i+2)=4;
            data(i+3)=-1;
            q=q+1;
  
            i=i+4;
            G(k)=g(x);
            k=k+1;
         elseif (mod(k,m)==0) &&(k>m*m-m)
            row(i:i+2)=k;
            col(i)=q-m;
            col(i+1:i+2)=q-1:1:q;
            data(i)=-1;
            data(i+1)=-1;
            data(i+2)=4;
            q=q+1;
         
            i=i+3;
            G(k)=g(x)+g(y);
            k=k+1;
     end
     
     
 end
 
 A=sparse(row,col,data,m*m,m*m);%sparse matrix
 
 u=A\(h^2*F+G);
 U=zeros(m,m);
 for j=1:m
     U(j,:)=u((j-1)*m+1:j*m);
 end
end