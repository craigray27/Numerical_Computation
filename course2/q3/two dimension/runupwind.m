a=[1,0.5];
deltax=2*pi/20;
 deltat=0.1;
 N=6*pi/deltax+1;
 M=1/deltat+1;
 uwind=zeros(N,N,M);%Upwind
 ureal=zeros(N,N,M);%Analytical solution
 x0=-2*pi:deltax:4*pi;
 y0=-2*pi:deltax:4*pi;
 uwind(:,:,1)=(1-cos(x0))'*(1-cos(y0));
 ureal(:,:,1)=(1-cos(x0))'*(1-cos(y0));
 lambda=deltat/deltax;
 
 for i=2:M
     for j=i:1:N-i+1
        for k=i:1:N-i+1
            lax1=-a(1)*lambda/2*(uwind(k+1,j,i-1)-uwind(k-1,j,i-1))+abs(a(1))*lambda/2*(uwind(k+1,j,i-1)-2*uwind(k,j,i-1)+uwind(k-1,j,i-1));
            lax2=-a(2)*lambda/2*(uwind(k,j+1,i-1)-uwind(k,j-1,i-1))+abs(a(2))*lambda/2*(uwind(k,j+1,i-1)-2*uwind(k,j,i-1)+uwind(k,j-1,i-1));
            uwind(k,j,i)=uwind(k,j,i-1)+lax1+lax2;
            
            ureal(k,j,i)=(1-cos(-2*pi+(k-1)*deltax-a(1)*i*deltat))*(1-cos(-2*pi+(j-1)*deltax-a(2)*i*deltat));
            
        end
     end
     
 end
 uwind=uwind(20:40,20:40,:);
 ureal=ureal(20:40,20:40,:);
 
 L2=zeros(11,1);
 k=1;
 for t=0:deltat:1
     L2(k)=deltax*norm(uwind(:,:,k),2);
     k=k+1;
 end
 
 result=uwind(:,:,end);
 real=ureal(:,:,end);
 
 x1=0:deltax:2*pi;
 y1=0:deltax:2*pi;
 [X,Y]=meshgrid(x1,y1);
 mesh(X,Y,real)
 xlabel('x')
 ylabel('y')
 zlabel('u(x,y,1)')
 title('t=1')
 
