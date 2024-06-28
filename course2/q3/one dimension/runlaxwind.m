 deltax=2*pi/20;
 deltat=0.1;
 a=1;
 N=6*pi/deltax+1;
 M=1/deltat+1;
 ulax=zeros(M,N);%Lax-Wendro?
 uwind=zeros(M,N);%Upwind
 ureal=zeros(M,N);%Analytical solution
 x0=-2*pi:deltax:4*pi;
 ulax(1,:)=1-cos(x0);
 uwind(1,:)=1-cos(x0);
 ureal(1,:)=1-cos(x0);
 lambda=deltat/deltax;
 
 for i=2:M
     for j=i:1:N-i+1
        ulax(i,j)=ulax(i-1,j)-a*lambda/2*(ulax(i-1,j+1)-ulax(i-1,j-1))+a^2*lambda^2/2*(ulax(i-1,j+1)-2*ulax(i-1,j)+ulax(i-1,j-1));
        
        uwind(i,j)=uwind(i-1,j)-a*lambda/2*(uwind(i-1,j+1)-uwind(i-1,j-1))+abs(a)*lambda/2*(uwind(i-1,j+1)-2*uwind(i-1,j)+uwind(i-1,j-1));
        
        ureal(i,j)=1-cos(-2*pi+(j-1)*deltax-a*i*deltat);
     end
 end
 ulax=ulax(:,20:40);%select data where x=0:deltax:2*pi
 uwind=uwind(:,20:40);
 ureal=ureal(:,20:40);
 
 k=1;
 L2lax=zeros(11,1);
 L2wind=zeros(11,1);
 for t=0:deltat:1
 L2lax(k)=sum(deltax*(ulax(k,:)).^2);
 L2wind(k)=sum(deltax*(uwind(k,:)).^2);
 k=k+1;
 end
 
 x1=0:deltax:2*pi;
 t1=0:deltat:1;
 [X,Y]=meshgrid(x1,t1);
 mesh(X,Y,uwind)
 xlabel('x')
 ylabel('t')
 zlabel('u(x,t)')
 
 
 