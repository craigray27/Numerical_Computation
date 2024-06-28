function F=Func(x)
 F=zeros(2,1);
 F(1)=x(1)^2-x(2);
 F(2)=x(1)^2+x(2)^2-2;
end