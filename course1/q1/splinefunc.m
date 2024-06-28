%Spline function
function B=splinefunc(x)
if (x<=-2)
    B=0;

elseif (-2<=x && x<=-1)
    B=(x+2)^3;

elseif (-1<=x && x<=0)
    B=1+3*(x+1)+3*(x+1)^2-3*(x+1)^3;

elseif (0<=x && x<=1)
    B=1+3*(1-x)+3*(1-x)^2-3*(1-x)^3;
    
elseif (1<=x && x<=2)
    B=(2-x)^3;
    
else
    B=0;
end
end