G=@Gauss;
F=@testguassfunc;
G10=G(F,1);
real=sqrt(pi)*erf(1);
err=real-G10;
