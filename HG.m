%-----------------------------------------------------------------------
%This function returns the complex field associated with an Hermite 
%Guassian beam. This function can be called in the main script using;
%
%HG(w0,x,y,m,n)
%
%where the parameters are;
%w0 is the beam waist of the generated mode
%x is the horizontal cartesian component along the screen,
%y is the verticle cartesian component up the screen,
%m is the characteristic number in x of the HG beam 
%n is the characteristic number in y if the HG beam
%
%-----------------------------------------------------------------------


function [psi] = HG(w,x,y,m,n)

psi = hermitePoly(m,(x.*sqrt(2.))/w).*(exp((-x.^2)/(w.^2))).*(hermitePoly(n,(y.*sqrt(2.))/(w))).*(exp((-y.^2)/(w.^2)));

end

