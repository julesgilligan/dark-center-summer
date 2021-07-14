%-----------------------------------------------------------------------
%This function returns the complex field associated with an Hermite 
%Guassian beam. This function can be called in the main script using;
%
%hgForm(x,y,m,n)
%
%where the parameters are;
%x is the x-axis component along the screen,
%y is the y-axis component around the screen,
%m is the number of nulls of the HG beam along the x axis
%n is the number of nulls of the HG beam along the y axis
%
%-----------------------------------------------------------------------

function [out] = hgFormpp(x,y,m,n)
global w0;

out = hermitePoly(m,(sqrt(2).*x)/w0) .* hermitePoly(n,(sqrt(2).*y)/w0) .* exp(-(x.^2+y.^2)/w0^2);

end