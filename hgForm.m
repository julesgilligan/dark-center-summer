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

function [out] = hgForm(x,y,m,n)
global w0 z lambda;

E0 = 1000/sqrt(2*pi);
k0=(2*pi)/lambda;
z0 = pi*w0.^2./lambda;
if z ~= 0
    R = z*(1+(z0/z)^2);         % Radius of curviture of phase front
    circPhase = -1i*k0*(x.^2+y.^2)/(2*R);
else
    circPhase = 0;
end
w = w0*sqrt(1+(z^2/z0^2));  %Beam waist at z.

% out = exp(atan(z/z0)); % Complete HG equation

out = E0 * w0/w * hermitePoly(m,(sqrt(2).*x)/w) .* hermitePoly(n,(sqrt(2).*y)/w) .* exp(-(x.^2+y.^2)/w^2) .* exp(circPhase) * exp(-1i*k0*z) * exp(1i*(m+n+1)*atan(z/z0));

end