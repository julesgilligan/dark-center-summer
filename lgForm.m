%-----------------------------------------------------------------------
%This function returns the complex field associated with an Laguerre 
%Guassian beam. This function can be called in the main script using;
%
%lgForm(r,phi,l,p)
%
%where the parameters are;
%r is the radial component along the screen,
%phi is the angular component around the screen,
%l is the azimuthal winding number of the LG beam 
%p is the radial number of the LG beam
%
%-----------------------------------------------------------------------

function out = lgForm(r,phi,l,p)
global w0 z lambda;

% E0 = 1000/sqrt(2*pi);
E0 = factorial(p)*sqrt(2/(w0^2*pi*factorial(p)*factorial(abs(l)+p)));
k0=(2*pi)/lambda;
z0 = pi*w0.^2./lambda;
R = z*(1+(z0/z)^2);         % Radius of curviture of phase front
w = w0*sqrt(1+(z^2/z0^2));  %Beam waist at z.

amplitude = E0 * w0/w .* (r*sqrt(2)/w).^(abs(l)) .* exp(-r.^2/w^2) .* (LaguerrePoly([p abs(l)],(2.*(r./w).^2)));

out = amplitude .* exp(-1i.*k0.*r.^2./(2.*R) + -1i*l.*phi + -1i*k0.*z + 1i*(2.*p+abs(l)+1).*atan(z/z0));
  
end