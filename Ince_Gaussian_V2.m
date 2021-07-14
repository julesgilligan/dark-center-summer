%% Ince-Gaussian Beam amplitude, phase plot, and cross talk script.

% Version 2  -  Jules Gilligan created July 24, 2017
% IG delta K histogram implementation based on HG 3.2.7 math.

clear vars;
close all;
%% initializers
steps = 892; %IDK why but at px=9e-6 and steps >892 there is a script error
% ig_gouy=zeros(1:2);
% tot_gouy=zeros(1:2);
%% Parameters
resx=steps; %xresolution
resy=steps; %yresolution
pixsize=9e-6; %SLM pixel si
apx=(resx*pixsize)/2; %aperture in x size based on slm pixel size
apy=(resy*pixsize)/2; %aperture in y size based on slm pixel size

[x,y] = meshgrid(-20e-3:40e-3/steps:20e-3,-20e-3:40e-3/steps:20e-3); % steps^2 grid (unit mmm)
lambda=1.5e-6;
w0=5e-3;
z=4;
parity=0;  %parity of the beam, 0 = EVEN  C;  1 = ODD  S;
ellip = 4; % ellipticity
%% Ince polynomial indices
p1 = 2;  p2 = 5;
m1 = 2;  m2 = 3;
%% Gaussian beam equation components.
E0 = 1000/sqrt(2*pi);
r = sqrt((x).^2+(y).^2);
k0 = (2*pi)/lambda;
z0 = pi*w0.^2./lambda;
R = z.*(1+(z0./z).^2);
w = w0.*sqrt(1+(z.^2./z0.^2));

plnwv = -1i*k0*z;
normterm = (w0/w);
circphase = -1i.*pi.*r.^2./(lambda*R);
gauss_gouy = 1i*atan(z/z0);

ig_gouy1 = p1 + 1;  % IG Gouy phase
tot_gouy1 = gauss_gouy*ig_gouy1; % Gaussian and IG Gouy phase together

ig_gouy2 = p2 + 1;  % IG Gouy phase
tot_gouy2 = gauss_gouy*ig_gouy2; % Gaussian and IG Gouy phase together

igmode1 = IG(apx,apy,parity,p1,m1,ellip,w0,k0,0,x,y,pixsize);
IG1 = E0.*igmode1;
IG1phase = angle(IG1);
I1 = abs(IG1).^2;

igmode2 = IG(apx,apy,parity,p2,m2,ellip,w0,k0,0,x,y,pixsize);
IG2 = E0.*normterm.*igmode2.*exp(circphase).*exp(plnwv).*exp(tot_gouy1);
IG2phase = angle(IG2);
I2 = abs(IG2).^2;

%%
% M1 is the 1/e^4 value of the maximum intensity value of mode 1.
% If intensity value of a point is less than M1 set it to NaN.
M1 = max(I1(:))*1/exp(4);
I1(I1<M1) = NaN;
Chopped1 = NaN(size(IG1phase));
Chopped1(I1>M1) = IG1phase(I1>M1);

% Counts total number of nonzero power grid points
II1 = isfinite(I1);
figure(4); subplot(2,3,1);
spy(sparse(II1));
scount1T = nnz(II1);

% M2 is the 1/e^4 value of the maximum intensity value of mode 2.
% If intensity value of a point is less than M2 set it to zero.
M2 = max(I2(:))*1/exp(4);
I2(I2<M2) = NaN;
Chopped2 = NaN(size(IG2phase));
Chopped2(I2>M2) = IG2phase(I2>M2);

% Counts total number of nonzero power grid points
II2 = isfinite(I2);
figure(4); subplot(2,3,3);
spy(sparse(II2));
scount2T = nnz(II2);

%% Delta k
% Calculate dot product between normal vectors from the two beam surfaces
% and calculate magnitude of each mode.  All needed for calculating angle
% between two vectors.

[Nx1,Ny1,Nz1] = surfnorm(x,y,Chopped1);

[Nx2,Ny2,Nz2] = surfnorm(x,y,Chopped2);

CosTheta = Nx1.*Nx2 + Ny1.*Ny2 + Nz1.*Nz2;
% Approximate floating point numbers almost 1 or 0 to be 1 or 0
% CosTheta(CosTheta > .999) = 1; % acosd(.999)=2.6 aprox=0=acosd(1)
% CosTheta(CosTheta < .05 ) = 0; % acosd(.05)=87 aprox=90=acosd(0)

deltak = acosd(CosTheta);
deltak(deltak > 100) = NaN;
deltakmax = max(deltak(:));

scatter = deltak;
K= zeros(steps + 1);
K(isfinite(scatter(:))) = 1 + scatter(isfinite(scatter(:)));
figure(4); subplot(2,3,2)
spy(K);
title('\Delta K');

% plot the HISTOGRAM
figure(3); %subplot(3,3,n(2));
H1 = histogram(deltak,20);
title(sprintf('$$IG_{(%d,%d)} \\Delta$$ k',p1,m1),'Interpreter','latex');
xlabel('\Delta k','FontSize',10);
ylabel('# \Delta k','FontSize',10);

figure(6);
surf(x,y,deltak);
axis vis3d; shading interp;

%% Cross Talk
Pure(1) = 0;
Pure(2) = 0;
spower(1)=0;
spower(2)=0;
cross(1)=0;
cross(2)=0;
comb(2)=0;
for i=1:steps + 1
    for j=1:steps + 1
        if isfinite(Chopped1(i,j)) && isnan(Chopped2(i,j))
            Pure(1) = Pure(1) + I1(i,j);
            comb(1) = comb(1) + I1(i,j);
        elseif isnan(Chopped1(i,j)) && isfinite(Chopped2(i,j))
            Pure(2) = Pure(2) + I2(i,j);
            comb(2) = comb(2) + I2(i,j);
        elseif isfinite(Chopped2(i,j)) && isfinite(Chopped1(i,j))
            % Winner takes all power separation method
            if I1(i,j) > I2(i,j)
                spower(1) = spower(1) + I1(i,j) - I2(i,j);
                cross(1) = cross(1) + I2(i,j);
                comb(1) = comb(1) + I1(i,j);
            elseif I2(i,j) > I1(i,j)
                spower(2) = spower(2) + I2(i,j) - I1(i,j);
                cross(2) = cross(2) + I1(i,j);
                comb(2) = comb(2) + I2(i,j);
            end
        end
    end
end

Phase1 = Chopped1; %Phase
Phase2 = Chopped2; %Phase
pplate1 = angle(igmode1)-3*pi; %Phase plate
pplate2 = angle(igmode2)-3*pi; %Phase plate

figure(1);  movegui('northwest');
subplot(1,2,1)
surf(x,y,I1);
title(strcat('G1 Intensity'));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
rotate3d;

subplot(1,2,2)
surf(x,y,I2);
title(strcat('G2 Intensity'));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
rotate3d;

figure(2);  movegui('northeast');

subplot(1,2,1)
surf(x,y,Phase1);
hold on;
surf(x,y,pplate1);
hold off;
title(strcat('G1 Phase'));
xlabel('x','FontSize',24); ylabel('y','FontSize',24); zlabel('Radians','FontSize',24);
shading flat
rotate3d;

subplot(1,2,2)
surf(x,y,Phase2);
hold on;
surf(x,y,pplate2);
hold off;
title(strcat('G2 Phase'));
xlabel('x','FontSize',24); ylabel('y','FontSize',24); zlabel('Radians','FontSize',24);
shading flat
rotate3d;
