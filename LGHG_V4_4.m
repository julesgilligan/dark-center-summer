%% Laguerre-Gaussian Beam amplitude and phase plot script.

% Version 4.4  -  Jules Gilligan
% created July 18, 2017
% Cross talk calculations based on HG V3.2.7

clearvars;
close all;

%% Parameters
global w0 z lambda;

lambda = 1.5e-6;
w0=4e-3;
z=40;
%% initializers
steps = 1000;
Rmax = 40e-3;
PHImax = 2*pi;
[R, PHI] = meshgrid(0:Rmax/(steps-1):Rmax,0:PHImax/(steps-1):PHImax);

%% Laguerre polynomial indices
l(1) = 1;  l(2) = 7;
p(1) = 0;  p(2) = 0;

%% ###### LG 1 #######
LG1 = LG_wave(l(1),p(1),R,PHI);

% total power on whole front < Rmax
tot_power(1) = sum(LG1.intensity(:),'omitnan');

% extract surfnorm vectors of phase front
[Nx1,Ny1,Nz1] = surfnorm(R.*cos(PHI),R.*sin(PHI),LG1.phaseC);

%% ###### LG 2 #######
for L = 1:2
    l(2) = L;   % allows incrementing second mode's l from 1 to 9
    LG2 = LG_wave(l(2),p(2),R,PHI);
    
    % total power on whole front < Rmax
    tot_power(2) = sum(LG2.intensity(:),'omitnan');
    
    % extract surfnorm vectors of phase front
    [Nx2,Ny2,Nz2] = surfnorm(R.*cos(PHI),R.*sin(PHI),LG2.phaseC);
    
    % Delta K - LG calc
    %   (1) extract phase norms (K vectors)
    %   (2) dot product (delK values)
    %   (3) plot
    
    % (2): Calculate Delta k dot product between normal vectors from the
    % two beam surfaces.
    
    % Fastest dot product method I can formulate .06s at (1,3:3) 3000stp
    CosTheta = Nx1.*Nx2 + Ny1.*Ny2 + Nz1.*Nz2;
    CosTheta(CosTheta > 1) = 1;
    
    % (3)
    if l(1) ~= l(2) % erratic behavior if the two modes' orders are equal
        deltak = acosd(CosTheta);
        figure('Name', 'Delta K plot','NumberTitle','off');
        surf(R.*cos(PHI),R.*sin(PHI),deltak); axis 'vis3d'; shading flat;
    end
    
    Pure(1,2) = 0;
    spower(1,2)=0;
    cross(1,2) =0;
    comb(1,2) = 0;
    
    for i=1:steps
        for j=1:steps
            if isfinite(LG1.intensity(i,j)) && isnan(LG2.intensity(i,j))
                Pure(1) = Pure(1) + LG1.intensity(i,j);
                comb(1) = comb(1) + LG1.intensity(i,j);
            elseif isnan(LG1.intensity(i,j)) && isfinite(LG2.intensity(i,j))
                Pure(2) = Pure(2) + LG2.intensity(i,j);
                comb(2) = comb(2) + LG2.intensity(i,j);
            elseif isfinite(LG2.intensity(i,j)) && isfinite(LG1.intensity(i,j))
                % Winner takes all method for spatially separating power
                if LG1.intensity(i,j) > LG2.intensity(i,j)
                    spower(1) = spower(1) + LG1.intensity(i,j) - LG2.intensity(i,j);
                    cross(1) = cross(1) + LG2.intensity(i,j);
                    comb(1) = comb(1) + LG1.intensity(i,j);
                elseif LG2.intensity(i,j) > LG1.intensity(i,j)
                    spower(2) = spower(2) + LG2.intensity(i,j) - LG1.intensity(i,j);
                    cross(2) = cross(2) + LG1.intensity(i,j);
                    comb(2) = comb(2) + LG2.intensity(i,j);
                end
            end
        end
    end
    
    % Checking to be sure acPow# = tot_power(#)
    acPow1 = Pure(1) + spower(1) + cross(1) + cross(2);
    acPow2 = Pure(2) + spower(2) + cross(1) + cross(2);
    
    % Cross talk calculation (ratio of noise against signal)
    CT(l(2),1) = (cross(1)/comb(1))*100;
    CT(l(2),2) = (cross(2)/comb(2))*100;
    
    % Separable power (%) by only the "pure" (non-overlapping) points.
    PureRatio(l(2),1) = Pure(1)/tot_power(1) * 100;
    PureRatio(l(2),2) = Pure(2)/tot_power(2) * 100;
end

%### Graph at the end
Intensity1 = LG1.intensityC;
Intensity2 = LG2.intensityC;
Phase1 = LG1.phaseC;
Phase2 = LG2.phaseC;
%###

figure(3);  movegui('southwest');
mesh(R.*cos(PHI) *1000,R.* sin(PHI) *1000, LG1.phase);
hold on;
surf(R.*cos(PHI) *1000,R.* sin(PHI) *1000, LG1.phaseC);
hold off;
title('Phase Cutoff by Intensity');
xlabel('x (mm)','FontSize',18); ylabel('y (mm)','FontSize',18); zlabel('Radians (\pi)','FontSize',18);
shading flat;
axis('vis3d');
rotate3d;

%% Graphing
Graphing(Intensity1,Intensity2,LG1.phaseC,LG2.phaseC,R,PHI);

function Graphing(A1, A2, P1, P2,R, PHI)

figure('Name','Intensity Plots','NumberTitle','off');  movegui('northwest');
subplot(1,2,1)
surf(R.*cos(PHI) *1000,R.* sin(PHI) *1000,A1);
title('LG1 Intensity');
xlabel('x (mm)','FontSize',18); ylabel('y (mm)','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
axis('vis3d');
subplot(1,2,2)
surf(R.*cos(PHI) *1000,R.* sin(PHI) *1000,A2);
title('LG2 Intensity');
xlabel('x (mm)','FontSize',18); ylabel('y (mm)','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
axis('vis3d');
rotate3d;

figure('Name','Phase Plots','NumberTitle','off');  movegui('northeast');
subplot(1,2,1)
surf(R.*cos(PHI) *1000,R.*sin(PHI) *1000,P1);
hold on;
surf(R.*cos(PHI) *1000,R.* sin(PHI) *1000,P2);
hold off;
title('LG1 Phase');
xlabel('x (mm)','FontSize',18); ylabel('y (mm)','FontSize',18); zlabel('Radians (\pi)','FontSize',18);
shading flat
axis('vis3d');
subplot(1,2,2)
surf(R.*cos(PHI) *1000,R.* sin(PHI) *1000,P2);
title('LG2 Phase');
xlabel('x (mm)','FontSize',18); ylabel('y (mm)','FontSize',18); zlabel('Radians (\pi)','FontSize',18);
shading flat
axis('vis3d');
rotate3d;

end
function moved = MidToZero(Input)
Zmax = max(Input(:));
Zmin = min(Input(:));
Zmid = (Zmax-Zmin)/2+Zmin;

moved = Input - Zmid;
end