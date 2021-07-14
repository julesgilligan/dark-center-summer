%% Hermite-Gaussian Beam amplitude, phase plot, and cross talk script.

% Version 3.2.8  -  Jules Gilligan created July 21, 2017
% HG_wave class implementation attempt and therefore code simplification.

%% initializers
clear all;
close all;

steps = 1500;
hg_gouy=zeros(1:2);
tot_gouy=zeros(1:2);
%% Parameters
global w0 z lambda;
[x,y] = meshgrid(-20e-3:40e-3/(steps-1):20e-3,-20e-3:40e-3/(steps-1):20e-3); % steps^2 grid (unit mmm)

lambda=1.5e-6;
w0=5e-3; % beam waist
z=4;  % At certain z values there is a suddent 2*pi phase shift on the
% phase front. ???
%% Hermite polynomial indices
m(1) = 1;  m(2) = 1;    m(3) = 3; % m,n(3) are not used in this version
n(1) = 1;  n(2) = 1;    n(3) = 2;
%% Gaussian beam equation components.
E0 = 1000/sqrt(2*pi);
r = sqrt((x).^2+(y).^2);
k0 = (2*pi)/lambda;
z0 = pi*w0.^2./lambda;
R = z.*(1+(z0./z).^2);
w = w0.*sqrt(1+(z.^2./z0.^2));

% Helpful reminder of the different phase parts.
%   plnwv = -1i*k0*z;
%   normterm = (w0/w);
%   circphase = -1i.*pi.*r.^2./(lambda*R);
%   gauss_gouy = 1i*atan(z/z0);

%% Some input checks
if steps > 1500
    input('Step count is high. Possible high runtime. \n');
end
if w0 <= 1e-6
    error('Beam waist is smaller than a single facet. Improper dimension.');
end
if mod(m(1),1)==0 && mod(n(1),1)==0 && mod(m(2),1)==0 && mod(n(2),1)==0 && mod(m(3),1)==0 && mod(n(3),1)==0
    m(1) = abs(m(1));   m(2) = abs(m(2));   m(3) = abs(m(3));
    n(1) = abs(n(1));   n(2) = abs(n(2));   n(3) = abs(n(3));
end

%% ###### HG 1 #######
HG1 = HG_wave(m(1),n(1),x,y);
I1 = HG1.intensity;     % Store intensity of HG1 beam.
tot_Power1 = sum(HG1.intensityC(:),'omitnan'); %Power of the 'chopped' beam

% Normalize HG1 intensity by dividing through by the entire beam's power.
I1 = I1/tot_Power1;
M1 = max(I1(:))/exp(4); % max/e^4 is an acceptable threshold
I1chopped = zeros(size(I1));
I1chopped(I1>M1) = I1(I1>M1);  % Where I1>M1 inside mode boundary, set I1chopped (which is zero) to the I1 value.

tot_Power1 = sum(I1chopped(:)); % Re-tally power; should equal 1 (100%) now

% Visualize all non-zero points
C1 = zeros(steps);
C1(isfinite(I1chopped(:))) = I1chopped(isfinite(I1chopped(:)));
figure(1); subplot(2,3,1);
spy(C1);
title('Full First Intensity Area');

% extract surfnorm vectors of phase front
%#%#%#%#%##%#%#%#%#%#%#%#%%#%#%#%#%##%#%#%#%#%#%#%#%
[Nx1,Ny1,Nz1] = surfnorm(x,y,HG1.phaseC); %#%#%#%#%#
%#%#%#%#%##%#%#%#%#%#%#%#%%#%#%#%#%##%#%#%#%#%#%#%#%
%% ###### HG 2 #######
for N = 1:9
    n(2) = N;
    
    if mod(n(2),1) == 0     % non-integer number screen
        n(2) = abs(n(2));
    else
        error('Non-integer n value.');
    end
        
        
    HG2 = HG_wave(m(2),n(2),x,y);
    I2 = HG2.intensity;     % Store intensity of HG1 beam.
    tot_Power2 = sum(HG2.intensityC(:),'omitnan'); %Power of the beam
    
  % Normalize HG1 intensity by dividing through by the entire beam's power.
    I2 = I2/tot_Power2;
    M2 = max(I2(:))/exp(4); % max/e^4 is an acceptable threshold
    I2chopped = zeros(size(I2));
    I2chopped(I2>M2) = I2(I2>M2);
    
    tot_Power2 = sum(I2chopped(:)); % Re-tally power should equal 100% now
    
    % Visualize all non-zero points
    C2= zeros(steps);
    C2(isfinite(I2chopped(:))) = I2chopped(isfinite(I2chopped(:)));
    figure(1); subplot(2,3,3)
    spy(C2);
    title('Full Second Intensity Area');
    C1(isfinite(I1chopped(:))) = I1chopped(isfinite(I1chopped(:)));
    
    % Calculate overlapping areas
    C3= zeros(steps);
    for e=1:steps
        for f=1:steps
            if C1(e,f) && C2(e,f)
                C1(e,f) = 0;
                C2(e,f) = 0;
                C3(e,f) = 1;
            end
        end
    end
    subplot(2,3,4);
    spy(C1); title('Pure 1');
    subplot(2,3,5)
    spy(C3); title('Overlap: ''Impure'' Power');
    subplot(2,3,6);
    spy(C2); title('Pure 2');
    
    
    % extract surfnorm vectors of phase front
    %#%#%#%#%##%#%#%#%#%#%#%#%%#%#%#%#%##%#%#%#%#%#%#%#%
    [Nx2,Ny2,Nz2] = surfnorm(x,y,HG2.phaseC); %#%#%#%#%#
    %#%#%#%#%##%#%#%#%#%#%#%#%%#%#%#%#%##%#%#%#%#%#%#%#%
    
    % Calculate Delta k dot product between normal vectors from the two
    % beam surfaces. Then round erroneous data towards reasonable values.
    
    % Fastest dot product method I can formulate .06s at (1,3:3) 3000stp
    CosTheta = Nx1.*Nx2 + Ny1.*Ny2 + Nz1.*Nz2;
    
    % Approximate floating point numbers almost 1 or 0 to be 1 or 0
    CosTheta(CosTheta > .999) = 1; % acosd(.999)=2.6 aprox=0=acosd(1)
    CosTheta(CosTheta < .05 ) = 0; % acosd(.05)=87 aprox=90=acosd(0)
    
    
    deltak = acosd(CosTheta);
    deltak(deltak > 100) = NaN;
    deltakmax = max(deltak(:));
    
    scatter = deltak;
    K= zeros(steps);
    K(isfinite(scatter(:))) = 1 + scatter(isfinite(scatter(:)));
    figure(1); subplot(2,3,2)
    spy(K);
    title('\Delta K');
    
    % plot the HISTOGRAM
    if n(2) > 0 && n(2) < 10
        figure(2); subplot(3,3,n(2));
        H1 = histogram(deltak,5);
        title(sprintf('$$HG_{(1,%d)} \\Delta$$ k',n(2)),'Interpreter','latex');
        xlabel('\Delta k','FontSize',10);
        ylabel('# \Delta k','FontSize',10);
    end
    
    % Calculate cross talk based on non-overlapping and winner takes
    % methods.
    Pure(1,2) = 0;
    spower(1,2)=0;
    cross(1,2)=0;
    comb(1,2)=0;
    for i=1:steps
        for j=1:steps
            if isfinite(I1chopped(i,j)) && isnan(I2chopped(i,j))
                Pure(1) = Pure(1) + I1chopped(i,j);
                comb(1) = comb(1) + I1chopped(i,j);
            elseif isnan(I1chopped(i,j)) && isfinite(I2chopped(i,j))
                Pure(2) = Pure(2) + I2chopped(i,j);
                comb(2) = comb(2) + I2chopped(i,j);
            elseif isfinite(I1chopped(i,j)) && isfinite(I2chopped(i,j))
                % Winner takes all method for spatially separating modes.
                if I1(i,j) > I2(i,j)
                    spower(1) = spower(1) + I1chopped(i,j) - I2chopped(i,j);
                    cross(1) = cross(1) + I2chopped(i,j);
                    comb(1) = comb(1) + I1chopped(i,j);
                elseif I2(i,j) > I1(i,j)
                    spower(2) = spower(2) + I2chopped(i,j) - I1chopped(i,j);
                    cross(2) = cross(2) + I1chopped(i,j);
                    comb(2) = comb(2) + I2chopped(i,j);
                end
            end
        end
    end
    
    % Checking to be sure acPow# = tot_power(#)
    acPow1 = Pure(1) + spower(1) + cross(1) + cross(2);
    acPow2 = Pure(2) + spower(2) + cross(1) + cross(2);
    
    % Cross talk calculation (ratio of noise against signal)
    CT(n(2),1) = (cross(1)/comb(1))*100;
    CT(n(2),2) = (cross(2)/comb(2))*100;
    
    % Separable power (%) by only the "pure" (non-overlapping) points.
    PureRatio(n(2),1) = Pure(1)/tot_Power1 *100;
    PureRatio(n(2),2) = Pure(2)/tot_Power2 *100;
    
end


%### Graphed at the end
Phase1 = HG1.phaseC;
pplate1 = angle(HG1.pplate)-3*pi;

Phase2 = HG2.phaseC;
pplate2 = angle(HG2.pplate)-3*pi;
%###

% Old code to write SimData to an excel document
%   xlswrite('Simulation Data for MUXed HG modes.xlsx', CT(1), 'A1:A10');
%   xlswrite('Simulation Data for MUXed HG modes.xlsx', CT(2), 'B1:B10');

%------------%
%  Graphing  %
%------------%

%Graphs based on x and y
%       2 Amplitude graphs A1 and A2 
%       2 Phase graphs P1 and P2 
%       2 Phase Plate graphs PP1 and PP2

figure(3);

subplot(2,2,1)
surf(x,y,I1);
title(sprintf('HG_{(%d,%d)} Intensity',m(1),n(1)));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
rotate3d;

subplot(2,2,2)
surf(x,y,I2);
title(sprintf('HG_{(%d,%d)} Intensity',m(2),n(2)));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
axis([-.02,.02, -.02,.02]);
rotate3d;

test = NaN(steps);
test(I1>I2)=I1(I1>I2);
test(I2>I1)=I2(I2>I1);

subplot(2,2,3)
surf(x,y,test);
title(sprintf('HG_{(%d,%d)} & HG_{(%d,%d)} Intensity',m(1),n(1),m(2),n(2)));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('W/m^2','FontSize',18);
shading flat
axis([-.02,.02, -.02,.02]);
rotate3d;

% phases = figure('Name','Phase Plots');
figure(4);

zm = max(pplate1(:));
subplot(1,2,1)
surf(x,y,Phase1);
hold on;
surf(x,y,pplate1);
hold off;
title(sprintf('HG_{(%d,%d)} Phase',m(1),n(1)));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('Radians','FontSize',18);
shading flat
axis vis3d; rotate3d;


subplot(1,2,2)
surf(x,y,Phase2);
hold on;
surf(x,y,pplate2);
hold off;
title(sprintf('HG_{(%d,%d)} Phase',m(2),n(2)));
xlabel('x','FontSize',18); ylabel('y','FontSize',18); zlabel('Radians','FontSize',18);
shading flat
axis vis3d;
rotate3d;

%      move gui is taking 12 seconds on first call; even on an empty figure
%figure(ints); movegui('southwest');
%figure(thetaPlot); movegui('south');
% figure(1); movegui('southeast');
% figure(2); movegui('southwest');
% figure(3); movegui('northwest');
% figure(4); movegui('northeast');
