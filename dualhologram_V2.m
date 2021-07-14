%% Genenerates hologram for LG, HG and mode sorter elements.
% Changed variable names.  Generating new dual hologram that have either of
% the two holograms in each region only.  Where ever one mode has higher
% power that mode's hologram is produced.
% clc;
clear;close all;
global w0 z lambda;

% General pamarmeters
w0=0.0008;%beamwaist
w01=0.0004;
lam=632.8*10^(-9); %wavelength of the light in meters
resx=1024; %xresolution
resy=768; %yresolution
pixsize=9.0*10^(-6); %SLM pixel si
apx=(resx*pixsize)/2; %aperture in x size based on slm pixel size
apy=(resy*pixsize)/2; %aperture in y size based on slm pixel size
lambda = lam;
z = 0;
z0 = pi*w0.^2./lambda;
w = w0.*sqrt(1+(z.^2./z0.^2));

%Choose the hologram type you wish to generate.
mode=0; % 0 for HG, 1 for LG,2 for IG, 3 for Modesorter element 1, 4 for Modesorter element 2.20

%gratingform=exp(i*x*h);

%HG mode specific parameters
m=1; mh=1; %m quantum number of the Hermite Gaussian Mode
n=1; nh=3; %n quantum number of the Hermite Gaussian Mode
% %
% % %LG mode specific parameters
% % l=1; ll=2; %azimthal quantum number (l) for Laguerre Gaussian Mode
% % p=0; pl=0;%radial quantum number (p) for Laguerre Gaussian Mode
% %
% % %IG mode specific parameters
% % p1=1; p1i=9;
% % m1=1; m1i=9;
% % parity=0;
% % ellip=4;

%Modesorter Specific parameters
%d=0.012; %defines the size of the width transformed line
%a=d/(2*pi); %the scaling parameter that is required for setting the width
%b=0.00477; %the parameter that define the height of the transformed beam
%f=0.5; %the focal length of the lens integrated, or between sorter elements
%lens=0; %defines wether the mode sorter elements have an integrated lens 0 for false, 1 for true.

%Parameters defining wedge function of grating, which seperates the beam
%into the first order
angx=0; %defines the angular tilt along the x axis
angy=1.5; %defines the angulat tilt along the y axis
angx1=1.5; %defines the angular tilt along the x axis
angy1=0; %defines the angulat tilt along the y axis

%% Build grid for computing the 2d data
[X,Y] = meshgrid(-apx:pixsize:apx,-apy:pixsize:apy);

%% Builds base grating resulting from a particular wedge function
wedge1=wedge(X,Y,angx,angy,apx,apy,lam);
wedge2=wedge(X,Y,angx1,angy1,apx,apy,lam);

%% Plots the desired hologram on second monitor if attached.

% Determine the screen resolution and number of monitors attached. If
% second monitor attached then hologram fills second monitor.
mp= get(groot,'MonitorPositions');
scrsz=mp(1,:); %resolution of of main monitor
nofs=1;     % no. of monitors attached

close all   % closes all figures' windows
if nofs(1) == 1
    
    if mode == 0
        %% HG Hologram Grating Creation
        
        hgmode1=HG(w0,X,Y,m,n); %builds hg mode array given the above parameters.
        
        hgmode2=HG(w0,X,Y,mh,nh);
        
        
        %plots the intensity profile for this partiular mode
        I1 = abs(hgmode1).^2;
        I2 = abs(hgmode2).^2;
        figure(4); movegui southeast
        subplot(1,2,1); surf(X,Y,I1); shading interp;
        subplot(1,2,2); surf(X,Y,I2); shading interp;
        
        %Assures the peak intensity of any mode is equal to 1.
        I1sc = I1/max(I1(:));
        I2sc = I2/max(I2(:));
        
        % Prune grating1 for increased contrast in significant portions
        grating1=angle(hgmode1.*wedge1);
        grating1=(3*grating1)+pi;   %picks only middle 3rd for max contrast
        grating1(grating1<0)=0;
        grating1(grating1>2*pi)=2*pi;
        holo1=I1sc.*grating1;
        
        % Prune grating1 for increased contrast in significant portions
        grating2=angle(hgmode2.*wedge2);
        grating2=(3*grating2)+pi;   %picks only middle 3rd for max contrast
        grating2(grating2<0)=0;
        grating2(grating2>2*pi)=2*pi;
        holo2=I2sc.*grating2;
       
        
        regions = NaN(resy+1,resx+1);
        for i=1:resy+1
            for j=1:resx+1
                if I2sc(i,j) > I1sc(i,j)
                      regions(i,j) = holo2(i,j);
                elseif I1sc(i,j) > I2sc(i,j)
                      regions(i,j) = holo1(i,j);
                end
            end
        end
        
        figure(6); 
        subplot(1,2,1);
        surf(X,Y,I1sc); hold on; 
        surf(X,Y,I2sc); hold off; shading interp; 
        subplot(1,2,2); surf(X,Y,regions); shading interp;
        
        
        c2=holo1+holo2;
        
        figure(1);
        subplot(2,2,1)
        imagesc(holo1);
        axis equal; axis off;
        subplot(2,2,2)
        imagesc(holo2);
        axis equal; axis off;
        colormap(gray);
        subplot(2,2,3)
        imagesc(c2);
        axis equal; axis off;
        colormap(gray);
        
        d=holo1(:,256:768);
        d1=holo2(:,256:767);
        d2=horzcat(d,d1);
        figure('units','normalized','outerposition',[0 0 1 1],'menubar', 'none')
        figure(2);
        set(gca,'Position',[0 0 1 1])
        imagesc(d2);
        colormap(gray);
        axis equal; axis off;
        saveas(gcf,'Double_Hologram','jpg');
        
        %plots the hologram for this mode, combined with the required grating.
        figure(3);
        imagesc(regions);colormap(gray); axis off; axis equal;
        saveas(gcf,'DeMUX_Hologram','jpg');
        %{
    elseif mode ==1
        
        %builds lg mode array given the above parameters.
        lgmode=LG(w0,X,Y,l,p);
        lgmode1=LG(w1,X,Y,ll,pl);
        
        %plots the intensity profile for this partiular mode
        
        x=abs(lgmode);
        x1=abs(lgmode1);
        
        %plots the hologram for this mode, combined with the required gating.
        y=angle(lgmode.*grating);
        z=angle(lgmode1.*grating1);
        y1=(3*y)+3.14;
        i=y1<0;
        y1(i)=0;
        o=y1>6.28;
        y1(o)=6.28;
        z1=(3*z)+3.14;
        i1=z1<0;
        z1(i1)=0;
        o1=z1>6.28;
        z1(o1)=6.28;
        c=x.*y1;
        %c=x.*y;
        c1=x1.*z1;
        %c1=x.*z;
        c2=c+c1;
        fig1=figure;
        imagesc(c2);
        colormap(gray);
        axis off;
        d=c(:,256:768);
        d1=c1(:,256:767);
        d2=horzcat(d,d1);
        fig1=figure;
        imagesc(d2);
        colormap(gray);
        axis off;
        
    elseif mode ==2
        
        %builds IG mode array given the above parameters.
        [IGB]=Ince_Gaussian2(apx,apy,parity,p1,m1,ellip,w0,(2*pi/632.8e-9),0,X,Y,pixsize);
        
        %plots the intensity profile for this particular mode
        
        
        MAX=IGB(1,1);
        for i=1:769
            for j=1:1025
                if MAX<= IGB(i,j);
                    MAX=IGB(i,j);
                end
            end
        end
        IGB2=IGB/MAX;
        
        [IGB1]=Ince_Gaussian2(apx,apy,parity,p1i,m1i,ellip,w1,(2*pi/632.8e-9),0,X,Y,pixsize);
        MAX1=IGB1(1,1);
        for i=1:769
            for j=1:1025
                if MAX1<= IGB1(i,j);
                    MAX1=IGB1(i,j);
                end
            end
        end
        IGB3=IGB1/MAX1;
        
        x=abs(IGB2);
        x1=abs(IGB3);
        %plots the hologram for this mode, combined with the required grating.
        fig2=figure;
        y=angle(IGB2.*grating);
        z=angle(IGB3.*grating1);
        y1=(3*y)+3.14;
        i=y1<0;
        y1(i)=0;
        o=y1>6.28;
        y1(o)=6.28;
        z1=(3*z)+3.14;
        i1=z1<0;
        z1(i1)=0;
        o1=z1>6.28;
        z1(o1)=6.28;
        c=x.*y1;
        %c=x.*y;
        c1=x1.*z1;
        %c1=x.*z;
        c2=c+c1;
        fig1=figure;
        imagesc(c2);
        colormap(gray);
        axis off;
        d=c(:,256:768);
        d1=c1(:,256:767);
        d2=horzcat(d,d1);
        fig2=figure;
        imagesc(d2);
        colormap(gray);
        axis off;
        
        %plots the hologram for tilt incidence amd higher power
        
        
    elseif mode ==3
        
        %plots the hologram for the first mode-sorter element, combined with the required gating.
        
        element1=ms1(X,Y,a,b,lam,f,0);
        
        fig2= figure;
        imagesc(angle((exp(1i*element1)).*grating));
        colormap(gray);
        axis off;
        
        
    elseif mode ==4
        
        %plots the hologram for the second mode-sorter element, combined with the required gating.
        
        
        element2=ms2(X,Y,a,b,lam,f,0);
        fig1= figure;
        imagesc(angle((exp(1i*element2)).*grating));
        colormap(gray);
        axis off;
        %}
    end
    %{
elseif nofs(1) == 2
    
    scrsz2=mp(2,:)
    
    if mode ==0
        
        hgmode=HG(w0,X,Y,m,n);
        
        %plots the intensity profile for this partiular mode on main screen
        fig1=figure;
        imagesc(abs(hgmode));
        
        colormap(gray);
        
        %plots the hologram for this mode, combined with the required gating,
        %filling the second monitor attached to the computer.
        fig2=figure;
        imagesc(angle(hgmode.*grating));
        set(fig2, 'menubar', 'none');
        set(fig2, 'position', [scrsz2(1),scrsz2(2),(scrsz2(3)),scrsz2(4)]);
        colormap(gray);
        truesize(fig2)
        axis off;
        
    elseif mode ==1
        
        lgmode=LG(w0,X,Y,l,p);
        %plots the intensity profile for this partiular mode on main screen
        fig1=figure;
        imagesc(abs(lgmode));
        colormap(gray);
        
        %plots the hologram for this mode, combined with the required gating,
        %filling the second monitor attached to the computer.
        fig2=figure;
        imagesc(angle(lgmode.*grating));
        set(fig2, 'menubar', 'none');
        set(fig2, 'position', [scrsz2(1),scrsz2(2),(scrsz2(3)),scrsz2(4)]);
        
        colormap(gray);
        axis off;
        
    elseif mode ==2
        
        %plots the hologram for the first mode-sorter element, combined with
        %the required gating. One displayed on each screen, second monitor is
        %filled with hologram.
        
        element1=ms1(X,Y,a,b,lam,f,0);
        fig1= figure;
        imagesc(angle((exp(1i*element1)).*grating));
        colormap(gray);
        axis off;
        
        fig2=figure;
        imagesc(angle((exp(1i*element1)).*grating));
        set(fig2, 'menubar', 'none');
        set(fig2, 'position', [scrsz2(1),scrsz2(2),(scrsz2(3)),scrsz2(4)]);
        
        colormap(gray);
        axis off;
        
    elseif mode ==3
        
        %plots the hologram for the second mode-sorter element, combined with
        %the required gating, combined with the required gating. One displayed
        %on each screen, second monitor is filled with hologram.
        
        element2=ms2(X,Y,a,b,lam,f,0);
        fig1= figure;
        imagesc(angle((exp(1i*element2)).*grating));
        colormap(gray);
        axis off;
        
        fig2=figure;
        imagesc(angle((exp(1i*element2)).*grating));
        set(fig2, 'menubar', 'none');
        set(fig2, 'position', [scrsz2(1),scrsz2(2),(scrsz2(3)),scrsz2(4)]);
        
        colormap(gray);
        axis off;
        
    end
    %}
end

