function [NPS] = nps(homogeneous_slice,center_hs,pixelsize,radius,N_ROI,Nx,Ny)
%NPS function
%Gijs van Praagh - g.d.van.praagh@umcg.nl
%April 2021

    %2D NPS
    NPS_vector = zeros(Nx,Ny);
    if N_ROI > 1
        for i = 1:N_ROI-1
            NPS_ROI = homogeneous_slice(floor(center_hs(1,2)+((radius/2)/pixelsize)*sind((360/(N_ROI-1))*i)-((Nx-1)/2)):floor(center_hs(1,2)+((radius/2)/pixelsize)*sind((360/(N_ROI-1))*i)+((Nx-1)/2)), floor(center_hs(1,1)+((radius/2)/pixelsize)*cosd((360/(N_ROI-1))*i)-((Ny-1)/2)):floor(center_hs(1,1)+((radius/2)/pixelsize)*cosd((360/(N_ROI-1))*i)+((Ny-1)/2)));
            I_mean = mean(mean(NPS_ROI));
            
            F = NPS_ROI-I_mean;
            DFT = fftshift(fft2(F));
            DFT2 = (abs(DFT).^2)*((pixelsize^2)/(Nx*Ny));
            NPS_vector = NPS_vector + DFT2;
        end
    elseif N_ROI < 1
        disp('The amount of ROI''s chosen is too low')
        return
    else
        
    end
    
    %center ROI
    NPS_ROI = homogeneous_slice(floor(center_hs(1,2)-((Nx-1)/2)):floor(center_hs(1,2)+((Nx-1)/2)), floor(center_hs(1,1)-((Ny-1)/2)):floor(center_hs(1,1)+((Ny-1)/2)));
    I_mean = mean(mean(NPS_ROI));
    
    F = NPS_ROI-I_mean;
    DFT = fftshift(fft2(F));
    DFT2 = (abs(DFT).^2)*((pixelsize^2)/(Nx*Ny));
    NPS_vector = NPS_vector + DFT2;
    
    NPS_2D = NPS_vector/(N_ROI);
    %interpolate
    factor = 8; %must be an odd number
    [xn,yn] = meshgrid(1:Nx,1:Ny);
    [xi,yi] = meshgrid(1:1/factor:Nx,1:1/factor:Ny);
    factor = size(xi,1)/size(xn,1);
    Nx = Nx*factor;
    Ny = Ny*factor;
    
    NPS_2D = interp2(xn, yn, NPS_2D, xi, yi, 'cubic');
    %plot NPS for check
%     figure(1)
%     imshow(NPS_2D,[])
    
    %1D NPS
    NPS_1D = NPS_2D(ceil(Ny/2),ceil(Nx/2):end)/360;
    for angle = 1:359
        rot_NPS = imrotate(NPS_2D,angle,'bicubic','crop');
        Fr = rot_NPS(ceil(Ny/2),ceil(Nx/2):end)/360;
        NPS_1D = NPS_1D + Fr;
    end
    NPS_1D = round(NPS_1D,3);
    NPS_x = linspace(0,Nx*factor,Nx*factor) / (Nx*factor*(pixelsize/factor));
    NPS_x = NPS_x(1:ceil(Nx/2));
    NPS = [NPS_x, NPS_1D];
    
    %plot NPS
    figure(2)
    plot(NPS_x,NPS_1D,'b--o')
    xlabel('Radial frequency, f_{r} / mm')
    ylabel('NPS / (HU^{2} mm^{2})')
end