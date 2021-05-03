function [TTF] = ttf(centerslice_calibration, center_calibration, pixelsize)
%MTF function
%Gijs van Praagh - gijsvanpraagh@live.nl
%October 2019

%slice = slice of CCI phantom where calcifications are present (central
%calcification slice). Named LHD_centerslice in agatstonscore_thesis.m
%center_calc is center of the large, high density calcification. Named
%LHD_centercoordinates in agatstonscore_thesis.m
    
    factor = 6; %REFERENTIE NAAR PAPER NIELS
    [x,y] = size(centerslice_calibration);
    [xn,yn] = meshgrid(1:x,1:y);
    [xi,yi] = meshgrid(1:1/factor:x,1:1/factor:y);
    factor = size(xi,1)/size(xn,1);
    centerslice_int = interp2(xn, yn, centerslice_calibration, xi, yi, 'linear');
    pixelsize_int = pixelsize/factor;
    new_center = [center_calibration(1)*factor, center_calibration(2)*factor];
    
    xy_size = 21 / pixelsize_int;
    calc_ROI = centerslice_int(ceil(new_center(2)-xy_size):ceil(new_center(2)+xy_size),ceil(new_center(1)-xy_size):ceil(new_center(1)+xy_size));
    size_ROI = size(calc_ROI);
    
    ESF = zeros(1,ceil(size_ROI(1)/2));
    not_used_angles = 360 - abs(floor(angle_to_water-45)) + abs(ceil(angle_to_water+45));
    for angle = 1:360
        if angle >= floor(angle_to_water-45) && angle <= ceil(angle_to_water+45)
            continue
        end
        rot_ROI = imrotate(calc_ROI,angle,'bicubic','crop');
        Fr = rot_ROI(ceil(size_ROI(2)/2),ceil(size_ROI(1)/2):end)/(360-not_used_angles);
        ESF = ESF + Fr;
    end
    
    %ESF logistic function
    %mean background on homogeneous slice
    background = homogeneous_slice(floor(center_hs(1,2)-round(25/pixelsize)):floor(center_hs(1,2)+round(25/pixelsize)), floor(center_hs(1,1)-round(25/pixelsize)):floor(center_hs(1,1)+round(25/pixelsize)));
    mean_background = mean(mean(background));
    %mean calibration rod
    BW = imbinarize(centerslice_calibration,129);
    CC = bwconncomp(BW,4);                                  %find all connected components in the slice, i.e. structures. The last number is the connectivity
    numobjects = size(CC.PixelIdxList);                 %amount of structures in the slice
    indcenter = sub2ind([x,y],round(center_calibration(2)),round(center_calibration(1)));
    for lesion = 1:numobjects(2)
        les = CC.PixelIdxList{lesion};
        if ismember(indcenter,les)
            HU = centerslice_calibration(les);                                   %HU values in the lesion
            mean_calrod = mean(HU);
            break
        end
    end
    
    d = mean_background;
    a = 0.75*(mean_calrod-mean_background);
    b = size(ESF,2)/2;
    c = 0.15;
    x0 = [d,a,b,c];
    fermi = @(x,xdata) x(1) + (x(2) ./ (exp((xdata - x(3)) / x(4)) + 1));
        
    xdata = linspace(1,size(ESF,2),size(ESF,2));
    x = lsqcurvefit(fermi,x0,xdata,ESF);
    
    ESF = zeros(1,size(ESF,2));
    for i = 1:size(ESF,2)
        output = fermi(x,i);
        ESF(1,i) = output;
    end
        
    %TOT HIER GEKOMEN
    
    LSF = abs(diff(ESF));
    size_LSF = size(LSF);
    LSF = LSF.*hann(size_LSF(2))';
    TTF = abs(fftshift(fft(LSF)))/sum(LSF);
    size_TTF = size(TTF);
    TTF_x = linspace(0, 1/pixelsize_int, size_TTF(2));
    TTF_x = TTF_x(:,1:ceil(size_TTF(2)/2));
    TTF_y = TTF(:,ceil(size_TTF(2)/2)+1:end);
    TTF = [TTF_x;TTF_y];
%     MTF_0_5 = round(interp(0.5, MTF_y, MTF_x),2);
%     MTF_0_1 = round(interp(0.1, MTF_y, MTF_x),2);
    figure(1)
    plot(TTF_x,TTF_y)
    xlim([0 1.2])
    ylim([0 1])
    ylabel('TTF')
    xlabel('f (mm^{-1})')
end