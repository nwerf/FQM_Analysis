function [Data_Ag, lesion_info_Ag, sumdata_Ag, pixelsize, TTF, NPS] = CACscore(scans, method)
%Calcium scoring script
%Authors: Gijs van Praagh (g.d.van.praagh@umcg.nl) /
%Niels van der Werf (n.vanderwerf@erasmusmc.nl)
%May 2021
%version 1.0

%This function will calculate the Agatston, volume, and mass score from a
%CT scan of the commercially available QRM phantom with cardiac
%calcification insert. Besides that it also calculates some image quality
%measurements. "Data_Ag" contains all general information, the total CAC
%scores, and the image quality measurements. "lesion_info_Ag" contains all
%information of every structure in all slice and can be used for
%troubleshooting. "sumdata_Ag" contains all scores per calcification.

%Possible answers for input argument 'method': 'literatureCAC' -> will use
%scoring methods as described in literature. This is also the default
%setting if no input arg is given. 'vendorspecific' -> gets vendor name
%from dicom header and uses the vendor specific scoring method. 'Philips'
%-> will use Philips specific method. 'Canon' -> will use Canon specific
%method. 'Siemens' -> will use Siemens specific method. 'GE' -> will use GE
%specific method. For more information about the methods and the validation
%of all methods, see
%https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.14912

%When using this script, please cite our article: van Praagh GD, van der
%Werf NR, Wang J, van Ommen F, Poelhekken K, Slart RH, Fleischmann D,
%Greuter MJ, Leiner T, Willemink MJ. Fully Automated Quantification Method
%(FQM) of Coronary Calcium in an Anthropomorphic Phantom. Med Phys. 2021
%May 1. doi: 10.1002/mp.14912. Epub ahead of print. PMID: 33932026.

    tic
    
    if ~exist('method','var')
        method = 'literatureCAC';
    end
    
    dicomfiles = dir(scans);
    try
        info = dicominfo(fullfile(scans,dicomfiles(4).name));   %load the dicom header into a struct
        info2 = dicominfo(fullfile(scans,dicomfiles(5).name));
        st = abs(info.SliceLocation - info2.SliceLocation);     %slice thickness depending on increment/overlapping slices
        slice_thickness = info.SliceThickness;                  %real slice thickness
        kVp = info.KVP;                                         %tube voltage used
        pixelsize = info.PixelSpacing(1);                       %pixelsize of the CT images in mm
        fn = info.Filename;                                     %filename
        mf = info.Manufacturer;                                 %manufacturer name
        description = info.SeriesDescription;                   %description of the scan
%For our research it is important to know the IR level used. To keep
%'sumdata_Ag' (later defined) a matrix and not a cell, the following if
%statement is used. Because in our scans the mAs and IR level was not found
%correctly in the dicomheader of GE we got it out of the description.
        if contains(mf,'GE MEDICAL SYSTEMS') || contains(mf,'Imatron')
            try
                mAs = (info.ExposureTime/1000)*info.XRayTubeCurrent;
            catch
                mAs = (info.ExposureTime/1000)*info.XrayTubeCurrent;
            end
            try
                level = info.Private_0053_1043;
                level = str2double(level);
            catch
                level = 0;
            end
        elseif contains(mf,'Philips')
            mAs = info.Exposure;
            CTDIvol = info.CTDIvol;
            level = info.Private_01f7_109b;
            CK = info.ConvolutionKernel;
            if isempty(level)
                if contains(CK,'XCA')
                    level = 0;
                elseif contains(CK,'IMR')
                    level = 12;
                end
            end
        else
            mAs = info.Exposure;
            CTDIvol = info.CTDIvol;
            kernel = info.ConvolutionKernel;
            if kernel == "Qr36d" || kernel == "Qr36f" || kernel == "FC12" || kernel == "FC11" || kernel == "XCA" || kernel == "Qr32d" || kernel == "Qr44d"
                level = 0;
            elseif kernel == "Qr36d\2"
                level = 2;
            elseif kernel == "Qr36d\3"
                level = 3;
            elseif kernel == "Qr36d\4"
                level = 4;
            elseif kernel == "Sa36f"
                level = 10;
            elseif kernel == "Sharp R62"
                level = 11;
            else
                level = 99;
            end
            if kernel == "FC12" || kernel == "CARDIAC"
                sort = 1;
            elseif kernel == "FC11"
                sort = 2;
            elseif kernel == "FC17"
                sort = 3;
            else
                sort = 99;
            end
        end
    catch
        info = dicominfo(fullfile(scans,dicomfiles(3).name));
        st = abs(info.PerFrameFunctionalGroupsSequence.Item_1.CTPositionSequence.Item_1.DataCollectionCenterPatient(3)-info.PerFrameFunctionalGroupsSequence.Item_2.CTPositionSequence.Item_1.DataCollectionCenterPatient(3));
        slice_thickness = info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
        kVp = info.SharedFunctionalGroupsSequence.Item_1.CTXRayDetailsSequence.Item_1.KVP;
        mAs = info.PerFrameFunctionalGroupsSequence.Item_1.CTExposureSequence.Item_1.ExposureInmAs;
        CTDIvol = info.PerFrameFunctionalGroupsSequence.Item_1.CTExposureSequence.Item_1.CTDIvol;
        kernel = info.SharedFunctionalGroupsSequence.Item_1.CTReconstructionSequence.Item_1.ConvolutionKernel;
        if kernel == "FC12" || kernel == "CARDIAC"
            sort = 1;
        elseif kernel == "FC11"
            sort = 2;
        elseif kernel == "FC17"
            sort = 3;
        else
            sort = 99;
        end
        IR = info.SharedFunctionalGroupsSequence.Item_1.CTReconstructionSequence.Item_1.ReconstructionAlgorithm;
        if contains(IR,'VOLUMEXACT')
            IR = info.SharedFunctionalGroupsSequence.Item_1.CTReconstructionSequence.Item_1.ImageFilter;
            if contains(IR,'ORG')
                level = 0;
            elseif contains(IR,'AIDR 3D MILD')
                level = 1;
            elseif contains(IR,'AIDR 3D eSTD')
                level = 4;
            elseif contains(IR,'AIDR 3D STD')
                level = 2;
            elseif contains(IR,'AIDR 3D STR')
                level = 3;
            end
        elseif contains(IR,'FIRST')
            level = 11;
        elseif contains(IR,'AICE')
            level = 21;
        else
            level = 99;
        end
        pixelsize = info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1);
        fn = info.Filename;                                     %filename
        mf = info.Manufacturer;                                 %manufacturer name
        description = info.SeriesDescription;                   %description of the scan
    end
    
    try
        [V,~,~] = dicomreadVolume(scans);                       %load the dicom files into 4D volumes
    catch
        newfolder = strcat(scans,'\dfile');
        mkdir(newfolder)
        file = strcat(scans,'\DIRFILE');
        movefile(file,newfolder);
        [V,~,~] = dicomreadVolume(scans);                       %load the dicom files into 4D volumes
    end
    threed = double(squeeze(V));                            %remove singleton dimension and than make a double instead of int for the interpolation
    try
        threed = threed*info.RescaleSlope+info.RescaleIntercept;%multiply the volume with the slope and subtract the intercept to get the window level correct
    catch
        threed = threed*info.SharedFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleSlope+info.SharedFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1.RescaleIntercept;
    end
    [x,y,z] = size(threed);                                 %calculate the size of the 3 dimensions of the CT image
    
    if strcmp(method,'vendorspecific')
        if contains(mf,'GE MEDICAL SYSTEMS')
            method = 'GE';
        elseif contains(mf,'Philips')
            method = 'Philips';
        elseif contains(mf,'SIEMENS')
            method = 'Siemens';
        elseif contains(mf,'TOSHIBA') || contains(mf,'Canon')
            method = 'Canon';
        end
    end
    model = info.ManufacturerModelName;                     %CT system model name
    
%Sometimes the difference between slice locations is larger than the slice
%thickness, which should not be possible. Therefore different slice
%locations are used.
    n = 1;
    while abs(st) > slice_thickness
        info = dicominfo(fullfile(scans,dicomfiles(4+n).name));
        info2 = dicominfo(fullfile(scans,dicomfiles(5+n).name));
        n = n + 1;
        st = abs(info.SliceLocation - info2.SliceLocation);
    end
    
    if contains(mf,'Imatron')
        time = info.ContentTime;
    else
        time = info.SeriesTime;                             %time when the CT scan was taken
    end
    time = str2double(time);
    
    [columnsinimage, rowsinimage] = meshgrid(1:x,1:y);      %create a meshgrid with the same size as the image
    circlemask = zeros(x,y);
    center_calibration = [];
    center_slices = [];
    
%The following loop looks for the large calibration rod in every slice and
%saves the slices and center of the rod, i.e. the center of the cardiac
%insert. When the rod is not found anymore, it defines the slices where the
%calcifications are positioned. It thereby takes into account how the
%phantom is positioned.
    for calibration_slice = 1:z
        im = uint8(255.*mat2gray(threed(:,:,calibration_slice),[100,300])); %make a grey image of the slice with the calibration rods
        J = imgaussfilt(im,4);                              %smoothen the image so imfindcircles can find the circles easier
        [center,~] = imfindcircles(J,[floor(5/pixelsize),ceil(12.5/pixelsize)]);    %find the center and the radius of the circle in the slice, which has a radius between 7.5 and 12.5 mm, to only get the circle with the calcifications
        
        if isempty(center) && isempty(center_calibration)
            continue
        elseif (isempty(center) || center(1,2) >= 275 || center(1,2) <= 150 || center(1,1) >= 300 || center(1,1) <= 200) && ~isempty(center_calibration) && size_cs(1) <= 3
            center_slices = [];
            center_calibration = [];
            continue
        elseif (isempty(center) || center(1,2) >= 275 || center(1,2) <= 150 || center(1,1) >= 300 || center(1,1) <= 200) && ~isempty(center_calibration) && size_cs(1) > 3
            if calibration_slice <= (z/2)
                circlemask((rowsinimage-center_slices(end-1,3)).^2+(columnsinimage-center_slices(end-1,2)).^2 <= (46.5/pixelsize).^2)=1;       %(x-a)^2+(y-b)^2=r^2 is the formula of a circle, where (a,b) is the center of the circle. The mask sets all the pixels inside the calcification circle on 1 and all other pixels on 0
                calcslices = threed(:,:,(round(calibration_slice+round(3*3/st))):(round(calibration_slice+round(14*3/st))));
                homogeneous_slicenr = calibration_slice+round(17*3/st);
                center_hs = [center_slices(end-1,2), center_slices(end-1,3)];
            else
                circlemask((rowsinimage-center_slices(2,3)).^2+(columnsinimage-center_slices(2,2)).^2 <= (46.5/pixelsize).^2)=1;       %(x-a)^2+(y-b)^2=r^2 is the formula of a circle, where (a,b) is the center of the circle. The mask sets all the pixels inside the calcification circle on 1 and all other pixels on 0
                calcslices = threed(:,:,(round(center_slices(2,1)-round(14*3/st))):(round(center_slices(2,1)-round(3*3/st))));
                homogeneous_slicenr = center_slices(1,1)-round(17*3/st);
                center_hs = [center_slices(2,2), center_slices(2,3)];
            end
            sc = size(calcslices);
            circlemask = repmat(circlemask,[1,1,sc(3)]);
            calccircle = circlemask.*calcslices;
            homogeneous_slice = threed(:,:,homogeneous_slicenr);
            break
        elseif center(1,2) <= 275 && center(1,2) >= 175 && center(1,1) <= 300 && center(1,1) >= 200     %to prevent errors because of circles localized on a different place than the real calibration rod, the center has to be positioned y-coordinate 150 and 275
            center_calibration = center;
            center_slices = [center_slices; calibration_slice, center_calibration(1,1), center_calibration(1,2)]; %#ok<*AGROW>
            size_cs = size(center_slices);
        end
    end
    
%For image quality measurements the center slice of the calibration rods is
%used.
    nr_center_slices = size(center_slices);
    center_calibration = [center_slices(round(nr_center_slices(1)/2),2) center_slices(round(nr_center_slices(1)/2),3)];
    centerslice_calibration = threed(:,:,center_slices(round(nr_center_slices(1)/2),1));
    
    lesion_info_Ag = [];
    leslocs = [];

    Agatston = 0;                                           %Total Agatston Score
    Volume = 0;                                             %Total Volume score
    Mass = 0;                                               %Total Mass score

    if strcmp(method,'literatureCAC') || strcmp(method,'GE') || strcmp(method,'TeraRecon')
        numpixthreshold = ceil(1/(pixelsize^2));            %a structure needs to be larger than 1 mm2 to be defined as a calcification
    elseif strcmp(method,'Siemens')
        numpixthreshold = 0;
    elseif strcmp(method,'Philips')
        numpixthreshold = ceil(0.5/(pixelsize^2));
    elseif strcmp(method,'Canon')
        numpixthreshold = 3;
    end
    
%     BW = imbinarize(calccircle,129);                        %make binary image of the 3D CT images. Everything above threshold 129 is 1, everything below is 0

    if kVp == 120
        BW = imbinarize(calccircle,129);                        %make binary image of the 3D CT images. Everything above threshold 129 is 1, everything below is 0
    elseif kVp == 100
        BW = imbinarize(calccircle,146);                        %100 kVp threshold
    end

%In this double for loop the script goes through every slice where the
%calcifications could be positioned and calculates the calcium scores, the
%position, the mean HU value, and the distance and angle to the center of
%every structure.
    for slice = 1:sc(3)
        CC = bwconncomp(BW(:,:,slice),4);                   %find all connected components in the slice, i.e. structures. The last number is the connectivity
        numPixels = cellfun(@numel, CC.PixelIdxList);       %count the amount of structures in the slice
        small = numPixels < numpixthreshold;                %find all structures smaller than 1 mm2
        CC.PixelIdxList(small) = [];                        %remove all the small structures
        numobjects = size(CC.PixelIdxList);                 %amount of structures in the slice
        sl = calcslices(:,:,slice);
        
        for lesion = 1:numobjects(2)
            les = CC.PixelIdxList{lesion};
            numvox = size(les);                             %amount of pixels in the lesion
            Ai = numvox(1)*pixelsize^2;                     %area of the lesion in mm2
            HU = sl(les);                                   %HU values in the lesion
            maxHU = max(HU);                                %maximum HU value in the lesion
            CTi = mean(HU);                                 %mean HU value of the calcification
            
            %calculate distance from each lesion to center of the cardiac
            %insert.
            [yy,xx] = ind2sub([x,y],les);
            middlex = (min(xx)+max(xx))/2;
            middley = (min(yy)+max(yy))/2;
            centerdistance = sqrt(abs((middlex-center_calibration(1,1))^2 + (middley-center_calibration(1,2))^2));
            
            %This if statement will define the weighting factor for the
            %Agatston score.
            if maxHU >= 130 && maxHU < 200
                Wi = 1;
            elseif maxHU >= 200 && maxHU < 300
                Wi = 2;
            elseif maxHU >= 300 && maxHU < 400
                Wi = 3;
            elseif maxHU >= 400
                Wi = 4;
            else
                Wi = 0;
%                 fprintf('Error, maxHU has a strange value, look into it: maxHU = %3.0f, [X,Y,Z] = [%3.0f, %3.0f, %3.0f], lesion = %2.0f\n', maxHU, middlex, middley, slice, lesion)
            end

            CSi = Wi*Ai;                                    %the Agatston score in the lesion
            Agatston = Agatston + CSi;                      %the total Agatston score at that moment

            Vi = numvox(1)*pixelsize^2*st;
            Volume = Volume + Vi;                           %the Volume score (is not used for our research, see volumescore_thesis.m)
            mi = Vi*CTi;                    %watch out! This is not the complete mass score, the calibration factor will be calculated later on in the script
            Mass = Mass + mi;
            
            location = [les, {"extra"}];
            leslocs = [leslocs; location];                  %the locations are saved for troubleshooting
            tdat = [slice, middlex, middley, CTi, CSi, Agatston, Ai, Vi, Volume, mi, Mass, centerdistance, 0, Wi];
            lesion_info_Ag = [lesion_info_Ag; tdat];
        end
    end
    
    r = size(lesion_info_Ag);
    slices = lesion_info_Ag(:,1);
    [amount,numbers] = hist(slices,unique(slices));
    label = 1;
    am = 0;
    w = 0;
    
%This for loop labels every calcification. The calcifications in the first
%slice do not have any calcifications to compare with, so they all get a
%unique label. 'Else' it gets to the following slices, it will compare
%with all calcifications in the previous slice. If any pixel is overlapping
%with a pixel in the calcification in the previous slice, it gets the same
%label and it won't compare with the other calcifications. If it doesn't
%overlap with calcifications in the previous slice, it's a new
%calcification and it gets a new label.
    for s = 1:r(1)
        v = amount(numbers==(lesion_info_Ag(s,1)));         %v is the amount of calcifications in the slice
        if lesion_info_Ag(s,1) == numbers(1)
            lesion_info_Ag(s,13) = label;
            label = label + 1;
        else
            overlap = amount(numbers==(lesion_info_Ag(s-w-1,1)));
            for u = 1:overlap
                if any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1})) && lesion_info_Ag(s,13) == 0
                    lesion_info_Ag(s,13) = lesion_info_Ag((s-(overlap+am)+(u-1)),13);
                elseif any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1}))
                    lesion_info_Ag(lesion_info_Ag(:,13) == lesion_info_Ag((s-(overlap+am)+(u-1)),13),13) = lesion_info_Ag(s,13);      %calcifications in noisy images can fall apart, so to be sure every "island" is included, this sentence gives all these "islands" the same label
                elseif u == overlap && lesion_info_Ag(s,13) == 0
                    lesion_info_Ag(s,13) = label;
                    label = label + 1;
                end
            end
        end
        w = w + 1;
        if w == v
            am = 0;
            w = 0;
        else
            am = am + 1;
        end
    end
    
    labels = lesion_info_Ag(:,13);                          %the column with labels
    groupedData = arrayfun(@(y)find(labels == y), unique(labels), 'UniformOutput',false);   %find unique labels, i.e. unique calcifications
    nrofcalc = size(groupedData);
    sumdata_Ag = zeros(nrofcalc(1),12);
    
%This for loop adds up all scores and measurements with the same label,
%i.e. the same calcification, defined by the previous loop and puts it in a
%matrix.
    for add = 1:nrofcalc(1)
        calcification = lesion_info_Ag(groupedData{add,1},:);
        A = calcification(:,[2:12 14]);
        middlex = median(A(:,1));
        middley = median(A(:,2));
        meanHU = mean(A(:,3));
        Ag = sum(A(:,4))/(3/st);                            %because traditional Agatston scores were acquired on an EBT scanner with 3 mm slice thickness and without overlapping slices, the Agatston score has to be corrected for thinner or overlapping slices
        vol = sum(A(:,7));
        mass = sum(A(:,9));
        meanWi = mean(A(:,12));

        mediandistance = median(A(:,11));
        maxarea = max(A(:,6));                              %the maximum area is used here, because the partial volume effect can have a large influence on the size of the area in the outer slices of a calcification
        
        sumdata_Ag(add,2) = middlex;
        sumdata_Ag(add,3) = middley;
        sumdata_Ag(add,4) = meanHU;
        sumdata_Ag(add,5) = Ag;
        sumdata_Ag(add,6) = vol;
        sumdata_Ag(add,7) = mass;
        
        sumdata_Ag(add,9) = mediandistance;
        sumdata_Ag(add,10) = maxarea;
        sumdata_Ag(add,11) = calcification(1,13);
        sumdata_Ag(add,12) = meanWi;
    end
    
    sortarea = sortrows(sumdata_Ag, 10, 'descend');         %largest areas must be the largest calcifications in this phantom
    sortdens = sortrows(sortarea(1:3,:), 4, 'descend');     %sort on density to find the different densities
    sumdata_Ag(find(sumdata_Ag(:,11)==sortdens(1,11),1),1) = 1; %label large size, high density (HD) calcification
    sumdata_Ag(find(sumdata_Ag(:,11)==sortdens(2,11),1),1) = 2; %label large size, medium density (MD) calcification
    sumdata_Ag(find(sumdata_Ag(:,11)==sortdens(3,11),1),1) = 3; %label large size, low density (LD) calcification
    
    r1 = sortdens(1,9);                                     %distance to center of large size, high density calcification
    r2 = sortdens(2,9);                                     %distance to center of large size, medium density calcification
    r3 = sortdens(3,9);                                     %distance to center of large size, low density calcification
    ravg = (r1+r2+r3)/3;                                    %thus the distance to the "real center" should be the average of those three
    
    middlex1 = sortdens(1,2);                               %x value of center of large size, high density calcification
    middlex2 = sortdens(2,2);                               %x value of center of large size, medium density calcification
    middlex3 = sortdens(3,2);                               %x value of center of large size, low density calcification
    middley1 = sortdens(1,3);                               %y value of center of large size, high density calcification
    middley2 = sortdens(2,3);                               %y value of center of large size, medium density calcification
    middley3 = sortdens(3,3);                               %y value of center of large size, low density calcification
    circcenter1 = [middlex1;middley1];
    circcenter2 = [middlex2;middley2];
    
%Find the intersections of two circles with the correct radius with the
%large HD and MD as centerpoints. The real centerpoint is the intersection
%closest to the center of the large LD calcification.
    distancecenters = sum((circcenter2-circcenter1).^2);
    coordinate0 = (circcenter1+circcenter2)/2+(ravg^2-ravg^2)/distancecenters/2*(circcenter2-circcenter1);
    overlap = ((ravg+ravg)^2-distancecenters)*(distancecenters-(ravg-ravg)^2);
    if overlap <= 0
        fprintf('The circles don''t intersect.\n')
    else
        coordinatediff = sqrt(overlap)/distancecenters/2*[0 -1;1 0]*(circcenter2-circcenter1);
        intersect1 = coordinate0 + coordinatediff;
        intersect2 = coordinate0 - coordinatediff;
    end
    distanceLD1 = sqrt(abs((middlex3-intersect1(1))^2 + (middley3-intersect1(2))^2));
    distanceLD2 = sqrt(abs((middlex3-intersect2(1))^2 + (middley3-intersect2(2))^2));
    
    if distanceLD1 < distanceLD2
        center_calcifications = [intersect1(1),intersect1(2)];
    else
        center_calcifications = [intersect2(1),intersect2(2)];
    end
    
    for replace = 1:nrofcalc(1)
        deltax = sumdata_Ag(replace,2)-center_calcifications(1,1);
        deltay = sumdata_Ag(replace,3)-center_calcifications(1,2);
        
        %Because the origin (0,0) of the image is in the upper left corner
        %in MATLAB, a lesion above the center has a negative angle and
        %below the center has a positive angle.
        angle = atan2d(deltay, deltax);                     %measure the angle to the center of the calcifications
        centerdistance = sqrt(abs((sumdata_Ag(replace,2)-center_calcifications(1,1))^2 + (sumdata_Ag(replace,3)-center_calcifications(1,2))^2));        %measures the distance to the center of the calcifications
        sumdata_Ag(replace,8) = angle;
        sumdata_Ag(replace,9) = centerdistance;
    end
    
    highdensangle = sumdata_Ag(sumdata_Ag(:,1)==1,8);       %angle to the center of the large size, high density calcification
    highdensdistance = sumdata_Ag(sumdata_Ag(:,1)==1,9);    %distance to the center of the large size, high density calcification
    meddensangle = sumdata_Ag(sumdata_Ag(:,1)==2,8);        %angle to the center of the large size, medium density calcification
    meddensdistance = sumdata_Ag(sumdata_Ag(:,1)==2,9);     %distance to the center of the large size, medium density calcification
    lowdensangle = sumdata_Ag(sumdata_Ag(:,1)==3,8);        %angle to the center of the large size, low density calcification
    lowdensdistance = sumdata_Ag(sumdata_Ag(:,1)==3,9);     %distance to the center of the large size, low density calcification
    
%This for loop labels the calcification with label 4-9 depending on the
%position of the calcification with respect to the large calcifications.
    for calcorder = 1:nrofcalc(1)
        if (highdensangle-3) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (highdensangle+3) && (8/pixelsize) <= (highdensdistance - sumdata_Ag(calcorder,9)) && (highdensdistance - sumdata_Ag(calcorder,9)) <= (16/pixelsize)
            sumdata_Ag(calcorder,1) = 4;                    %medium size, high density calcification
        elseif (meddensangle-3) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (meddensangle+3) && (8/pixelsize) <= (meddensdistance - sumdata_Ag(calcorder,9)) && (meddensdistance - sumdata_Ag(calcorder,9)) <= (16/pixelsize)
            sumdata_Ag(calcorder,1) = 5;                    %medium size, medium density calcification
        elseif (lowdensangle-3) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (lowdensangle+3) && (8/pixelsize) <= (lowdensdistance - sumdata_Ag(calcorder,9)) && (lowdensdistance - sumdata_Ag(calcorder,9)) <= (16/pixelsize)
            sumdata_Ag(calcorder,1) = 6;                    %medium size, low density calcification
        elseif (highdensangle-5) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (highdensangle+5) && (21/pixelsize) <= (highdensdistance - sumdata_Ag(calcorder,9)) && (highdensdistance - sumdata_Ag(calcorder,9)) <= (28.5/pixelsize)
            sumdata_Ag(calcorder,1) = 7;                    %small size, high density calcification
        elseif (meddensangle-5) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (meddensangle+5) && (21/pixelsize) <= (meddensdistance - sumdata_Ag(calcorder,9)) && (meddensdistance - sumdata_Ag(calcorder,9)) <= (28.5/pixelsize)
            sumdata_Ag(calcorder,1) = 8;                    %small size, medium density calcification
        elseif (lowdensangle-5) <= sumdata_Ag(calcorder,8) && sumdata_Ag(calcorder,8) <= (lowdensangle+5) && (21/pixelsize) <= (lowdensdistance - sumdata_Ag(calcorder,9)) && (lowdensdistance - sumdata_Ag(calcorder,9)) <= (28.5/pixelsize)
            sumdata_Ag(calcorder,1) = 9;                    %small size, low density calcification
        end
    end
    
%Remove all slices before the first slice with a large calcification and
%after the last slice with a large calcification. Because there is a chance
%that some noise will be on exactly the same spot as a medium sized
%calcification in the slice before the first or after the last, there is a
%one slice extra range on both sides.
    large = find(ismember(labels, sortarea(1:3,11)));
    beginslice = lesion_info_Ag(large(1),1);
    endslice = lesion_info_Ag(large(end),1);
    slices = lesion_info_Ag(:,1);
    remove1 = find(slices<beginslice-1);
    remove2 = find(slices>endslice+1);
    if isempty(remove1) && isempty(remove2)
        
    elseif isempty(remove1)
        L2 = find(sumdata_Ag(:,11)==lesion_info_Ag(remove2(1),13),1);
        sumdata_Ag(L2:end,:) = [];
    elseif isempty(remove2)
        L1 = find(sumdata_Ag(:,11)==lesion_info_Ag(remove1(end),13),1);
        sumdata_Ag(1:L1,:) = [];
    else
        L1 = find(sumdata_Ag(:,11)==lesion_info_Ag(remove1(end),13),1);
        L2 = find(sumdata_Ag(:,11)==lesion_info_Ag(remove2(1),13),1);
        sumdata_Ag(L2:end,:) = [];
        sumdata_Ag(1:L1,:) = [];
    end
    
    noise = sumdata_Ag(:,1) == 0;                           %find all noise in the slices
    nrofnoise = size(noise);                                %nr of false positives, i.e. noise
    sumdata_Ag(noise,:) = [];                               %remove all noise from sumdata
    
    sumdata_Ag = sortrows(sumdata_Ag,1);                    %sort sumdata in order of calcification number
    
%If any large or medium sized calcification is not found or the max area of
%all calcifications is smaller than 2 mm2 and thus probably only noise, it
%will give an error and you should look into it.
    large = [1,2,3];
    medium = [4,5,6];
    if any(~ismember(large,sumdata_Ag(:,1)))
        fprintf('Error, one or more calcifications of large size have been removed accidentally in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
    elseif any(~ismember(medium,sumdata_Ag(:,1)))
        fprintf('Error, one or more calcifications of medium size are not found in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
    elseif all(sumdata_Ag(sumdata_Ag(:,1)==medium(1),10)<2)
        fprintf('Error, the medium size, high density calcification (nr4) is not properly found in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
    elseif all(sumdata_Ag(sumdata_Ag(:,1)==medium(2),10)<2)
        fprintf('Error, the medium size, medium density calcification (nr5) is not properly found in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
    elseif all(sumdata_Ag(sumdata_Ag(:,1)==medium(3),10)<2)
        fprintf('Error, the medium size, low density calcification (nr6) is not properly found in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
    end
    
%Image quality and Mass score calibration factor measurements
    calibration = zeros(x,y);
    water = zeros(x,y);
    calibration((rowsinimage-center_calibration(1,2)).^2+(columnsinimage-center_calibration(1,1)).^2 <= (6.9/pixelsize).^2)=1;        %use the center of the calccircle to make an ROI for the calibration rod in the middle. I used a radius to make an ROI of 1.5 cm2 to make sure the ROI will keep inside the rod
    calcrod = centerslice_calibration(logical(calibration));
    CTc = mean(calcrod);
    
    angle_to_water = lowdensangle-180;                      %the water-equivalent rod is placed exactly 180 degrees from the low density calcifications
    x_water = (22.5/pixelsize)*cosd(angle_to_water);        %the water-equivalent rod is placed at 22.5 mm from the center of the cardiac insert
    y_water = (22.5/pixelsize)*sind(angle_to_water);
    center_waterrod = [(center_calibration(1,1)+x_water), (center_calibration(1,2)+y_water)];   %saves the center of the water-equivalent rod
    water((rowsinimage-center_waterrod(1,2)).^2+(columnsinimage-center_waterrod(1,1)).^2 <= (6.9/pixelsize).^2)=1;        %use the center of the calccircle to make an ROI for the calibration rod in the middle. I used a radius to make an ROI of 1.5 cm2 to make sure the ROI will keep well within distance of the edges of the rod
    waterrod = centerslice_calibration(logical(water));     %all pixel values within the ROI
    CTw = mean(waterrod);                                   %the mean HU value in the ROI
    SD_noise = std(waterrod);                               %the standard deviation in the ROI
    
%     %if you want to see how the ROIs are placed, run the following:
%     imshow(centerslice_calibration,[-50,100])
%     h = viscircles([center_waterrod; center_calibration],[(6.9/pixelsize);(6.9/pixelsize)],'LineWidth',1);
%     h.Children(1).Color = 'r';
%     h.Children(2).Color = 'r';
% %     print(gcf,'ROIs_calibration_rods','-dtiff','-r600')
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Volumescore with interpolation
    if strcmp(method,'literatureCAC') %|| strcmp(method,'Siemens') || strcmp(method,'TeraRecon')
        Zint = pixelsize/st;
        [X,Y,Z] = meshgrid(1:x,1:y,1:sc(3));                            %defines the amount of X, Y and Z values for the original dicomfiles
        [Xq,Yq,Zq] = meshgrid(1:x,1:y,1:Zint:sc(3));                    %defines the amount of X, Y and Z values you want after interpolating
        calccircle_int = interp3(X,Y,Z,calccircle,Xq,Yq,Zq);
%         [X,Y,Z] = meshgrid(1:x,1:y,1:z);                            %defines the amount of X, Y and Z values for the original dicomfiles
%         [Xq,Yq,Zq] = meshgrid(1:x,1:y,1:Zint:z);                    %defines the amount of X, Y and Z values you want after interpolating
%         Vint = interp3(X,Y,Z,threed,Xq,Yq,Zq);
%         circlemask = zeros(x,y);
%         if calibration_slice <= (z/2)
%             circlemask((rowsinimage-center_slices(end-1,3)).^2+(columnsinimage-center_slices(end-1,2)).^2 <= (46.5/pixelsize).^2)=1;       %(x-a)^2+(y-b)^2=r^2 is the formula of a circle, where (a,b) is the center of the circle. The mask sets all the pixels inside the calcification circle on 1 and all other pixels on 0
%             calcslices_int = Vint(:,:,(round(calibration_slice*3/pixelsize+round(3*3/pixelsize))):(round(calibration_slice*3/pixelsize+round(14*3/pixelsize))));
%         else
%             circlemask((rowsinimage-center_slices(2,3)).^2+(columnsinimage-center_slices(2,2)).^2 <= (46.5/pixelsize).^2)=1;       %(x-a)^2+(y-b)^2=r^2 is the formula of a circle, where (a,b) is the center of the circle. The mask sets all the pixels inside the calcification circle on 1 and all other pixels on 0
%             calcslices_int = Vint(:,:,(round(center_slices(2,1)*3/pixelsize-round(14*3/pixelsize))):(round(center_slices(2,1)*3/pixelsize-round(3*3/pixelsize))));
%         end
%         sc = size(calcslices_int);
%         circlemask = repmat(circlemask,[1,1,sc(3)]);
%         calccircle_int = circlemask.*calcslices_int;

        lesion_info_Vol = [];
        leslocs = [];
        Volume = 0;                                             %Total Volume score
        
        if strcmp(method,'literatureCAC') || strcmp(method,'TeraRecon')
            numpixthreshold = ceil(1/(pixelsize^2));            %a structure needs to be larger than 1 mm2 to be defined as a calcification
        elseif strcmp(method,'Siemens')
            numpixthreshold = 0;
        end
        
        BW = imbinarize(calccircle_int,129);                        %make binary image of the 3D CT images. Everything above threshold 129 is 1, everything below is 0

%         if kVp == 120
%             BW = imbinarize(calccircle_int,129);                        %make binary image of the 3D CT images. Everything above threshold 129 is 1, everything below is 0
%         elseif kVp == 100
%             BW = imbinarize(calccircle_int,146);                        %100 kVp threshold
%         end
        [~,~,z_int] = size(calccircle_int);
        
        %In this double for loop the script goes through every slice where the
        %calcifications could be positioned and calculates the calcium scores, the
        %position, the mean HU value, and the distance and angle to the center of
        %every structure.
        for slice = 1:z_int
            CC = bwconncomp(BW(:,:,slice),4);                   %find all connected components in the slice, i.e. structures
            numPixels = cellfun(@numel, CC.PixelIdxList);       %count the amount of structures in the slice
            small = numPixels < numpixthreshold;                %find all structures smaller than 1 mm2
            CC.PixelIdxList(small) = [];                        %remove all the small structures
            numobjects = size(CC.PixelIdxList);                 %amount of structures in the slice
            sl = calccircle_int(:,:,slice);
            
            for lesion = 1:numobjects(2)
                les = CC.PixelIdxList{lesion};
                numvox = size(les);                             %amount of pixels in the lesion
                Ai = numvox(1)*pixelsize^2;                     %area of the lesion in mm2
                HU = sl(les);                                   %HU values in the lesion
                CTi = mean(HU);                                 %mean HU value of the calcification
                
                %calculate distance from each lesion to center of the cardiac
                %insert.
                [yy,xx] = ind2sub([x,y],les);
                middlex = (min(xx)+max(xx))/2;
                middley = (min(yy)+max(yy))/2;
                centerdistance = sqrt(abs((middlex-center_calcifications(1,1))^2 + (middley-center_calcifications(1,2))^2));
                deltax = middlex-center_calcifications(1,1);
                deltay = middley-center_calcifications(1,2);
                angle = atan2d(deltay, deltax);                     %measure the angle to the center of the calcifications

                Vi = numvox(1)*pixelsize^3;
                Volume = Volume + Vi;
                
                location = [les, {"extra"}];
                leslocs = [leslocs; location];                  %the locations are saved for troubleshooting
                tdat = [slice, middlex, middley, CTi, Ai, Vi, Volume, centerdistance, 0, angle];
                lesion_info_Vol = [lesion_info_Vol; tdat];
            end
        end
        
        r = size(lesion_info_Vol);
        slices = lesion_info_Vol(:,1);
        [amount,numbers] = hist(slices,unique(slices));
        label = 1;
        am = 0;
        w = 0;
        
        %This for loop labels every calcification. The calcifications in the first
        %slice do not have any calcifications to compare with, so they all get a
        %unique label. 'Else' it gets to the following slices, it will compare
        %with all calcifications in the previous slice. If any pixel is overlapping
        %with a pixel in the calcification in the previous slice, it gets the same
        %label and it won't compare with the other calcifications. If it doesn't
        %overlap with calcifications in the previous slice, it's a new
        %calcification and it gets a new label.
        for s = 1:r(1)
            v = amount(numbers==(lesion_info_Vol(s,1)));         %v is the amount of calcifications in the slice
            if lesion_info_Vol(s,1) == numbers(1)
                lesion_info_Vol(s,9) = label;
                label = label + 1;
            else
                overlap = amount(numbers==(lesion_info_Vol(s-w-1,1)));
                for u = 1:overlap
                    if any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1})) && lesion_info_Vol(s,9) == 0
                        lesion_info_Vol(s,9) = lesion_info_Vol((s-(overlap+am)+(u-1)),9);
                    elseif any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1}))
                        lesion_info_Vol(lesion_info_Vol(:,9) == lesion_info_Vol((s-(overlap+am)+(u-1)),9),9) = lesion_info_Vol(s,9);      %calcifications in noisy images can fall apart, so to be sure every "island" is included, this sentence gives all these "islands" the same label
                    elseif u == overlap && lesion_info_Vol(s,9) == 0
                        lesion_info_Vol(s,9) = label;
                        label = label + 1;
                    end
                end
            end
            w = w + 1;
            if w == v
                am = 0;
                w = 0;
            else
                am = am + 1;
            end
        end
        
        labels = lesion_info_Vol(:,9);                          %the column with labels
        groupedData = arrayfun(@(y)find(labels == y), unique(labels), 'UniformOutput',false);   %find unique labels, i.e. unique calcifications
        nrofcalc = size(groupedData);
        sumdata_Vol = zeros(nrofcalc(1),9);
        
        %This for loop adds up all scores and measurements with the same label,
        %i.e. the same calcification, defined by the previous loop and puts it in a
        %matrix.
        for add = 1:nrofcalc(1)
            calcification = lesion_info_Vol(groupedData{add,1},:);
            A = calcification(:,[2:8 10]);
            middlex = median(A(:,1));
            middley = median(A(:,2));
            meanHU = mean(A(:,3));
            maxarea = max(A(:,4));
            vol = sum(A(:,5));
            medianangle = median(A(:,8));
            mediandistance = median(A(:,7));
            
            sumdata_Vol(add,2) = middlex;
            sumdata_Vol(add,3) = middley;
            sumdata_Vol(add,4) = meanHU;
            sumdata_Vol(add,5) = vol;
            sumdata_Vol(add,6) = medianangle;
            sumdata_Vol(add,7) = mediandistance;
            sumdata_Vol(add,8) = maxarea;
            sumdata_Vol(add,9) = calcification(1,9);
        end
        
        sortarea = sortrows(sumdata_Vol, 8, 'descend');         %largest areas must be the largest calcifications in this phantom
        sortdens = sortrows(sortarea(1:3,:), 4, 'descend');     %sort on density to find the different densities
        sumdata_Vol(find(sumdata_Vol(:,9)==sortdens(1,9),1),1) = 1; %label large size, high density (HD) calcification
        sumdata_Vol(find(sumdata_Vol(:,9)==sortdens(2,9),1),1) = 2; %label large size, medium density (MD) calcification
        sumdata_Vol(find(sumdata_Vol(:,9)==sortdens(3,9),1),1) = 3; %label large size, low density (LD) calcification
        
        highdensangle = sumdata_Vol(sumdata_Vol(:,1)==1,6);       %angle to the center of the large size, high density calcification
        highdensdistance = sumdata_Vol(sumdata_Vol(:,1)==1,7);    %distance to the center of the large size, high density calcification
        meddensangle = sumdata_Vol(sumdata_Vol(:,1)==2,6);        %angle to the center of the large size, medium density calcification
        meddensdistance = sumdata_Vol(sumdata_Vol(:,1)==2,7);     %distance to the center of the large size, medium density calcification
        lowdensangle = sumdata_Vol(sumdata_Vol(:,1)==3,6);        %angle to the center of the large size, low density calcification
        lowdensdistance = sumdata_Vol(sumdata_Vol(:,1)==3,7);     %distance to the center of the large size, low density calcification
        
        %This for loop labels the calcification with label 4-9 depending on the
        %position of the calcification with respect to the large calcifications.
        for calcorder = 1:nrofcalc(1)
            if (highdensangle-3) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (highdensangle+3) && (8/pixelsize) <= (highdensdistance - sumdata_Vol(calcorder,7)) && (highdensdistance - sumdata_Vol(calcorder,7)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 4;                    %medium size, high density calcification
            elseif (meddensangle-3) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (meddensangle+3) && (8/pixelsize) <= (meddensdistance - sumdata_Vol(calcorder,7)) && (meddensdistance - sumdata_Vol(calcorder,7)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 5;                    %medium size, medium density calcification
            elseif (lowdensangle-3) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (lowdensangle+3) && (8/pixelsize) <= (lowdensdistance - sumdata_Vol(calcorder,7)) && (lowdensdistance - sumdata_Vol(calcorder,7)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 6;                    %medium size, low density calcification
            elseif (highdensangle-5) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (highdensangle+5) && (21/pixelsize) <= (highdensdistance - sumdata_Vol(calcorder,7)) && (highdensdistance - sumdata_Vol(calcorder,7)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 7;                    %small size, high density calcification
            elseif (meddensangle-5) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (meddensangle+5) && (21/pixelsize) <= (meddensdistance - sumdata_Vol(calcorder,7)) && (meddensdistance - sumdata_Vol(calcorder,7)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 8;                    %small size, medium density calcification
            elseif (lowdensangle-5) <= sumdata_Vol(calcorder,6) && sumdata_Vol(calcorder,6) <= (lowdensangle+5) && (21/pixelsize) <= (lowdensdistance - sumdata_Vol(calcorder,7)) && (lowdensdistance - sumdata_Vol(calcorder,7)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 9;                    %small size, low density calcification
            end
        end
        
        %Remove all slices before the first slice with a large calcification and
        %after the last slice with a large calcification. Because there is a chance
        %that some noise will be on exactly the same spot as a medium sized
        %calcification in the slice before the first or after the last, there is a
        %one slice extra range on both sides.
        large = find(ismember(labels, sortarea(1:3,9)));
        beginslice = lesion_info_Vol(large(1),1);
        endslice = lesion_info_Vol(large(end),1);
        slices = lesion_info_Vol(:,1);
        remove1 = find(slices<beginslice-1);
        remove2 = find(slices>endslice+1);
        if isempty(remove1) && isempty(remove2)
            
        elseif isempty(remove1)
            L2 = find(sumdata_Vol(:,9)==lesion_info_Vol(remove2(1),9),1);
            sumdata_Vol(L2:end,:) = [];
        elseif isempty(remove2)
            L1 = find(sumdata_Vol(:,9)==lesion_info_Vol(remove1(end),9),1);
            sumdata_Vol(1:L1,:) = [];
        else
            L1 = find(sumdata_Vol(:,9)==lesion_info_Vol(remove1(end),9),1);
            L2 = find(sumdata_Vol(:,9)==lesion_info_Vol(remove2(1),9),1);
            sumdata_Vol(L2:end,:) = [];
            sumdata_Vol(1:L1,:) = [];
        end
        noise = sumdata_Vol(:,1) == 0;                            %find all noise in the slices
        sumdata_Vol(noise,:) = [];                                %remove all noise from sumdata
        sumdata_Vol = sortrows(sumdata_Vol,1);                    %sort sumdata in order of calcification number
        
        large = [1,2,3];
        medium = [4,5,6];
        if any(~ismember(large,sumdata_Vol(:,1)))
            fprintf('Error, one or more calcifications of large size have been removed accidentally at the interpolation part in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
        elseif any(~ismember(medium,sumdata_Vol(:,1)))
            fprintf('Error, one or more calcifications of medium size are not found at the interpolation part in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
        elseif all(sumdata_Vol(sumdata_Vol(:,1)==medium(1),8)<2)
            fprintf('Error, the medium size, high density calcification (nr4) is not properly found at the interpolation part in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
        elseif all(sumdata_Vol(sumdata_Vol(:,1)==medium(2),8)<2)
            fprintf('Error, the medium size, medium density calcification (nr5) is not properly found at the interpolation part in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
        elseif all(sumdata_Vol(sumdata_Vol(:,1)==medium(3),8)<2)
            fprintf('Error, the medium size, low density calcification (nr6) is not properly found at the interpolation part in the following scan: %3.0f, %3.0f, %1.1f, %1.0f (%s)\n', kVp, mAs, slice_thickness, level, scans)
        end
        
        size_Ag = size(sumdata_Ag);
        size_Vol = size(sumdata_Vol);
        matrix_diff = size_Vol(1) - size_Ag(1);
        if matrix_diff > 0
            sumdata_Ag(size_Ag(1)+matrix_diff,:) = zeros;
        elseif matrix_diff < 0
            sumdata_Vol(size_Vol(1)-matrix_diff,:) = zeros;
        end
        sumdata_Ag(:,6) = sumdata_Vol(:,5);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume and Mass score with 100 mgHA as threshold
    elseif strcmp(method,'Philips')
        lesion_info_Vol = [];
        leslocs = [];
        Volume = 0;                                             %Total Volume score
        Mass = 0;
        
        numpixthreshold = ceil(0.5/(pixelsize^2));
        BW = imbinarize(calccircle,(99/0.88));                  %threshold for Philips volume and mass score
        for slice = 1:sc(3)
            CC = bwconncomp(BW(:,:,slice),4);                   %find all connected components in the slice, i.e. structures
            numPixels = cellfun(@numel, CC.PixelIdxList);       %count the amount of structures in the slice
            small = numPixels < numpixthreshold;                %find all structures smaller than 1 mm2
            CC.PixelIdxList(small) = [];                        %remove all the small structures
            numobjects = size(CC.PixelIdxList);                 %amount of structures in the slice
            sl = calcslices(:,:,slice);
            
            for lesion = 1:numobjects(2)
                les = CC.PixelIdxList{lesion};
                numvox = size(les);                             %amount of pixels in the lesion
                Ai = numvox(1)*pixelsize^2;                     %area of the lesion in mm2
                HU = sl(les);                                   %HU values in the lesion
                CTi = mean(HU);                                 %mean HU value of the calcification
                
                %calculate distance from each lesion to center of the cardiac
                %insert.
                [yy,xx] = ind2sub([x,y],les);
                middlex = (min(xx)+max(xx))/2;
                middley = (min(yy)+max(yy))/2;
                centerdistance = sqrt(abs((middlex-center_calibration(1,1))^2 + (middley-center_calibration(1,2))^2));
                
                Vi = numvox(1)*pixelsize^2*st;
                Volume = Volume + Vi;
                mi = Vi*CTi;                    %watch out! This is not the complete mass score, the calibration factor will be calculated later on in the script
                Mass = Mass + mi;
                
                location = [les, {"extra"}];
                leslocs = [leslocs; location];                  %the locations are saved for troubleshooting
                tdat = [slice, middlex, middley, CTi, Ai, Vi, Volume, mi, Mass, centerdistance, 0];
                lesion_info_Vol = [lesion_info_Vol; tdat];
            end
        end
        
        r = size(lesion_info_Vol);
        slices = lesion_info_Vol(:,1);
        [amount,numbers] = hist(slices,unique(slices));
        label = 1;
        am = 0;
        w = 0;
        
        %This for loop labels every calcification. The calcifications in the first
        %slice do not have any calcifications to compare with, so they all get a
        %unique label. 'Else' it gets to the following slices, it will compare
        %with all calcifications in the previous slice. If any pixel is overlapping
        %with a pixel in the calcification in the previous slice, it gets the same
        %label and it won't compare with the other calcifications. If it doesn't
        %overlap with calcifications in the previous slice, it's a new
        %calcification and it gets a new label.
        for s = 1:r(1)
            v = amount(numbers==(lesion_info_Vol(s,1)));         %v is the amount of calcifications in the slice
            if lesion_info_Vol(s,1) == numbers(1)
                lesion_info_Vol(s,11) = label;
                label = label + 1;
            else
                overlap = amount(numbers==(lesion_info_Vol(s-w-1,1)));
                for u = 1:overlap
                    if any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1})) && lesion_info_Vol(s,11) == 0
                        lesion_info_Vol(s,11) = lesion_info_Vol((s-(overlap+am)+(u-1)),11);
                    elseif any(ismember(leslocs{s,1},leslocs{(s-(overlap+am)+(u-1)),1}))
                        lesion_info_Vol(lesion_info_Vol(:,11) == lesion_info_Vol((s-(overlap+am)+(u-1)),11),11) = lesion_info_Vol(s,11);      %calcifications in noisy images can fall apart, so to be sure every "island" is included, this sentence gives all these "islands" the same label
                    elseif u == overlap && lesion_info_Vol(s,11) == 0
                        lesion_info_Vol(s,11) = label;
                        label = label + 1;
                    end
                end
            end
            w = w + 1;
            if w == v
                am = 0;
                w = 0;
            else
                am = am + 1;
            end
        end
        
        labels = lesion_info_Vol(:,11);                          %the column with labels
        groupedData = arrayfun(@(y)find(labels == y), unique(labels), 'UniformOutput',false);   %find unique labels, i.e. unique calcifications
        nrofcalc = size(groupedData);
        sumdata_Vol = zeros(nrofcalc(1),10);
        
        %This for loop adds up all scores and measurements with the same label,
        %i.e. the same calcification, defined by the previous loop and puts it in a
        %matrix.
        for add = 1:nrofcalc(1)
            calcification = lesion_info_Vol(groupedData{add,1},:);
            A = calcification(:,2:10);
            middlex = median(A(:,1));
            middley = median(A(:,2));
            meanHU = mean(A(:,3));
            maxarea = max(A(:,4));
            vol = sum(A(:,5));
            mass = sum(A(:,7));
            mediandistance = median(A(:,9));
            
            sumdata_Vol(add,2) = middlex;
            sumdata_Vol(add,3) = middley;
            sumdata_Vol(add,4) = meanHU;
            sumdata_Vol(add,5) = vol;
            sumdata_Vol(add,6) = mass;
            
            sumdata_Vol(add,8) = mediandistance;
            sumdata_Vol(add,9) = maxarea;
            sumdata_Vol(add,10) = calcification(1,11);
        end
        
        sortarea = sortrows(sumdata_Vol, 9, 'descend');         %largest areas must be the largest calcifications in this phantom
        sortdens = sortrows(sortarea(1:3,:), 4, 'descend');     %sort on density to find the different densities
        sumdata_Vol(find(sumdata_Vol(:,10)==sortdens(1,10),1),1) = 1; %label large size, high density (HD) calcification
        sumdata_Vol(find(sumdata_Vol(:,10)==sortdens(2,10),1),1) = 2; %label large size, medium density (MD) calcification
        sumdata_Vol(find(sumdata_Vol(:,10)==sortdens(3,10),1),1) = 3; %label large size, low density (LD) calcification
        
        r1 = sortdens(1,8);                                     %distance to center of large size, high density calcification
        r2 = sortdens(2,8);                                     %distance to center of large size, medium density calcification
        r3 = sortdens(3,8);                                     %distance to center of large size, low density calcification
        ravg = (r1+r2+r3)/3;                                    %thus the distance to the "real center" should be the average of those three
        
        middlex1 = sortdens(1,2);                               %x value of center of large size, high density calcification
        middlex2 = sortdens(2,2);                               %x value of center of large size, medium density calcification
        middlex3 = sortdens(3,2);                               %x value of center of large size, low density calcification
        middley1 = sortdens(1,3);                               %y value of center of large size, high density calcification
        middley2 = sortdens(2,3);                               %y value of center of large size, medium density calcification
        middley3 = sortdens(3,3);                               %y value of center of large size, low density calcification
        circcenter1 = [middlex1;middley1];
        circcenter2 = [middlex2;middley2];
        
        %Find the intersections of two circles with the correct radius with the
        %large HD and MD as centerpoints. The real centerpoint is the intersection
        %closest to the center of the large LD calcification.
        distancecenters = sum((circcenter2-circcenter1).^2);
        coordinate0 = (circcenter1+circcenter2)/2+(ravg^2-ravg^2)/distancecenters/2*(circcenter2-circcenter1);
        overlap = ((ravg+ravg)^2-distancecenters)*(distancecenters-(ravg-ravg)^2);
        if overlap <= 0
            fprintf('The circles don''t intersect.\n')
        else
            coordinatediff = sqrt(overlap)/distancecenters/2*[0 -1;1 0]*(circcenter2-circcenter1);
            intersect1 = coordinate0 + coordinatediff;
            intersect2 = coordinate0 - coordinatediff;
        end
        distanceLD1 = sqrt(abs((middlex3-intersect1(1))^2 + (middley3-intersect1(2))^2));
        distanceLD2 = sqrt(abs((middlex3-intersect2(1))^2 + (middley3-intersect2(2))^2));
        
        if distanceLD1 < distanceLD2
            center_calcifications = [intersect1(1),intersect1(2)];
        else
            center_calcifications = [intersect2(1),intersect2(2)];
        end
        
        for replace = 1:nrofcalc(1)
            deltax = sumdata_Vol(replace,2)-center_calcifications(1,1);
            deltay = sumdata_Vol(replace,3)-center_calcifications(1,2);
            
            %Because the origin (0,0) of the image is in the upper left corner
            %in MATLAB, a lesion above the center has a negative angle and
            %below the center has a positive angle.
            angle = atan2d(deltay, deltax);                     %measure the angle to the center of the calcifications
            centerdistance = sqrt(abs((sumdata_Vol(replace,2)-center_calcifications(1,1))^2 + (sumdata_Vol(replace,3)-center_calcifications(1,2))^2));        %measures the distance to the center of the calcifications
            sumdata_Vol(replace,7) = angle;
            sumdata_Vol(replace,8) = centerdistance;
        end
        
        highdensangle = sumdata_Vol(sumdata_Vol(:,1)==1,7);       %angle to the center of the large size, high density calcification
        highdensdistance = sumdata_Vol(sumdata_Vol(:,1)==1,8);    %distance to the center of the large size, high density calcification
        meddensangle = sumdata_Vol(sumdata_Vol(:,1)==2,7);        %angle to the center of the large size, medium density calcification
        meddensdistance = sumdata_Vol(sumdata_Vol(:,1)==2,8);     %distance to the center of the large size, medium density calcification
        lowdensangle = sumdata_Vol(sumdata_Vol(:,1)==3,7);        %angle to the center of the large size, low density calcification
        lowdensdistance = sumdata_Vol(sumdata_Vol(:,1)==3,8);     %distance to the center of the large size, low density calcification
        
        %This for loop labels the calcification with label 4-9 depending on the
        %position of the calcification with respect to the large calcifications.
        for calcorder = 1:nrofcalc(1)
            if (highdensangle-3) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (highdensangle+3) && (8/pixelsize) <= (highdensdistance - sumdata_Vol(calcorder,8)) && (highdensdistance - sumdata_Vol(calcorder,8)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 4;                    %medium size, high density calcification
            elseif (meddensangle-3) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (meddensangle+3) && (8/pixelsize) <= (meddensdistance - sumdata_Vol(calcorder,8)) && (meddensdistance - sumdata_Vol(calcorder,8)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 5;                    %medium size, medium density calcification
            elseif (lowdensangle-3) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (lowdensangle+3) && (8/pixelsize) <= (lowdensdistance - sumdata_Vol(calcorder,8)) && (lowdensdistance - sumdata_Vol(calcorder,8)) <= (16/pixelsize)
                sumdata_Vol(calcorder,1) = 6;                    %medium size, low density calcification
            elseif (highdensangle-5) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (highdensangle+5) && (21/pixelsize) <= (highdensdistance - sumdata_Vol(calcorder,8)) && (highdensdistance - sumdata_Vol(calcorder,8)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 7;                    %small size, high density calcification
            elseif (meddensangle-5) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (meddensangle+5) && (21/pixelsize) <= (meddensdistance - sumdata_Vol(calcorder,8)) && (meddensdistance - sumdata_Vol(calcorder,8)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 8;                    %small size, medium density calcification
            elseif (lowdensangle-5) <= sumdata_Vol(calcorder,7) && sumdata_Vol(calcorder,7) <= (lowdensangle+5) && (21/pixelsize) <= (lowdensdistance - sumdata_Vol(calcorder,8)) && (lowdensdistance - sumdata_Vol(calcorder,8)) <= (28.5/pixelsize)
                sumdata_Vol(calcorder,1) = 9;                    %small size, low density calcification
            end
        end
        
        %Remove all slices before the first slice with a large calcification and
        %after the last slice with a large calcification. Because there is a chance
        %that some noise will be on exactly the same spot as a medium sized
        %calcification in the slice before the first or after the last, there is a
        %one slice extra range on both sides.
        large = find(ismember(labels, sortarea(1:3,10)));
        beginslice = lesion_info_Vol(large(1),1);
        endslice = lesion_info_Vol(large(end),1);
        slices = lesion_info_Vol(:,1);
        remove1 = find(slices<beginslice-1);
        remove2 = find(slices>endslice+1);
        if isempty(remove1) && isempty(remove2)
            
        elseif isempty(remove1)
            L2 = find(sumdata_Vol(:,10)==lesion_info_Vol(remove2(1),11),1);
            sumdata_Vol(L2:end,:) = [];
        elseif isempty(remove2)
            L1 = find(sumdata_Vol(:,10)==lesion_info_Vol(remove1(end),11),1);
            sumdata_Vol(1:L1,:) = [];
        else
            L1 = find(sumdata_Vol(:,10)==lesion_info_Vol(remove1(end),11),1);
            L2 = find(sumdata_Vol(:,10)==lesion_info_Vol(remove2(1),11),1);
            sumdata_Vol(L2:end,:) = [];
            sumdata_Vol(1:L1,:) = [];
        end
        noise = sumdata_Vol(:,1) == 0;                            %find all noise in the slices
        sumdata_Vol(noise,:) = [];                                %remove all noise from sumdata
        sumdata_Vol = sortrows(sumdata_Vol,1);                    %sort sumdata in order of calcification number
        size_Ag = size(sumdata_Ag);
        size_Vol = size(sumdata_Vol);
        matrix_diff = size_Vol(1) - size_Ag(1);
        if matrix_diff > 0
            sumdata_Ag(size_Ag(1)+matrix_diff,:) = zeros;
        elseif matrix_diff < 0
            sumdata_Ag(size_Ag(1)+matrix_diff:end,6) = 0;
        end
        sumdata_Ag(:,6) = sumdata_Vol(:,5);
        sumdata_Ag(:,7) = sumdata_Vol(:,6);
    end
    
    %Proceed with CAC scores:
    c = 0.2/(CTc-CTw);                                      %calibrationfactor for the mass score with 200 mg/cm3 = 0.2 mg/mm3 as known density of the calibration rod
    lesion_info_Ag(:,10:11) = lesion_info_Ag(:,10:11)*c;    %multiply the calibrationfactor with the mass scores calculated before to get the final mass score
    sumdata_Ag(:,7) = sumdata_Ag(:,7)*c;
    
    Agatston = round(sum(sumdata_Ag(:,5)));                 %sum to get total Agatston score
    Volume = round(sum(sumdata_Ag(:,6)));                   %sum to get total Volume score
    Mass = round(sum(sumdata_Ag(:,7)));                     %sum to get total Mass score
    nroflesions = size(sumdata_Ag);                         %number of calcifications visible
    nrofnoise = nrofnoise(1)-nroflesions(1)+sum(small==1);  %number of false positives
    
    if contains(mf,'TOSHIBA') || contains(mf,'Canon')
        scantype = [kVp, mAs, slice_thickness, level, time, c, sort, 0, 0, 0, 0, 0];
    else
        scantype = [kVp, mAs, slice_thickness, level, time, c, 0, 0, 0, 0, 0, 0];
    end
    sumdata_Ag = [scantype; sumdata_Ag];
    
    if contains(mf,'GE MEDICAL SYSTEMS') || contains(mf,'Imatron')
        Data_Ag = [{fn}, {model}, {description}, kVp, mAs, slice_thickness, level, time, Agatston, Volume, Mass, nroflesions(1), nrofnoise, SD_noise];
    elseif contains(mf,'TOSHIBA') || contains(mf,'Canon')
        Data_Ag = [{fn}, {model}, sort, kVp, mAs, slice_thickness, level, CTDIvol, Agatston, Volume, Mass, nroflesions(1), nrofnoise, SD_noise];
    else
        Data_Ag = [{fn}, {model}, {description}, kVp, mAs, slice_thickness, level, CTDIvol, Agatston, Volume, Mass, nroflesions(1), nrofnoise, SD_noise];
    end
    
    
    
    %TTF:
    [TTF] = ttf(centerslice_calibration, center_calibration, pixelsize);
    
    %NPS
    [NPS] = nps(homogeneous_slice,center_hs,pixelsize,32,15,15,15);
    
    toc
end