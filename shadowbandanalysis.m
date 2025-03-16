% Import raw video file of eclipse shadow bands taken with 
% a phone camera on a stand
clear all; 
close all force;
block_size = 3; % number of frames to load from video
start_frame = 2;
video1 = VideoReader('eclipse-shadow-bands-video.avi');
for n = 1:block_size
    videoBlock(:,:,:,n) = read(video1,start_frame+n); %#ok<*SAGROW>
end

% Display last frame of video block in figure
figure
imshow(videoBlock(:,:,:,end)); 
title('Raw Video')

% STEP 1: Perspective Transformation
% The first step is to transform the video from a 3-D perspective 
% to a planar 2-D view. This transformation ensures subsequent 
% shadow band geometry data is collected in a spatially consistent manner. 
xform_rotate = -5;
xform_corners = [935 330; 2560 345; 3100 1170; 620 1200;];
xform_crop = [2900 1940 3890 3160];
xFormBlock = uint8([]);
for n = 1:block_size
    sz = size(videoBlock(:,:,:,n));
    fixedPoints = [0 0; sz(2) 0; sz(2) sz(1); 0 sz(1)];
    xform_perspective = fitgeotform2d(xform_corners,fixedPoints,'projective');
    temp_image = imrotate(videoBlock(:,:,:,n),xform_rotate);
    temp_uncropped_image = imwarp(temp_image,xform_perspective);
    xFormBlock(:,:,:,n) = imcrop(temp_uncropped_image,xform_crop);
end

% Display in figure
figure
imshow(xFormBlock(:,:,:,end));
title('2-D Transformation')

% STEP 2: Grayscale Transformation
% The next step is to convert the video from RGB to grayscale 
% of intensity values. Overall intensity is needed in order to 
% simplify extraction of shadow band features. 
videoBlockGray = uint8([]);
for n = 1:block_size
    videoBlockGray(:,:,n) = rgb2gray(xFormBlock(:,:,:,n));
end

% STEP 3: Masking
% Mask pixels outside of the white sheet with NaNs
for n = 1:block_size
    temp_image = double(videoBlockGray(:,:,n));
    mask_ind = find(temp_image<=200);
    temp_image(mask_ind) = nan;
    videoBlockGrayMasked(:,:,n) = temp_image;
end

% STEP 4: Temporal Flat Field Correction
% Transform all image frames within video block to have the same
% median intensity. This accounts for ambient light changing over 
% time as the eclipse approaches. Median is used, instead of mean, 
% since the image intensity values were originally captured 
% as uint8 (255 intensity values).
medianBaseline = median(videoBlockGrayMasked(:,:,end),'all','omitnan');
medianAdjustedVideoBlockGray = [];
for n = 1:block_size
    medianOffset = median(videoBlockGrayMasked(:,:,n),'all','omitnan') - medianBaseline;
    medianAdjustedVideoBlockGray(:,:,n) = videoBlockGrayMasked(:,:,n) - medianOffset;
end
flatFieldBaselineImg = median(double(medianAdjustedVideoBlockGray(:,:,:)),3);

% STEP 5: Temporal Differential
% Subtract flat field corrected frames from the video with the 
% flat field baseline. 
for n = 1:block_size
    timediff(:,:,n) = medianAdjustedVideoBlockGray(:,:,n)-flatFieldBaselineImg;
end

% STEP 6: Smoothing
% Use 2-D convolution to smooth image. Nans are replaced with 
% baseline to avoid aliasing effects from smoothing
timediff_smoothed = [];
timediff_smoothed_nan = [];
for n = 1:block_size
    temp_image = timediff(:,:,n);
    replacement_value = median(timediff(:,:,n),'all','omitnan');
    nanind = find(isnan(temp_image));
    temp_image(nanind) = replacement_value;
    s = 100;
    temp_image_smoothed = conv2(temp_image,ones(s,s)./(s*s),'same');
    timediff_smoothed(:,:,n) = temp_image_smoothed;
    temp_image_smoothed(nanind) = nan;
    timediff_smoothed_nan(:,:,n) = temp_image_smoothed;
end

% Display in figure
figure
imshow(timediff_smoothed(:,:,end));
title('Temporal Differential');

% STEP 7: Spectral Analysis
% Calculate the scintillation index for the last frame of
% the loaded video
avg_intensity = mean(timediff_smoothed(:,:,end),'all');
avg_intensity_squared = mean(timediff_smoothed(:,:,end).*timediff_smoothed(:,:,n),'all');
scintillation_index = (avg_intensity_squared - avg_intensity^2);

% Transform video block to binary 
timediff_binary = [];
for n = 1:block_size
    med_img = median(timediff_smoothed(:,:,n),'all');
    std_img = std(timediff_smoothed(:,:,n),[],'all');
    timediff_binary(:,:,n) = imbinarize(timediff_smoothed(:,:,n),med_img-std_img);
end

% STEP 8: Geometry Analysis
% The next step is to extract geometric properties of shadow 
% bands on the image. By modeling each shadow band as an eclipse 
% we can compute the major and minor axis as well the orientation.
figure
imshow(timediff_binary(:,:,end));
title('Geometry Analysis')
BWP = bwperim(timediff_binary(:,:,end),8);
s = regionprops(BWP,'centroid');
stats = regionprops("table",BWP,"Solidity","Centroid","Eccentricity","MajorAxisLength","MinorAxisLength","Area","BoundingBox","Orientation");
band_angle_sum = [];

% Sort shadow bands from biggest to smallest
sortedstats = sortrows(stats,'Area','descend');
for k = 2:height(sortedstats)
    
    x = sortedstats.Centroid(k,1);
    y = sortedstats.Centroid(k,2);
    d = sortedstats.Orientation(k);
    m1 = sortedstats.MajorAxisLength(k);
    m2 = sortedstats.MinorAxisLength(k);

    band_angle_sum(end+1) = d;

    x1 = x-m2/2*cos(deg2rad(d+90));
    y1 = y+m2/2*sin(deg2rad(d+90));
    x2 = x+m2/2*cos(deg2rad(d+90));
    y2 = y-m2/2*sin(deg2rad(d+90));

    x3 = x-m1/2*cos(deg2rad(d));
    y3 = y+m1/2*sin(deg2rad(d));
    x4 = x+m1/2*cos(deg2rad(d));
    y4 = y-m1/2*sin(deg2rad(d));

    % Draw a red dot at the center
    line('XData',x,'Ydata',y,'Color','red','Marker','o','MarkerSize',10);
    % Draw a red line on the minor axis of each shadow band
    line('XData',[x1,x2],'Ydata',[y1,y2],'Color','red','LineWidth',2);
    % Draw a blue line on the major axis of each shadow band
    line('XData',[x3,x4],'Ydata',[y3,y4],'Color','blue','LineWidth',2);
end

band_angle_median = median(band_angle_sum);
band_angle_std = std(band_angle_sum);

% STEP 9: Velocity Analysis
% Determine velocity of shadow bands through image registration 
% across image frames. Image registration can tell us how
% much the image is translating in pixels. Pixel translation
% values are then converted to meters per second since we
% know the size of the white sheet in meters and the frame
% rate of the video in frames per second.
currImg = timediff_binary(:,:,end);
prevImg = timediff_binary(:,:,end-1);
fixedRefObj = imref2d(size(prevImg));
movingRefObj = imref2d(size(currImg));
[tform, peakCorr]= imregcorr(currImg,movingRefObj,prevImg,fixedRefObj,'Method','gradcorr','transformtype','rigid');
regImg = imwarp(currImg, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
vx = tform.Translation(1); 
vy = tform.Translation(2); 
shadowBandVelocity = sqrt(vx^2+vy^2); % pixels per frame
sheet_height = (34 * 5.5/2)/100; % meters
pixels_per_meter = 3062/sheet_height;
frames_per_second = 60;
shadowBandVelocityMetersPerSec = (shadowBandVelocity/pixels_per_meter) * frames_per_second;
figure
imshowpair(prevImg,regImg,"falsecolor");
title('Velocity Analysis');

% Transform extracted shadow bands back to 3-D perspective and overlay 
% onto raw video
inv_xform_perspective = invert(xform_perspective);
temp_image = imwarp(double(~timediff_binary(:,:,end)),inv_xform_perspective);
temp_image = imrotate(temp_image,-xform_rotate);
temp_raw_frame = double(videoBlock(:,:,:,end));
temp_overlay = cat(3, temp_image(:,:,end),temp_image(:,:,end),temp_image(:,:,end));
yshift = 40;
xshift = 150;
sz = size(temp_overlay);
temp_combined = temp_raw_frame(yshift+1:yshift+sz(1),xshift+1:xshift+sz(2),:) - temp_overlay*150;

% Display in figure
figure
imshow(uint8(temp_combined));
title('Shadow Band Overlay');