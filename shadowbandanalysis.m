% This script demonstrates a framework for analyzing eclipse shadow 
% band videos. For more information, see the following reference:
% Conti, Joseph. "A quantitative image-processing framework for eclipse 
% shadow band analysis." Experimental Astronomy 62.3 (2026): 25.

% Clear out any variables and figures
clear all;
close all force;

% Load configuration used for this demo
config = common_config;

% Load block of images from video file
block_size = config.block_size; 
video1 = VideoReader(config.videofile);
for n = 1:block_size
    videoBlock(:,:,:,n) = read(video1,n); 
end

% Display last frame of video block in figure
figure
imshow(videoBlock(:,:,:,end));
title('Raw Video Frame')

% Transform the video from a 3-D perspective
% to an orthographic 2-D view. This transformation ensures subsequent
% shadow band geometry data is collected in a spatially consistent manner.
xform_rotate = config.xform_rotate;
xform_corners = config.xform_corners;
xform_crop = config.xform_crop;
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
imshow(temp_uncropped_image);
title('Orthographic Transformation of Image')

% Display in figure
figure
imshow(xFormBlock(:,:,:,end));
title('Cropped Image')

% Grayscale transformation
% Convert the video from RGB to grayscale
% of intensity values. Overall intensity is needed in order to
% simplify extraction of shadow band features.
videoBlockGray = uint8([]); 
for n = 1:block_size
    videoBlockGray(:,:,n) = rgb2gray(xFormBlock(:,:,:,n));
end

% Mask pixels outside of the white sheet with NaNs
for n = 1:block_size
    temp_image = double(videoBlockGray(:,:,n));
    mask_ind = find(temp_image<=200);
    temp_image(mask_ind) = nan;
    videoBlockGrayMasked(:,:,n) = temp_image;
end

% Display in figure
figure
imshow(videoBlockGrayMasked(:,:,:,end));
title('Grayscale Transformation')

% Temporal normalization
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

% Image Differencing
% Subtract flat field corrected frames from the video with the
% flat field baseline.
for n = 1:block_size
    timediff(:,:,n) = medianAdjustedVideoBlockGray(:,:,n)-flatFieldBaselineImg;
end

% Display in figure
figure
temp_im = timediff(:,:,end);
mn = mean(temp_im,'all','omitnan');
temp_im_normalized = temp_im-mn+.5;
imshow(temp_im_normalized);
title('Image Differencing')

% Smoothing
% Use Gaussian convolution to smooth image. Nans are replaced with
% baseline to avoid aliasing effects from smoothing
timediff_smoothed = [];
for n = 1:block_size
    temp_image = timediff(:,:,n);
    temp_image_smoothed = imgaussfilt(temp_image,config.gaussfilt);
    timediff_smoothed(:,:,n) = temp_image_smoothed;
end

% Normalize data
mn = mean(temp_image_smoothed,'all','omitnan');
timediff_smoothed_normalized = temp_image_smoothed-mn+.5;

% Display in figure
figure
imshow(timediff_smoothed_normalized);
title('Image Smoothing');

% Scintillation index
I = timediff_smoothed_normalized;
I_mean = mean(I(:),'omitnan');
I_var = var(I(:),1,'omitnan');
si = I_var / (I_mean^2);

% Calculate Orientation Distribution Function (ODF)
orientations = [];
im = timediff_smoothed_normalized;

% Edge detection using gradient technique
imBW = edge(im,'Canny',.1,30);
edgeBW = imBW;

% Display in figure
figure
imshow(imBW);
title('Edge Segmentation');

% Get boundaries
B = bwboundaries(imBW,8,'noholes');
B = B;

% Get stats on boundaries
statsB = [];
for k = 1:length(B)
    statsB(k) = length(B{k});
end

% Get quantile of boundary length, this will be used to
% calculate histogram
contour_spline_size = config.splinesize;

% Get boundary data for histogram
for k = 1:length(B)
    b = B{k};
    indb = 1;
    while indb+contour_spline_size <= size(b,1)
        if size(b,1) >= indb+contour_spline_size
            dx = b(indb+contour_spline_size,2)-b(indb,2);
            dy = b(indb+contour_spline_size,1)-b(indb,1);
            indb = indb + contour_spline_size;
        end
        theta = atan2d(-dy,dx);
        orientations(end+1) = theta;
    end
end

% Calculate histogram
orientations(orientations < 0) = orientations(orientations < 0) + 180;
edges =  0:2:180;
counts = histcounts(orientations,'BinEdges',edges);

% Get x,y data from histogram
edgexy = B;
histy = counts';
histysmooth = movmean(counts',config.movmean);
histx = edges(1:length(counts))';
x = histx;
y = histy;
ysmooth = histysmooth;
origlen = size(x,1);
wraplen = 40;
s = 2;
xprime = (-wraplen*s:s:(90+wraplen)*s)';
ywrap = [y(origlen-wraplen:1:origlen); y(1:origlen); y(1:wraplen)];
ysmoothwrap = [ysmooth(origlen-wraplen:1:origlen); ysmooth(1:origlen); ysmooth(1:wraplen)];

% Fit data to gaussian and constant
ftype = fittype('a1*exp(-((x-b1)/c1)^2) + a1*exp(-((x-b1-180)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a2*exp(-((x-b2+180)/c2)^2) + d','independent', 'x');
options = fitoptions(ftype);

% Determine approximate location of first peak
[~,loc,~,~] = findpeaks(ysmoothwrap(wraplen:end-wraplen),xprime(wraplen:end-wraplen),'MinPeakDistance',45,'NPeaks',2,'SortStr','descend');

% Determine approximate location of first peak
if ~isempty(loc)>0
    bcenter1 = loc(1);
    [~, indmax] = max(ywrap);
    if indmax>wraplen
        bmax = xprime(indmax-wraplen);
    else
        bmax = xprime(indmax+wraplen);
    end
    if abs(bcenter1-bmax) < 20
        bcenter1 = bmax;
    end
    blower1 = bcenter1-3;
    bupper1 = bcenter1+3;
else
    blower1 = 0;
    bupper1 = 180;
    bcenter1 = 90;
end

% Determine approximate location of second peak
if ~isempty(loc) && length(loc)>1
    bcenter2 = loc(2);
    blower2 = bcenter2-20;
    bupper2 = bcenter2+20;
else
    if bcenter1 <= 90
        blower2 = bcenter1+45;
        bupper2 = 180;
    else
        blower2 = 0;
        bupper2 = bcenter1-45;
    end
    bcenter2 = blower2+(bupper2-blower2)/2;
end
options.Lower = [0, 0, blower1, blower2, 10, 10, -inf];
options.Upper = [max(ysmoothwrap)-min(ysmoothwrap), ...
    max(ysmoothwrap)-min(ysmoothwrap), ...
    bupper1, bupper2, inf, inf, inf];
options.StartPoint = [50, 50, bcenter1, bcenter2, 10, 10, 0];

[ft] = fit(xprime, ysmoothwrap, ftype, options);
cof = confint(ft);

% Ignore boundary conditions
if ft.b1>180 || ft.b1<0
    a1 = nan;
    a1ci = nan;
    b1 = nan;
    b1ci = nan;
    c1 = nan;
    c1ci = nan;
    snr1 = nan;
    snr1ci = nan;
else
    a1 = ft.a1;
    a1ci = cof(2,1)-cof(1,1);
    b1 = ft.b1;
    b1ci= cof(2,3)-cof(1,3);
    c1 = ft.c1;
    c1ci = cof(2,5)-cof(1,5);
    snr1 = ft.a1/ft.c1;
    da1 = a1ci;
    dc1 = c1ci;
    if a1 ~= 0 && c1 ~= 0 && ~isnan(da1) && ~isnan(dc1)
        snr1ci = snr1 * sqrt( (da1/a1)^2 + (dc1/c1)^2 );
    else
        snr1ci = nan;
    end
end

% Ignore boundary conditions
if ft.b2>180 || ft.b2<0
    a2 = nan;
    a2ci= nan;
    b2 = nan;
    b2ci = nan;
    c2 = nan;
    c2ci = nan;
    snr2 = nan;
    snr2ci = nan;
else
    a2 = ft.a2;
    a2ci = cof(2,2)-cof(1,2);
    b2 = ft.b2;
    b2ci = cof(2,4)-cof(1,4);
    c2 = ft.c2;
    c2ci= cof(2,6)-cof(1,6);
    snr2 = ft.a2/ft.c2;
    da2 = a2ci;
    dc2 = c2ci;
    if a2 ~= 0 && c2 ~= 0 && ~isnan(da2) && ~isnan(dc2)
        snr2ci = snr2 * sqrt( (da2/a2)^2 + (dc2/c2)^2 );
    else
        snr2ci = nan;
    end
end

indstart = find(xprime==0);
indstop = find(xprime==180);
x = xprime(indstart:indstop);
y = ywrap(indstart:indstop);

% Evaluate fitted curve on the displayed x-range
yft = ft(xprime(indstart:indstop));

% Calculate 95% confidence bounds for the fit over that range
[yfit] = predint(ft, xprime(indstart:indstop), 0.95, 'functional', 'on');
yfit_ci_lower = yfit(:,1);
yfit_ci_upper = yfit(:,2);

% Display in figure
figure
plot(x, y, 'k-', 'DisplayName','Data');
hold on;
hfit = plot(x, yft, 'r-', 'LineWidth',1.5, 'DisplayName','Fit');
% Ensure the fit line is on top
uistack(hfit,'top')
legend('show','Location','best')
xlabel('Orientation (degrees)')
ylabel('Response')
hold off;
title('Orientation Distribution Function')

if isnan(b1)
    b1 = 0;
else
    b1 = 90-b1;
end

% Rotate image so that shadow bands are vertical and convert to double
Irot = imrotate(im, theta, 'bilinear');
irot = Irot;

% Display in figure
figure
imshow(Irot)
title('Image Rotated Along Band Orientation')

% Collapse central rows to a single 1-D profile (mean across rows)
x = mean(Irot, 1, 'omitnan'); % units = pixel intensity

% Apply Hann window across profile to reduce spectral leakage
N = numel(x);
w = hann(N)'; % units = dimensionless
xw = x .* w; % units = pixel intensity

% FFT
X = fft(xw);

% one-sided power spectral density (averaged)
P = (abs(X)/N).^2; % units = power per bin

% Correct for window power
U = sum(w.^2) / N; % units = dimensionless

% Single-sided PSD (real signal)
if rem(N,2)==0
    % N even: include Nyquist
    P1 = P(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
else
    P1 = P(1:(N+1)/2);
    P1(2:end) = 2*P1(2:end);
end

% Compensate for window power
psd = P1 / U; % units = power per bin

% Plot PSD vs image spatial frequency (cycles per meter)
fs = 1; % sampling: 1 sample per pixel -> cycles/pixel (fs/N used below)
f_pixels = (0:numel(P1)-1)' * (fs / N); % frequency vector (cycles/pixel)
% convert cycles/pixel to cycles/meter using image scaling
imageWidthPerMeter = 1;
pixels_per_meter = (N * imageWidthPerMeter) .* (N / size(im,2));
f = f_pixels * pixels_per_meter; % cycles per meter

% Single-sided complex spectrum (keep complex values for phase)
if rem(N,2)==0
    X1 = X(1:N/2+1);
else
    X1 = X(1:(N+1)/2);
end

% Phase (radians) and unwrapped phase
phase_rad = angle(X1);  % principal value in [-pi, pi]

% Return results
f = f';
psd = psd;
phase = phase_rad;

% Display in figure
figure
loglog(f, real(10*log10(psd)));
axis tight;
title('Power Spectral Density along Primary Band Orientation');
