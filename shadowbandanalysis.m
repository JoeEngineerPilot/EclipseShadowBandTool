% Import raw video file of eclipse shadow bands taken with
% a phone camera on a stand
clear all;
close all force;
block_size = 3; % number of frames to load from video
video1 = VideoReader('eclipse-shadow-bands-video.avi');
for n = 1:block_size
    videoBlock(:,:,:,n) = read(video1,n); %#ok<*SAGROW>
end

% Display last frame of video block in figure
figure
imshow(videoBlock(:,:,:,end));
title('Raw Video Frame')

% Transform and Crop
% Transform the video from a 3-D perspective
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
title('Perspective and Crop Transformation')

% Grayscale Transformation
% Convert the video from RGB to grayscale
% of intensity values. Overall intensity is needed in order to
% simplify extraction of shadow band features.
videoBlockGray = uint8([]); %#ok<UNRCH>
for n = 1:block_size
    videoBlockGray(:,:,n) = rgb2gray(xFormBlock(:,:,:,n));
end

% Masking: Mask pixels outside of the white sheet with NaNs
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


% Temporal Normalization
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

% Temporal Differential
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
    temp_image_smoothed = imgaussfilt(temp_image,15);
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
si = I_var / (I_mean^2)

% Calculate Orientation Distribution Function (ODF)
orientations = [];
im = timediff_smoothed_normalized;

% Edge detection using gradient technique
imBW = edge(im,'Canny',.1,30);
retval.edgeBW = imBW;

% Remove laptop bag effects from image
imBWMasked = imBW;
imBWMasked(1:900,1:900) = 0;

% Display in figure
figure
imshow(imBWMasked);
title('Edge Segmentation');

% Get boundaries
B = bwboundaries(imBWMasked,8,'noholes');
retval.B = B;

% Get stats on boundaries
statsB = [];
for k = 1:length(B)
    statsB(k) = length(B{k});
end

% Get quantile of boundary length, this will be used to
% calculate histogram
contour_spline_size = 40;

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
histysmooth = movmean(counts',10);
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
    retval.a1 = nan;
    retval.a1ci = nan;
    retval.b1 = nan;
    retval.b1ci = nan;
    retval.c1 = nan;
    retval.c1ci = nan;
    retval.snr1 = nan;
    retval.snr1ci = nan;
else
    retval.a1 = ft.a1;
    retval.a1ci = cof(2,1)-cof(1,1);
    retval.b1 = ft.b1;
    retval.b1ci= cof(2,3)-cof(1,3);
    retval.c1 = ft.c1;
    retval.c1ci = cof(2,5)-cof(1,5);
    retval.snr1 = ft.a1/ft.c1;
    % Calculate approximate confidence interval (error) for snr1 using
    % standard propagation of uncertainty assuming independent errors:
    % snr = a1 / c1
    % d(snr) = snr * sqrt( (da1/a1)^2 + (dc1/c1)^2 )
    da1 = retval.a1ci;
    dc1 = retval.c1ci;
    if retval.a1 ~= 0 && retval.c1 ~= 0 && ~isnan(da1) && ~isnan(dc1)
        retval.snr1ci = retval.snr1 * sqrt( (da1/retval.a1)^2 + (dc1/retval.c1)^2 );
    else
        retval.snr1ci = nan;
    end
end

% Ignore boundary conditions
if ft.b2>180 || ft.b2<0
    retval.a2 = nan;
    retval.a2ci= nan;
    retval.b2 = nan;
    retval.b2ci = nan;
    retval.c2 = nan;
    retval.c2ci = nan;
    retval.snr2 = nan;
    retval.snr2ci = nan;
else
    retval.a2 = ft.a2;
    retval.a2ci = cof(2,2)-cof(1,2);
    retval.b2 = ft.b2;
    retval.b2ci = cof(2,4)-cof(1,4);
    retval.c2 = ft.c2;
    retval.c2ci= cof(2,6)-cof(1,6);
    retval.snr2 = ft.a2/ft.c2;
    % Calculate approximate confidence interval (error) for snr1 using
    % standard propagation of uncertainty assuming independent errors:
    % snr = a1 / c1
    % d(snr) = snr * sqrt( (da1/a1)^2 + (dc1/c1)^2 )
    da2 = retval.a2ci;
    dc2 = retval.c2ci;
    if retval.a2 ~= 0 && retval.c2 ~= 0 && ~isnan(da2) && ~isnan(dc2)
        retval.snr2ci = retval.snr2 * sqrt( (da2/retval.a2)^2 + (dc2/retval.c2)^2 );
    else
        retval.snr2ci = nan;
    end
end

indstart = find(xprime==0);
indstop = find(xprime==180);
retval.x = xprime(indstart:indstop);
retval.y = ywrap(indstart:indstop);
% Evaluate fitted curve on the displayed x-range
retval.yft = ft(xprime(indstart:indstop));
% Calculate 95% confidence bounds for the fit over that range
[yfit] = predint(ft, xprime(indstart:indstop), 0.95, 'functional', 'on');
retval.yfit_ci_lower = yfit(:,1);
retval.yfit_ci_upper = yfit(:,2);

% Display in figure
figure
plot(retval.x, retval.y, 'k-', 'DisplayName','Data');
hold on;
hfit = plot(retval.x, retval.yft, 'r-', 'LineWidth',1.5, 'DisplayName','Fit');
% Ensure the fit line is on top
uistack(hfit,'top')
legend('show','Location','best')
xlabel('Orientation (degrees)')
ylabel('Response')
hold off;
title('Orientation Distribution Function')

if isnan(retval.b1)
    b1 = 0;
else
    b1 = 90-retval.b1;
end

% Rotate image so that shadow bands are vertical and convert to double
Irot = imrotate(im, theta, 'bilinear');
retval.irot = Irot;

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

% Compenstate for window power
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
retval.f = f';
retval.psd = psd;
retval.phase = phase_rad;

% Display in figure
figure
loglog(f, real(10*log10(psd)));
axis tight;
title('Power Spectral Density along Primary Orientation');
