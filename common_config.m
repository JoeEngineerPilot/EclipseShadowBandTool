function config = common_config

% Configuration parameters
% These need to be updated for each video recording
config = [];
config.block_size = 3;
config.videofile = 'eclipse-shadow-bands-video.avi';
config.xform_rotate = -5;
config.xform_corners = [935 330; 2560 345; 3100 1170; 620 1200;];
config.xform_crop = [2900 1940 3890 3160];
config.gaussfilt = 15;
config.movmean = 10;
config.fontsize=14;
config.fontname='Times New Roman';
config.splinesize = 40;
