% clear all; close all; clc;

% ***** 01_VOR_0 *************************
% AUTHOR:       Zachary-Jacques Gray
% DATE CREATED: 14 Jan 2022
% DESCRIPTION:  MATLAB script to detect streams and droplets from image via
%               clustering and edge-detection algorithms.
%               Analyzes for injector flow qualities, e.g. atomisation and
%               oscillations.
% REQUIREMENTS: MATLAB Image Processing Toolbox
%               MATLAB Signal Processing Toolbox

tic;
% BACKGROUND ---------------------------//
imgB = imread('01_VOR_9B_20211217_1535230001.tif');
greyB = im2gray(imgB);
L = 2^16; % 16-bit Max Intensity
greyB = (L - 1) - greyB; % negative
% figure()
% imshow(greyB);
% title('VOR 0 BACKGROUND');


% LOOP FOR EACH IMAGE ===========================================//
n_images = 2400; % number of frames to run starting from image 1

WelchSolve = zeros(129,n_images); % may be different for each new data set
WelchSolveSmooth = WelchSolve; % may be different for each new data set
coneAngles = zeros(1,n_images);
amp = zeros(1,n_images);
freq = zeros(1,n_images);

for frame = 1:n_images
    imgstr = '01_VOR_9T_20211217_153652';
    imgnum = num2str(frame);
    imgtyp = '.tif';
    if frame < 10
        zeroes = '000';
    elseif frame < 100
        zeroes = '00';
    elseif frame < 1000
        zeroes = '0';
    elseif frame < 10000
        zeroes = '';
    else
        warning('Error in filename zeroes')
    end
img = strcat(imgstr,zeroes,imgnum,imgtyp);


% IMAGE N ------------------------------//
% img1 = imread('01_VOR_0T_20211217_1438460166.tif');
img1 = imread(img);
grey1 = im2gray(img1);
grey1 = (L - 1) - grey1; % negative
% figure()
% imshow(grey1);
% title('VOR 0 1 RAW');

% IMAGE 1 BACKGROUND SUBTRACTION -------//
Z = imsubtract(grey1,greyB);
% figure()
% imshow(Z);

% remove common bad sections from image
for i = 1:1024
    for j = 1:1024
        % top left and bottom right corners
        if (i + j) > 1024*1.55 || (i+j) < 1024*0.45
            Z(i,j) = 0;
        end
        % far left middle
        if (j<120)
            Z(i,j) = 0;
        end
        % small bubble hanging from injector
        if (j>850)
            Z(i,j) = 0;
        end
    end
end
% adjust image for clustering algorithm
Z = imadjust(Z,[0.16,0.7],[0,1],0.7); % gamma and contrast correction for viewing
% figure()
% imshow(Z);

for i = 1:1024
    for j = 1:1024
        if Z(i,j) > 0
      % if Z(i,j) > 10000 % for use if contrast and gamma correction is not used
            Z(i,j) = L;
        end
    end
end
% figure()
% imshow(Z);
% title('VOR 0 1');


% CLUSTERING ALGORITHM ------------------//
% MATLAB bwlabel > > >

[label,N] = bwlabel(Z);
label_count = zeros(1,N);
[m,n] = size(Z);
large_cluster = [0,0]; % to store index of the main stream clusters
iter = 1;

% Separate out largest cluster
for c = 1:N
    label_count(c) = sum(label(:) == c);
    if label_count(c) > 0.025*m*n % assumes main streams are at least 2.5% of image
        large_cluster(iter) = c;
        iter = iter + 1;
    end
end
if large_cluster(2) == 0
    large_cluster = large_cluster(1);
end
% [a,I] = max(label_count); % size and index of largest cluster


Z_cluster = Z; % Stream & Droplets
Cluster = Z; % Main Stream(s)
for i = 1:m
    for j = 1:n
        % Visually separate out largest cluster
        Cluster(i,j) = 0; % clear image data
        for c = 1:length(large_cluster)
            if label(i,j) == large_cluster(c)
                label(i,j) = N;
                Cluster(i,j) = L; % Max intensity for main stream
            end
        end
        % Adjust to better visualise clustering on black background
        if label(i,j) ~= 0
            Z_cluster(i,j) = round(label(i,j)*(L./N))/2 + L/2;
        end
    end
end
% figure()
% imshow(Z_cluster); % Streams and Droplets
% figure()
% imshow(Cluster); % Streams


% EDGE DETECTION ------------------------//
% Canny Edge > > >

% First rotate image
rotateAngle = 46.4; % degrees
Cluster45 = imrotate(Cluster, rotateAngle,'bilinear','crop');

% Edge Detection
CannyEdge = edge(Cluster45,'canny');
% figure()
% imshow(CannyEdge);


% FOURIER TRANSFORM ---------------------//
% Welch's PSD Estimate > > >

% First convert segment of stream into (x,y)-function
m45 = [250,900]; % first and last x-range
x = m45(1):m45(2); % initialise x-range
x = x - m45(1);
y = zeros(size(x)); % initialise y-range
for i = m45(1):m45(2)
    for j = 300:530 % top range
        if CannyEdge(j,i) == 1 % if line exists at that pixel
            if y(i-m45(1)+1) == 0
                y(i-m45(1)+1) = j; %bottom most part of line (most likely stream)
            else
                y(i-m45(1)+1) = (y(i-m45(1)+1) + j)/2;
            end
        end
    end
end
y = -y + mean(y); % flip to match same way as edge

ysmooth = smooth(x,y,0.05,'rloess');

% figure()
% plot(x,y,x,ysmooth)

% Perform Sine FFT Fitting
SineParameters = sineFit(x,y);
amp(frame) = SineParameters(2);
freq(frame) = SineParameters(3);

% Perform Welch analysis
% figure()
% WelchSolve(:,frame) = pwelch(y);
% WelchSolveSmooth(:,frame) = pwelch(y,200,50);
% pwelch(y) % single frame plot with automatic details/labels



% DETERMINE CONE ANGLE ------------------//

% First fit a linear slope to the data
[p,S] = polyfit(x,y,1);
[y_fit,delta] = polyval(p,x,S);
slopeAngle = rad2deg(atan(p(1)));
coneAngles(frame) = rotateAngle - slopeAngle;

% figure
% plot(x,ysmooth,'k-')
% hold on
% plot(x,y_fit,'r-')
% plot(x,y_fit+delta,'m--',x,y_fit-delta,'m--')
% title('Linear Fit of Stream Edge')
% legend('x','y','Standard Error')



end 
% *******************// image loop end //***********************

% COMBINED DATA ANALYSIS ----------------//
fprintf('\n> > > VOR 9 > > >')

% FREQUENCY AND AMPLITUDE
amplitude = rmoutliers(amp);
frequency = rmoutliers(freq);
MeanAmp = mean(amplitude);
MeanFreq = mean(frequency);
stdAmp = std(amplitude);
stdFreq = std(frequency);
fprintf('\nAmplitude: %g pts \tFrequency: %g Hz',MeanAmp,MeanFreq)
fprintf('\nSTDEV: %g pts \tSTDEV: %g Hz\n----//----\n',stdAmp,stdFreq)

figure % amplitude scatter
plot(1:length(amplitude),amplitude,'kx')
hold on
plot([1,length(amplitude)],MeanAmp.*[1,1],'r-')
plot([1,length(amplitude)],(MeanAmp+stdAmp).*[1,1],'m--')
plot([1,length(amplitude)],(MeanAmp-stdAmp).*[1,1],'m--')
title('VOR9 Stream Oscillation AMPLITUDE Scatter Plot')
xlabel('Frame Number')
ylabel('Amplitude / pts')
hold off

figure % frequency scatter
plot(1:length(frequency),frequency,'kx')
hold on
plot([1,length(frequency)],MeanFreq.*[1,1],'r-')
plot([1,length(frequency)],(MeanFreq+stdFreq).*[1,1],'m--')
plot([1,length(frequency)],(MeanFreq-stdFreq).*[1,1],'m--')
title('VOR9 Stream Oscillation FREQUENCY Scatter Plot')
xlabel('Frame Number')
ylabel('Frequency / Hz')
hold off

% CONE ANGLE
coneAngle = rmoutliers(coneAngles);
MeanCone = mean(coneAngle);
stdCone = std(coneAngle);
fprintf('Cone Angle: %.3g deg\tSTDEV: %g deg\n----//----\n',MeanCone,stdCone)

figure % cone angle scatter
plot(1:length(coneAngle),coneAngle,'kx')
hold on
plot([1,length(coneAngle)],MeanCone.*[1,1],'r-')
plot([1,length(coneAngle)],(MeanCone+stdCone).*[1,1],'m--')
plot([1,length(coneAngle)],(MeanCone-stdCone).*[1,1],'m--')
title('VOR9 CONE ANGLE Scatter Plot')
xlabel('Frame Number')
ylabel('Cone Angle / deg')
hold off


% Combined Welch Analysis
% WelchMean = mean(WelchSolve,2); % average across rows
% WelchMeanSmooth = mean(WelchSolveSmooth,2); % average across rows
% figure()
% subplot(2,1,1)
% plot(10*log10(WelchMean))
% subplot(2,1,2)
% plot(10*log10(WelchMeanSmooth))

fprintf('Number of frames analysed: %g.\n',n_images)
toc