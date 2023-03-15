% clear all; close all; clc;

% ***** 02_SH_2 *************************
% AUTHOR:       Zachary-Jacques Gray
% DATE CREATED: 14 Jan 2022
% DESCRIPTION:  MATLAB script to detect streams and droplets from image via
%               clustering and edge-detection algorithms.
%               Analyzes for injector flow qualities, e.g. atomisation and
%               oscillations.
% REQUIREMENTS: MATLAB Image Processing Toolbox
%               MATLAB Signal Processing Toolbox

tic;
startTime = clock;
fprintf('\nTime Started: %g:%g.%0.2g\n',startTime(4),startTime(5),startTime(6))
% BACKGROUND ---------------------------//
imgB = imread('02_VOR_9_20211217_1547010001.tif');
greyB = im2gray(imgB);
L = 2^16; % 16-bit Max Intensity
greyB = (L - 1) - greyB; % negative
% figure()
% imshow(greyB);
% title('VOR 0 BACKGROUND');


% LOOP FOR EACH IMAGE ===========================================//
n_images = 6400; % number of frames to run starting from image 1

xrange = 501;

% Estimated time of completion
if exist('runTime','var')
    estTime = runTime*n_images/runFrame;
else
    estTime = 7.8*n_images/100;
end
estMin = startTime(5)+round((startTime(6)+estTime)/60);
if estMin >= 60
    estHour = startTime(4) + 1;
    estMin = estMin - 60;
else
    estHour = startTime(4);
end
fprintf('Estimated Completion: %g:%g\n',estHour,estMin)


WelchSolve = zeros(129,n_images); % may be different for each new data set
WelchSolveSmooth = WelchSolve; % may be different for each new data set
coneAngles = zeros(1,n_images);
amp = zeros(1,n_images);
freq = zeros(1,n_images);
yWelch = zeros(xrange,n_images); % 551 is xvalue range given by m45 further down
yWelchT = zeros(xrange,n_images); % 551 is xvalue range given by m45 further down

for frame = 1:n_images
    imgstr = '02_VOR_9T_20211217_155036';
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
        warning('Can only scan between 1 to 9,999 images')
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
Z = imadjust(Z,[0.025,0.7],[0,1],0.7);  % gamma and contrast correction for viewing
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
rotateAngle = 120; % degrees
Cluster45 = imrotate(Cluster, rotateAngle,'bilinear','crop');
rotateAngle = 40;

% Edge Detection
CannyEdge = edge(Cluster45,'canny');
% figure()
% imshow(CannyEdge);


% FOURIER TRANSFORM ---------------------//
% Welch's PSD Estimate > > >

% First convert segment of stream into (x,y)-function
m45 = [500,1000]; % first and last x-range
x = m45(1):m45(2); % initialise x-range
x = x - m45(1);
y = zeros(size(x)); % initialise y-range
for i = m45(1):m45(2)
    for j = 550:750 % top range
        if CannyEdge(j,i) == 1 % if line exists at that pixel
            if y(i-m45(1)+1) == 0
                y(i-m45(1)+1) = j; %bottom most part of line (most likely stream)
            end
        end
    end
end

yWelchT(:,frame) = -flip(y);
y = -y + mean(y); % flip to match same way as edge

% figure
% plot(x,y)

% Perform Sine FFT Fitting
% SineParameters = sineFit(x,y);
% amp(frame) = SineParameters(2);
% freq(frame) = SineParameters(3);


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


% Prepare y for Welch/FFT
y = detrend(y); % remove any constant trend
y = flip(y); % define y(0) -> y(end) as downstream direction
yWelch(:,frame) = y;


end 
% *******************// image loop end //***********************

%% COMBINED DATA ANALYSIS ----------------//
fprintf('\n> > > VOR 9_2 > > >\n')

% Import Data of past results
warning('You are loading previously saved results!')
load VOR9_2_results.mat

% Contour of Edge Data
tsmall = 1:100; xvect = 1:xrange; 
[Time,Xs] = meshgrid(tsmall,xvect);
figure()
contourf(Time,Xs,yWelch(:,1:100))
xlabel('Time'); ylabel('X');
caxis([-50,50])
figure()
[Time,Xs] = meshgrid(tsmall,xvect(1:3:end));
surf(Time,Xs,yWelch(1:3:end,1:100))
xlabel('Time'); ylabel('X');
caxis([-50,50])

% WELCH ANALYSIS
yWelchTime = yWelchT - mean(yWelchT,2);
yWelchTime = transpose(detrend(yWelchTime'));
tvect = 1:n_images; xvect = 1:xrange;
window = 200;

[ptt,ft] = pwelch(yWelchTime,window,window/2,[]); 
[pxx,fx] = pwelch(yWelch',window,window/2,[]);

Fs = 3200; % sampling rate
% wave number axis
fx = fx/pi; % /pixel
% frequency axis
ft = ft.*Fs./(2*pi); % Hz

figure()
[Wave,X] = meshgrid(fx,xvect);
contourf(1./Wave,X,log10(pxx'));
title('Welch Power Spectral Density - Wavelength')
xlabel('Wavelength / pixels')
ylabel('Stream Displacement / pixels')
zlabel('log_{10}(PSD)')
colorbar

figure()
[Freq,T] = meshgrid(ft,tvect);
contourf(log10(Freq),T./Fs,log10(ptt'),'LineColor','none');
title('Welch Power Spectral Density - Frequency')
xlabel('log_{10}(Frequency / Hz)')
ylabel('Time / s')
zlabel('log_{10}(PSD)')
colorbar

figure()
[Wave,X] = meshgrid(fx,xvect);
surf(1./Wave,X,log10(pxx'));
title('Welch Power Spectral Density - Wavelength')
xlabel('Wavelength / pixels')
ylabel('Stream Displacement / pixels')
zlabel('log_{10}(PSD)')

figure()
[Freq,T] = meshgrid(ft,tvect(10:10:end));
surf(log10(Freq),T./Fs,log10(transpose(ptt(:,10:10:end))));
title('Welch Power Spectral Density - Frequency')
xlabel('log_{10}(Frequency / Hz)')
ylabel('Time / s')
zlabel('log_{10}(PSD)')

% MEAN PLOTS
% figure()
% plot(fx,10*log10(mean(pxx,2)))
% figure()
% plot(ft,10*log10(mean(ptt,2)))

% Mean Frequency
figure()
intFreq = trapz(tvect, ptt'); 
plot(ft,intFreq);
[MaxFreq, FreqIndex] = max(intFreq);
FreqPSD = ft(FreqIndex);
title('Welch Power Spectral Density - Mean Frequency')
xlabel('Frequency / Hz)')
ylabel('PSD')
xlim([0 120])

% Mean Wavelength
figure()
intWave = trapz(xvect, pxx'); 
semilogx(1./fx,intWave);
[MaxWave, WaveIndex] = max(intWave);
WavePSD = 1./fx(WaveIndex);
title('Welch Power Spectral Density - Mean Wavelength')
xlabel('Wavelength / pixels')
ylabel('PSD')

% Amplitude
intRMS = trapz(fx, sqrt(pxx)); 
figure()
intAmp = intRMS./0.7071; % RMS to amplitude calculation
plot(xvect,intAmp);
AmpPSD = mean(rmoutliers(intAmp));
title('Welch Power Spectral Density - Amplitude')
xlabel('Stream Displacement / pixels')
ylabel('Amplitude / pixels')

% Print to command window
fprintf('Amplitude: %g pts \tFrequency: %g Hz \tWavelength: %g pixels\n----//----\n',AmpPSD,FreqPSD,WavePSD)
% fprintf('\nSTDEV: %g pts \tSTDEV: %g Hz\n----//----\n',stdAmp,stdFreq)


% % FREQUENCY AND AMPLITUDE
% % for use with Sine Fitting
% freqNormal = freq .*Fs./(2*pi);
% amplitude = rmoutliers(amp);
% frequency = rmoutliers(freqNormal);
% MeanAmp = mean(amplitude);
% MeanFreq = mean(frequency);
% stdAmp = std(amplitude);
% stdFreq = std(frequency);
% 
% figure % amplitude scatter
% plot(1:length(amplitude),amplitude,'kx')
% hold on
% plot([1,length(amplitude)],MeanAmp.*[1,1],'r-')
% plot([1,length(amplitude)],(MeanAmp+stdAmp).*[1,1],'m--')
% plot([1,length(amplitude)],(MeanAmp-stdAmp).*[1,1],'m--')
% title('SH STREAM OSCILLATION AMPLITUDE')
% xlabel('Frame Number')
% ylabel('Amplitude / pts')
% hold off
% 
% figure % frequency scatter
% plot(1:length(frequency),frequency,'kx')
% hold on
% plot([1,length(frequency)],MeanFreq.*[1,1],'r-')
% plot([1,length(frequency)],(MeanFreq+stdFreq).*[1,1],'m--')
% plot([1,length(frequency)],(MeanFreq-stdFreq).*[1,1],'m--')
% title('SH STREAM OSCILLATION FREQUENCY')
% xlabel('Frame Number')
% ylabel('Frequency / Hz')
% hold off

% CONE ANGLE
coneAngle = rmoutliers(-coneAngles);
MeanCone = mean(coneAngle);
stdCone = std(coneAngle);
fprintf('Cone Angle: %.3g deg\tSTDEV: %g deg\n----//----\n',MeanCone,stdCone)

figure % cone angle scatter
plot(1:length(coneAngle),coneAngle,'kx')
hold on
plot([1,length(coneAngle)],MeanCone.*[1,1],'r-')
plot([1,length(coneAngle)],(MeanCone+stdCone).*[1,1],'m--')
plot([1,length(coneAngle)],(MeanCone-stdCone).*[1,1],'m--')
title('SH CONE ANGLE')
xlabel('Frame Number')
ylabel('Cone Angle / deg')
hold off


fprintf('Number of frames analysed: %g.\n',n_images)
toc
runTime = toc;
runFrame = n_images;