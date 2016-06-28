%% Supplemental figures
figDir = ['/Users/abock/Dropbox-Aguirre-Brainard-Lab/bock/bock_aguirre_manuscripts/' ...
    'retinotopy_templates/Figures/Figure7_Suppl/raw'];
sessions = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014/' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014/' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
subNames = {'AEK' 'ASB' 'GKA'};
subRuns = {...
    [3,4,6] ...
    [2,4,6] ...
    [2,4,6] ...
    };
threshs = {...
    [0.648624 0.619279] ...
    [0.916082 0.871358] ...
    [0.903127 0.863659] ...
    }; % 0.9 quantile for pRF correlation
template = 'coarse';
hemis = {'lh' 'rh'};

src_surf = '0.1.inflated';
trg_surf = 'inflated';
func = 's5.wdrf.tf';
cond = 'Movie';
templateSize = 'full';
%% Plot the peak HRFs
% Mean time-to-peak of the HRF models found by the pRF approach in each of
%   the visual areas (V1, V2, V3) for each of the 6 figures
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        rawHRF = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.' func '.cortex.copeakt.avg.prfs.nii.gz']));
        rawco = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.' func '.cortex.co.avg.prfs.nii.gz']));
        rawind = rawco.vol>threshs{ss}(hh); % 0.9 quantile threshold for pRF data
        pRFareas = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            [hemi '.areas.pRF.nii.gz']));
        V1ind = abs(pRFareas.vol)==1 & rawind;
        V2ind = abs(pRFareas.vol)==2 & rawind;
        V3ind = abs(pRFareas.vol)==3 & rawind;
        peak.V1(ss,hh) = mean(rawHRF.vol(V1ind));
        peak.V2(ss,hh) = mean(rawHRF.vol(V2ind));
        peak.V3(ss,hh) = mean(rawHRF.vol(V3ind));
    end
end
HRFpeak(1) = mean(peak.V1(:));
HRFstd(1) = std(peak.V1(:));
HRFpeak(2) = mean(peak.V2(:));
HRFstd(2) = std(peak.V2(:));
HRFpeak(3) = mean(peak.V3(:));
HRFstd(3) = std(peak.V3(:));
fullFigure;
errorbar(HRFpeak,HRFstd);
xlim([0 4]);
set(gca,'XTick',[0 1 2 3 4],'XTickLabel',{'' 'V1' 'V2' 'V3' ''},'FontSize',15);
xlabel('Visual Area','FontSize',20);
ylabel('HRF peak (seconds)','FontSize',20);
ylim([0 10]);
axis square;
cd(figDir);
savefigs('pdf','HRFpeak');
close all;
p = anova1([peak.V1(:),peak.V2(:),peak.V3(:)]);
disp(['p = ' num2str(p)]);
disp(['means (V1, V2, V3) = ' num2str(mean([peak.V1(:),peak.V2(:),peak.V3(:)]))]);
disp(['stds (V1, V2, V3) = ' num2str(std([peak.V1(:),peak.V2(:),peak.V3(:)]))]);
%% Movie stimulus correlation
% Load movie and set default variables
inputMovie = ['/Users/abock/Dropbox-Aguirre-Brainard-Lab/bock/bock_aguirre_manuscripts/' ...
    'retinotopy_templates/Figures/Figure7_Suppl/ALL_0_600_TRs.mat']; %movie must be grayscale
load(inputMovie); %mov
timeStamps = size(mov,3); % get TRs
downRes = round(([size(mov,1) size(mov,2)])/10); % downsample movie
sampledMovie = double(imresize(mov,[downRes(1) downRes(2)])); % downsample movie
clear mov;
%% Find distances to other pixels
distMat = zeros(downRes(1)*downRes(2));
location = zeros(downRes(1)*downRes(2),2);
xPos = 1:downRes(2);
yPos = 1:downRes(1);
ct = 0;
% Find x, y location for each pixel
for x = 1:length(xPos);
    for y = 1:length(yPos)
        ct = ct+1;
        location(ct,:) = [x y];
    end
end
% Calculate distance from each pixel to all other pixels
progBar = ProgressBar(length(distMat),'Calculating distances...');
for d = 1:length(distMat);
    distMat(d,:) = sqrt(...
        abs( (location(d,1) - location(:,1)).^2 ) + ...
        abs( (location(d,2) - location(:,2)).^2 )...
        );
    if ~mod(d,1000);progBar(d);end
end
%% Create correlation matrix
% Create 2D movie
twoDmov = (reshape(sampledMovie,[downRes(1)*downRes(2) 300]))';
corrMat = zeros(size(twoDmov,2));
[vals,Tasks] = break_up_matrix(length(twoDmov),1000);
progBar = ProgressBar(Tasks,'Creating correlation matrix...');
for p = 1:Tasks
    corrMat(vals{p},:) = corr(twoDmov(:,vals{p}),twoDmov);
    progBar(p);
end
%% Flatten data and sort
flatcorr = reshape(corrMat,[size(corrMat,1)*size(corrMat,2) 1]);
flatdist = reshape(distMat,[size(distMat,1)*size(distMat,2) 1]);
% sort the values (this takes a couple minutes)
[sortdist,ind] = sort(flatdist);
sortcorr = flatcorr(ind);
clear flatcorr flatdist
%% Bin the distance values, average correlation for these bins
[vals,Tasks] = break_up_matrix(length(sortdist),100000);
mean_dists = zeros(Tasks,1);
mean_corrs = zeros(Tasks,1);
progBar = ProgressBar(Tasks,'Computing means...');
for i = 1:Tasks
    mean_dists(i) = mean(sortdist(vals{i}));
    mean_corrs(i) = mean(sortcorr(vals{i}));
   if ~mod(i,100);progBar(i);end
end
%% Plot correlation by distance
cd(figDir);
fullFigure;plot(mean_dists*10,mean_corrs,'.');
xlim([0 max(mean_dists*10)]);
xlabel('Distance (pixels)','FontSize',20);
ylabel('Correlation (Pearson''s r)','FontSize',20);
title('Spatial correlation of movie stimulus','FontSize',20);
axis square;
savefigs('pdf','Movie_stimulus_correlation');
close all;