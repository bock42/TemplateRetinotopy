%% Supplemental figures
figDir = ['/Users/abock/Dropbox-Aguirre-Brainard-Lab/bock/bock_aguirre_manuscripts/' ...
    'retinotopy_templates/Figures/Figure6_Suppl/raw'];
sessions = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014/' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014/' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014/' ...
    };
threshs = {...
    [0.648624 0.619279] ...
    [0.916082 0.871358] ...
    [0.903127 0.863659] ...
    }; % 0.9 quantile for pRF correlation
hemis = {'lh' 'rh'};
func = 's5.wdrf.tf';
% Error in Visual angle
templateType = 'coarse';
eLims = [1 5]; % eccentricity lims [low high];
fitType = 'V2V3';
badAreas = {...
    [2 3 4 5] ...
    [1 3 4 5] ...
    [1 2 4 5]
    };
visualAreas = {'V1' 'V2' 'V3'};
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
% Find distances to other pixels
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
% Create correlation matrix
% Create 2D movie
twoDmov = (reshape(sampledMovie,[downRes(1)*downRes(2) 300]))';
corrMat = zeros(size(twoDmov,2));
[vals,Tasks] = break_up_matrix(length(twoDmov),1000);
progBar = ProgressBar(Tasks,'Creating correlation matrix...');
for p = 1:Tasks
    corrMat(vals{p},:) = corr(twoDmov(:,vals{p}),twoDmov);
    progBar(p);
end
% Flatten data and sort
flatcorr = reshape(corrMat,[size(corrMat,1)*size(corrMat,2) 1]);
flatdist = reshape(distMat,[size(distMat,1)*size(distMat,2) 1]);
% sort the values (this takes a couple minutes)
[sortdist,ind] = sort(flatdist);
sortcorr = flatcorr(ind);
clear flatcorr flatdist
% Bin the distance values, average correlation for these bins
[vals,Tasks] = break_up_matrix(length(sortdist),100000);
mean_dists = zeros(Tasks,1);
mean_corrs = zeros(Tasks,1);
progBar = ProgressBar(Tasks,'Computing means...');
for i = 1:Tasks
    mean_dists(i) = mean(sortdist(vals{i}));
    mean_corrs(i) = mean(sortcorr(vals{i}));
    if ~mod(i,100);progBar(i);end
end
% Plot correlation by distance
cd(figDir);
fullFigure;plot(mean_dists*10,mean_corrs,'.');
xlim([0 max(mean_dists*10)]);
xlabel('Distance (pixels)','FontSize',20);
ylabel('Correlation (Pearson''s r)','FontSize',20);
title('Spatial correlation of movie stimulus','FontSize',20);
axis square;
savefigs('pdf','Movie_stimulus_correlation');
close all;
%% Plot the error in visual angle separately for V1, V2, and V2 (similar to Figure 5)

%% Compute split error
splitComb = combnk(1:6,3);
splitError = nan(length(visualAreas),length(sessions),63,length(hemis));
progBar = ProgressBar(length(visualAreas),'computing error...');
for vv = 1:length(visualAreas)
    for ss = 1:length(sessions)
        session_dir = sessions{ss};
        pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
        modelTdir = fullfile(session_dir,'pRFs','coarse_model_templates');
        aDirs = fullfile(session_dir,'pRFs',templateType,func);
        pDirs = listdir(fullfile(aDirs,'Movie_half*'),'dirs');
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
            Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
            verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas{vv});
            for p = 1:length(pDirs)
                tDir = fullfile(aDirs,pDirs{p});
                [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,[],[],[],fitType);
                dotinds = strfind(sorted_templates{1},'.');
                tempEcc = fullfile(modelTdir,...
                    [hemi '.ecc.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
                tempPol = fullfile(modelTdir,...
                    [hemi '.pol.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
                clear tmp
                for i = 1:length(splitComb)
                    inEcc = fullfile(session_dir,'pRFs',...
                        [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                    inPol = fullfile(session_dir,'pRFs',...
                        [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                    [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
                    tmp(i) = nanmean(degerror);
                end
                splitError(vv,ss,p,hh) = mean(tmp(:));
            end
        end
    end
    progBar(vv);
end
save(fullfile(figDir,'splitError.mat'),'splitError');
%% Compute error in anat templates (deformed)
splitComb = combnk(1:6,3);
progBar = ProgressBar(length(visualAreas),'computing error...');
for vv = 1:length(visualAreas)
    ct = 0;
    for ss = 1:length(sessions)
        session_dir = sessions{ss};
        pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
        modelTdir = fullfile(session_dir,'pRFs','anat_templates');
        for hh = 1:length(hemis)
            ct = ct + 1;
            hemi = hemis{hh};
            Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
            Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
            verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas{vv});
            tempEcc = fullfile(modelTdir,...
                [hemi '.ecc.anat.nii.gz']);
            tempPol = fullfile(modelTdir,...
                [hemi '.pol.anat.nii.gz']);
            clear tmp
            for i = 1:length(splitComb)
                inEcc = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                inPol = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
                tmp(i) = nanmean(degerror);
            end
            anatsplitError(vv).mean(ct) = mean(tmp(:));
        end
    end
    progBar(vv);
end
%% Compute error in 'best' coarse templates
splitComb = combnk(1:6,3);
progBar = ProgressBar(length(visualAreas),'computing error...');
for vv = 1:length(visualAreas)
    ct = 0;
    for ss = 1:length(sessions)
        session_dir = sessions{ss};
        pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
        tDir = fullfile(session_dir,'pRFs',templateType,func,'Movie',fitType);
        modelTdir = fullfile(session_dir,'pRFs','coarse_model_templates');
        for hh = 1:length(hemis)
            ct = ct + 1;
            hemi = hemis{hh};
            Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
            Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
            verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas{vv});
            [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,[],[],[],fitType);
            dotinds = strfind(sorted_templates{1},'.');
            tempEcc = fullfile(modelTdir,...
                [hemi '.ecc.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
            tempPol = fullfile(modelTdir,...
                [hemi '.pol.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']);
            clear tmp
            for i = 1:length(splitComb)
                inEcc = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                inPol = fullfile(session_dir,'pRFs',...
                    [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
                [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
                tmp(i) = nanmean(degerror);
            end
            coarsesplitError(vv).mean(ct) = mean(tmp(:));
        end
    end
    progBar(vv);
end
%% Compute error in anat templates (non-deformed)
func = 's5.wdrf.tf';
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    modelTdir = fullfile(session_dir,'pRFs','coarse_model_templates');
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        tempEcc = fullfile(modelTdir,...
            [hemi '.ecc.4.4.4.nii.gz']);
        tempPol = fullfile(modelTdir,...
            [hemi '.pol.4.4.4.nii.gz']);
        clear tmp
        for i = 1:length(splitComb)
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        nondanatsplitError.mean(ct) = mean(tmp(:));
    end
end