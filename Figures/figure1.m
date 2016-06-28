% Example images of the retinotopic mapping data and pixar movies.
%   These images were taken as screen shots from the code used to produce
%   them.
%
%   Also takes an example session directory and run, and plots timecourse from a
%   single vertex in V1, V2, V3.
%
%   Currently, the subject is GKA, the run is #4 (middle movie run of the
%   session), and the template is the pRF template ('gold standard');
%
%   Vertices were chosen in V1, V2, V3 that were the closest in visual
%   space, near the fovea, and in the dorsal areas (ventral visual field)
%
%   V1 vertex - V1ecc = 3.3857; V1pol = 104.0787;
%   V2 vertex - V2ecc = 2.9540; V2pol = 110.4599;
%   V3 vertex - V3ecc = 3.6613; V3pol = 106.6110;

%% Set defaults
outDir = ['/Users/abock/Dropbox (Aguirre-Brainard Lab)/Bock/'...
    'Bock_Aguirre_manuscripts/Retinotopy_Templates/Figures/' ...
    'Figure1_experimental_overview/raw'];
cd(outDir);
tDir = '~/data';
subject_name = 'fsaverage';
dataEcc = fullfile(tDir,'2014-10-29.eccen-template.nii.gz');
dataPol = fullfile(tDir,'2014-10-29.angle-template.RADS.nii.gz');
dataAreas = fullfile(tDir,'2014-10-29.areas-template.nii.gz');
session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014';
hemi = 'lh';
template='pRF';
decimation_level='0.1';
src_surf = 'inflated';
trg_surf = [decimation_level '.' src_surf];
func = 's5.dbrf.tf.surf';
cond = 'Movie';
runNum = 4;
%% Play movie for retinotopic mapping
%   screenshots were taken to create example images
images = make_bars;
play_pRF_movie('foo',1,images);
%% Play pixar movies
%   screenshots were taken to create example images
block=0;
fix=0;
subj='foo';
run=999;
fromTime=0;
nBlocks=1;
blockDur=600;
soundVol=0;
indexisFrames=0;
moviename=['/Users/abock/Dropbox-Aguirre-Brainard-Lab/bock/Bock_Aguirre_manuscripts/' ...
    'Retinotopy_Templates/Figures/Figure1_experimental_overview/raw/ALL.mov'];
playmovie_block(block,fix,subj,run,fromTime,nBlocks,blockDur,...
    soundVol,indexisFrames,moviename);
%% Convert movie images to surface images
pRF_images = {'pRF_1' 'pRF_2' 'pRF_3'};
pixar_images = {'pixar_1' 'pixar_2' 'pixar_3'};
centerPixs={[1050 1680] [525 2520] [1575 2520] [525 840] [1575 840]};
centerNames = {'center' 'upper_right' 'lower_right' 'upper_left' 'lower_left'};
tmp = load_nifti(dataAreas);
thresh = abs(tmp.vol)<=3;
% Ecc
surface_plot('ecc',dataEcc,subject_name,'lh','inflated',thresh,0.85,'hemi');
savename = 'fsaverage_ecc';
savefigs('png',savename);
close all;
% Pol
surface_plot('pol',dataPol,subject_name,'lh','inflated',thresh,0.85,'hemi');
savename = 'fsaverage_pol';
savefigs('png',savename);
close all;
% Areas
surface_plot('blueareas',dataAreas,subject_name,'lh','inflated',thresh,0.85,'hemi');
savename = 'fsaverage_areas';
savefigs('png',savename);
close all;
for i = 1:length(pRF_images);
    % pRF
    inputImg = fullfile(outDir,[pRF_images{i} '.png']);
    outputNifti = fullfile('~',[pRF_images{i} '.nii.gz']);
    convert_image2surf(inputImg,outputNifti,dataEcc,dataPol,dataAreas);
    surface_plot('image',outputNifti,subject_name,'lh','inflated',thresh);
    savename = ['fsaverage_' pRF_images{i}];
    savefigs('png',savename);
    close all;
    system(['rm ' outputNifti]);
    for j = 1:length(centerPixs)
        centerPix = centerPixs{j};
        % pixar movie
        inputImg = fullfile(outDir,[pixar_images{i} '.png']);
        outputImg = fullfile(outDir,[pixar_images{i} '_' centerNames{j} '.png']);
        save_centerPix_images(inputImg,outputImg,centerPix);
        outputNifti = fullfile('~',[pixar_images{i} '.nii.gz']);
        convert_image2surf(inputImg,outputNifti,dataEcc,dataPol,dataAreas,centerPix);
        surface_plot('image',outputNifti,subject_name,'lh','inflated',thresh);
        savename = ['fsaverage_' pixar_images{i} '_' centerNames{j}];
        savefigs('png',savename);
        close all;
        system(['rm ' outputNifti]);
    end
end
%% Plot example timecourse from V1,V2,V3 from an example subject
%   Takes an example session directory and run, and plots timecourse from a
%   single vertex in V1, V2, V3.
%
%   Currently, the subject is GKA, the run is #4 (middle run of the
%   session), and the template is the pRF template ('gold standard');
%
% Find bold run directories
d = find_bold(session_dir);
% Decimated pRF templates
switch template
    case 'V1'
        % need to add this, currently don't have a decimated directory
    case 'pRF'
        areas = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.areas.pRF.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.ecc.pRF.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','pRF_templates',...
            'decimated_templates',[hemi '.pol.pRF.nii.gz']));
    case 'fine'
        tdir = fullfile(session_dir,'pRFs','fine',func(1:end-5),cond);
        [~,~,sorted_templates] = find_best_template(template,tdir,hemi);
        bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
        areas = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.areas.' bestTemplate '.nii.gz']));
        ecc = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.ecc.' bestTemplate '.nii.gz']));
        pol = load_nifti(fullfile(session_dir,'pRFs','fine_model_templates',...
            'decimated_templates',[hemi '.pol.' bestTemplate '.nii.gz']));
end
% Get target indices, ecc, pol for V1
V1trgind = find(abs(areas.vol)==1);
V1trgecc = ecc.vol(V1trgind);
V1trgpol = pol.vol(V1trgind);
% Get target indices, ecc, pol for V2
V2trgind = find(abs(areas.vol)==2);
V2trgecc = ecc.vol(V2trgind);
V2trgpol = pol.vol(V2trgind);
% Get target indices, ecc, pol for V3
V3trgind = find(abs(areas.vol)==3);
V3trgecc = ecc.vol(V3trgind);
V3trgpol = pol.vol(V3trgind);
% Decimate target file
in_vol = fullfile(session_dir,d{runNum},[func '.' hemi '.nii.gz']);
out_vol = fullfile(session_dir,d{runNum},[func '.' trg_surf '.' hemi '.nii.gz']);
if ~exist(out_vol,'file')
    decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol);
end
% Load target surface file
trgfile = out_vol;
% Get distances in target ROI in visual space (from V1 centers)
% V1
[V1trgX,V1trgY] = pol2cart(V1trgpol,V1trgecc);
for i = 1:length(V1trgind)
    V1trgdists(i,:) = sqrt( (V1trgX(i) - V1trgX).^2 + (V1trgY(i) - V1trgY).^2);
end
% V2
[V2trgX,V2trgY] = pol2cart(V2trgpol,V2trgecc);
for i = 1:length(V2trgind)
    V2trgdists(i,:) = sqrt( (V2trgX(i) - V1trgX).^2 + (V2trgY(i) - V1trgY).^2);
end
% V3
[V3trgX,V3trgY] = pol2cart(V3trgpol,V3trgecc);
for i = 1:length(V3trgind)
    V3trgdists(i,:) = sqrt( (V3trgX(i) - V1trgX).^2 + (V3trgY(i) - V1trgY).^2);
end
% Load timecourses
disp('Loading timecourses...');
% load files
trg = load_nifti(trgfile);
trgdims = size(trg.vol);
trgtc = reshape(trg.vol,trgdims(1)*trgdims(2)*trgdims(3),trgdims(4))';
% Pull out relevant timecourses (V1)
V1trgtc = trgtc(:,V1trgind);
V1trgtc = (V1trgtc ./ repmat(mean(V1trgtc),size(V1trgtc,1),1))*100; % scale the data
V1trgtc = V1trgtc - repmat(mean(V1trgtc),size(V1trgtc,1),1); % convert to percent signal change
% Pull out relevant timecourses (V2)
V2trgtc = trgtc(:,V2trgind);
V2trgtc = (V2trgtc ./ repmat(mean(V2trgtc),size(V2trgtc,1),1))*100; % scale the data
V2trgtc = V2trgtc - repmat(mean(V2trgtc),size(V2trgtc,1),1); % convert to percent signal change
% Pull out relevant timecourses (V3)
V3trgtc = trgtc(:,V3trgind);
V3trgtc = (V3trgtc ./ repmat(mean(V3trgtc),size(V3trgtc,1),1))*100; % scale the data
V3trgtc = V3trgtc - repmat(mean(V3trgtc),size(V3trgtc,1),1); % convert to percent signal change
% Set timecourses with very little variation (var<.1) to flat
V1trgtc = set_to_flat(V1trgtc);
V2trgtc = set_to_flat(V2trgtc);
V3trgtc = set_to_flat(V3trgtc);
disp('done.');
%% Plot and save timecourses
V1vtx = 468;
[~,V1ind] = min(V1trgdists(:,V1vtx));
[~,V2ind] = min(V2trgdists(:,V1vtx));
[~,V3ind] = min(V3trgdists(:,V1vtx));
V1ecc = V1trgecc(V1ind);
V2ecc = V2trgecc(V2ind);
V3ecc = V3trgecc(V3ind);
V1pol = rad2deg(V1trgpol(V1ind)) + 90;
V2pol = rad2deg(V2trgpol(V2ind)) + 90;
V3pol = rad2deg(V3trgpol(V3ind)) + 90;
% dashed line
dashedLine = zeros(size(V1trgtc(:,V1ind)));
dashedLine(1:10:end) = nan;
figure('units','normalized','position',[0 0 1 1]);
plot(V1trgtc(:,V1ind)); hold on;
plot(dashedLine,'k','LineWidth',1);
ylim([-10 10]);
title(['V1ecc = ' num2str(V1ecc) '; V1pol = ' num2str(V1pol)],'FontSize',20);
xlabel('TR','FontSize',20);
ylabel('Percentage Signal Change','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf','V1tc_movie');
close all
figure('units','normalized','position',[0 0 1 1]);
plot(V2trgtc(:,V2ind)); hold on;
plot(dashedLine,'k','LineWidth',1);
ylim([-10 10]);
title(['V2ecc = ' num2str(V2ecc) '; V2pol = ' num2str(V2pol)],'FontSize',20);
xlabel('TR','FontSize',20);
ylabel('Percentage Signal Change','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf','V2tc_movie');
close all
figure('units','normalized','position',[0 0 1 1]);
plot(V3trgtc(:,V3ind)); hold on;
plot(dashedLine,'k','LineWidth',1);
ylim([-10 10]);
title(['V3ecc = ' num2str(V3ecc) '; V3pol = ' num2str(V3pol)],'FontSize',20);
xlabel('TR','FontSize',20);
ylabel('Percentage Signal Change','FontSize',20);
savefigs('pdf','V3tc_movie');
set(gca,'FontSize',15);
close all
% Plot together
figure('units','normalized','position',[0 0 1 1]);
plot(V1trgtc(:,V1ind),'b'); hold on;
plot(V2trgtc(:,V2ind),'g');
plot(V3trgtc(:,V3ind),'r');
plot(dashedLine,'k','LineWidth',1);
ylim([-10 10]);
[r] = ICC(zscore([V1trgtc(:,V1ind),V2trgtc(:,V2ind),V3trgtc(:,V3ind)]),'1-1');
title([...
    'V1-V2 R = ' num2str(corr2(V1trgtc(:,V1ind),V2trgtc(:,V2ind))) '; ' ...
    'V1-V3 R = ' num2str(corr2(V1trgtc(:,V1ind),V3trgtc(:,V3ind))) '; ' ...
    'V2-V3 R = ' num2str(corr2(V2trgtc(:,V2ind),V3trgtc(:,V3ind))) '; ' ...
    'ICC = ' num2str(r)],'FontSize',20);
xlabel('TR');
ylabel('Percentage Signal Change');
legend({'V1' 'V2' 'V3'},'FontSize',15);
xlabel('TR','FontSize',20);
ylabel('Percentage Signal Change','FontSize',20);
set(gca,'FontSize',15);
savefigs('pdf','V1_V2_V3_tc_movie');
close all
%% Correlation matrix
% The correlation matrix in B is produced in the Figure 2 code, copy from
% that directory