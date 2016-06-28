%% Setup the directories and variables
dataDir = '/data/jet/abock/data';
sessions = {...
    fullfile(dataDir,'Retinotopy_Templates','AEK','10012014') ...
    fullfile(dataDir,'Retinotopy_Templates','ASB','10272014') ...
    fullfile(dataDir,'Retinotopy_Templates','GKA','10152014') ...
    };
hemis = {'lh' 'rh'};
badAreas = [4,5]; % visual areas to exclude;
eLims = [1 5]; % eccentricity lims [low high];
saveDir = ['/Users/abock/Dropbox (Aguirre-Brainard Lab)/bock/Bock_Aguirre_manuscripts/' ...
    'Retinotopy_Templates/Figures/Figure6_bar_plots_error/raw'];
cd(saveDir);
%% Compute error in split-halves
func = 's5.wdrf.tf';
splitComb = combnk(1:6,3);
ct = 0;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    pRFDir = fullfile(session_dir,'pRFs','pRF_templates');
    for hh = 1:length(hemis)
        ct = ct + 1;
        hemi = hemis{hh};
        Ecc = load_nifti(fullfile(pRFDir,[hemi '.ecc.pRF.nii.gz']));
        Areas = load_nifti(fullfile(pRFDir,[hemi '.areas.pRF.nii.gz']));
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
        clear tmp
        for i = 1:length(splitComb)/2
            inEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            inPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(i,:),'%1d') '.prfs.nii.gz']);
            tempEcc = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.coecc.avg.' num2str(splitComb(end-(i-1),:),'%1d') '.prfs.nii.gz']);
            tempPol = fullfile(session_dir,'pRFs',...
                [hemi '.' func '.cortex.copol.avg.' num2str(splitComb(end-(i-1),:),'%1d') '.prfs.nii.gz']);
            [degerror] = computepRFerror(inEcc,inPol,tempEcc,tempPol,verts);
            tmp(i) = nanmean(degerror);
        end
        splitError.mean(ct) = mean(tmp(:));
    end
end
%% Compute error in anat templates (deformed)
func = 's5.wdrf.tf';
splitComb = combnk(1:6,3);
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
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
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
        anatsplitError.mean(ct) = mean(tmp(:));
    end
end
%% Compute error in 'best' coarse templates
templateType = 'coarse';
func = 's5.wdrf.tf';
fitType = 'V2V3';
splitComb = combnk(1:6,3);
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
        verts = Ecc.vol<eLims(2) & Ecc.vol>eLims(1) & ~ismember(abs(Areas.vol),badAreas);
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
        coarsesplitError.mean(ct) = mean(tmp(:));
    end
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
%% Plot error
dim = [.5 .75 .1 .1];
split.mean = mean(splitError.mean);
split.std = std(splitError.mean);
anat.mean = mean(anatsplitError.mean);
anat.std = std(anatsplitError.mean);
coarse.mean = mean(coarsesplitError.mean);
coarse.std = std(coarsesplitError.mean);
nodanat.mean = mean(nondanatsplitError.mean);
nodanat.std = std(nondanatsplitError.mean);
fullFigure;bar([split.mean anat.mean coarse.mean nodanat.mean]);
hold on;
errorbar([split.mean anat.mean coarse.mean nodanat.mean],...
    [split.std anat.std coarse.std nodanat.std],'.','MarkerSize',0.01);
xlabel('Map Type','FontSize',20);
ylabel('Mean error (degress visual angle)','FontSize',20);
set(gca,'XTickLabel',{'split-half','anat-deformed','coarse','anat-non-deformed'},'FontSize',15);
axis square
annotation('textbox',dim,'String','error bars = SD','FontSize',15,'HorizontalAlignment','center');
title('Error in retinotopy prediction','FontSize',20);
%savefigs('pdf','templates_error');
%close all;