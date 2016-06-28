%% Error plots
% Map that shows the Sum-squared error in the residuals for fitting the
%   movie data using the best fine-scale search Shira model. This can be
%   done for each subject, but we should also create an average image
%   across hemispheres and subjects.
%
% Map that shows the error in prediction of polar angle and eccentricity
%   That is, a comparison of the best fine-scale Schira model and the
%   "gold standard" deformed template to the pRF data. This can be created
%   for each subject individually and also for an average across
%   subjects/hemispheres.

%% Setup the directories and variables
figDir = ['/Users/abock/Dropbox (Aguirre-Brainard Lab)/bock/' ...
    'Bock_Aguirre_manuscripts/Retinotopy_Templates/Figures/Figure3_error_in_fitting_and_prediction/raw'];
tmpDir = '~/tmpDropbox';
templateDir = '~/data';
if ~exist(tmpDir,'dir')
    mkdir(tmpDir);
end
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
templateSize = 'V2_V3';
fitType = 'V2V3';
%% Plot the correlation between the data matrix and the best prediction matrix
% Benson Retinotopic template
Tareas = load_nifti('/data/jet/abock/data/2014-10-29.areas-template.nii.gz');
Tecc = load_nifti('/data/jet/abock/data/2014-10-29.eccen-template.nii.gz');
Tpol = load_nifti('/data/jet/abock/data/2014-10-29.angle-template.nii.gz');
Tpol.vol = deg2rad(Tpol.vol) - pi/2; % convert to radians
[Tx,Ty] = pol2cart(Tpol.vol,Tecc.vol);
TV1ind = find(abs(Tareas.vol) == 1);
TV2ind = find(abs(Tareas.vol) == 2);
TV3ind = find(abs(Tareas.vol) == 3);
% Loop through subjects
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = subRuns{ss};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        % Setup for new volumes
        Ntmp = Tareas;
        Ntmp.vol = nan(size(Ntmp.vol));
        Nareas = Ntmp;
        Necc = Ntmp;
        Npol = Ntmp;
        Nsse = Ntmp;
        switch template
            case 'pRF'
                tDir = fullfile(session_dir,'pRFs','pRF_templates');
                pRF_dir = fullfile(session_dir,'pRFs','pRF_templates','decimated_templates');
                bestTemplate = 'pRF';
            case 'coarse'
                tDir = fullfile(session_dir,'pRFs',template,func,cond,fitType);
                pRF_dir = fullfile(session_dir,'pRFs','coarse_model_templates','decimated_templates');
                [~,~,sorted_templates] = find_best_template(template,tDir,hemi);
                bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
            case 'fine'
                pRF_dir = fullfile(session_dir,'pRFs','fine_model_templates',fitType,'decimated_templates');
                tDir = fullfile(session_dir,'pRFs','fine',func,cond,fitType);
                [~,~,sorted_templates] = find_best_template(template,tDir,hemi);
                bestTemplate = sorted_templates{1}(4:(strfind(sorted_templates{1},'varexp')-2));
            case 'anat'
                tDir = fullfile(session_dir,'pRFs','anat_templates');
                pRF_dir = fullfile(session_dir,'pRFs','anat_templates','decimated_templates');
                bestTemplate = 'anat';
        end
        % Create the matrices
        [obs_all,pred_all,srcind,trgind,rowV1inds,rowV2inds,rowV3inds,colV1inds,colV2inds,colV3inds] = ...
            create_template_matrices(session_dir,pRF_dir,bestTemplate,...
            hemi,func,runs,templateSize);
        % Compute the correlation
        % V1
        V1corr = nan(length(rowV1inds),1);
        for i = 1:length(rowV1inds)
            tmpobs = obs_all(rowV1inds(i),[colV2inds,colV3inds]);
            tmppred = pred_all(rowV1inds(i),[colV2inds,colV3inds]);
            tmpobs(isnan(tmpobs)) = [];
            tmppred(isnan(tmppred)) = [];
            V1corr(i) = fisher_z_corr(corr(tmpobs',tmppred'));
        end
        % V2
        V2corr = nan(length(rowV2inds),1);
        for i = 1:length(rowV2inds)
            tmpobs1 = obs_all(rowV2inds(i),colV3inds);
            tmppred1 = pred_all(rowV2inds(i),colV3inds);
            tmpobs1(isnan(tmpobs1)) = [];
            tmppred1(isnan(tmppred1)) = [];
            tmpobs2 = obs_all(rowV1inds,colV2inds(i));
            tmppred2 = pred_all(rowV1inds,colV2inds(i));
            tmpobs2(isnan(tmpobs2)) = [];
            tmppred2(isnan(tmppred2)) = [];
            %V2corr(i) = fisher_z_corr(corr(tmpobs2,tmppred2));
            V2corr(i) = fisher_z_corr(corr([tmpobs1';tmpobs2],[tmppred1';tmppred2]));
        end
        % V3
        V3corr = nan(length(colV3inds),1);
        for i = 1:length(colV3inds)
            %             tmpobs = obs_all(rowV1inds,colV3inds(i));
            %             tmppred = pred_all(rowV1inds,colV3inds(i));
            tmpobs = obs_all([rowV1inds,rowV2inds],colV3inds(i));
            tmppred = pred_all([rowV1inds,rowV2inds],colV3inds(i));
            tmpobs(isnan(tmpobs)) = [];
            tmppred(isnan(tmppred)) = [];
            V3corr(i) = fisher_z_corr(corr(tmpobs,tmppred));
        end
        % Save the surface files
        outFile = fullfile(tmpDir,[hemi '.' template '.' subject_name '.Zcorr.0.1.nii.gz']);
        tmp = load_nifti(fullfile(pRF_dir,[hemi '.areas.' bestTemplate '.nii.gz']));
        tmp.vol = nan(size(tmp.vol));
        tmp.vol([trgind{1};trgind{2};trgind{3}]) = [V1corr;V2corr;V3corr];
        save_nifti(tmp,outFile);
        % Save on fsaverage surface (using visual space)
        Sareas = load_nifti(fullfile(pRF_dir,[hemi '.areas.' bestTemplate '.nii.gz']));
        Secc = load_nifti(fullfile(pRF_dir,[hemi '.ecc.' bestTemplate '.nii.gz']));
        Spol = load_nifti(fullfile(pRF_dir,[hemi '.pol.' bestTemplate '.nii.gz']));
        NZcorr = load_nifti('/data/jet/abock/data/2014-10-29.eccen-template.nii.gz');
        NZcorr.vol = nan(size(NZcorr.vol));
        SZcorr = load_nifti(outFile);
        [Sx,Sy] = pol2cart(Spol.vol,Secc.vol);
        if strcmp(hemi,'rh')
            Sx = -Sx; % flip right hemisphere
        end
        SV1ind = find(abs(Sareas.vol) == 1);
        SV2ind = find(abs(Sareas.vol) == 2);
        SV3ind = find(abs(Sareas.vol) == 3);
        % V3
        for i = 1:length(TV3ind)
            [~,tmp] = min(sqrt( (Tx(TV3ind(i)) - Sx(SV3ind)).^2 + (Ty(TV3ind(i)) - Sy(SV3ind)).^2));
            NZcorr.vol(TV3ind(i)) = SZcorr.vol(SV3ind(tmp));
        end
        % V2
        for i = 1:length(TV2ind)
            [~,tmp] = min(sqrt( (Tx(TV2ind(i)) - Sx(SV2ind)).^2 + (Ty(TV2ind(i)) - Sy(SV2ind)).^2));
            NZcorr.vol(TV2ind(i)) = SZcorr.vol(SV2ind(tmp));
        end
        % V1
        for i = 1:length(TV1ind)
            [~,tmp] = min(sqrt( (Tx(TV1ind(i)) - Sx(SV1ind)).^2 + (Ty(TV1ind(i)) - Sy(SV1ind)).^2));
            NZcorr.vol(TV1ind(i)) = SZcorr.vol(SV1ind(tmp));
        end
        symFile = fullfile(tmpDir,[hemi '.' template '.' subject_name '.Zcorr.sym.nii.gz']);
        save_nifti(NZcorr,symFile);
    end
end
%% Save individual plots
cd(figDir);
for ss = 1:length(subjects)
    subject_name = subjects{ss};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        inFile = fullfile(tmpDir,[hemi '.' template '.' subject_name '.Zcorr.0.1.nii.gz']);
        surface_plot('rbco',inFile,subject_name,hemi,'0.1.inflated');
        savefigs('png',[hemi '.' subject_name '.Zcorr']);
        close all
    end
end
%% Average on fsaverag_sym
ct = 0;
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for ss = 1:length(subjects)
        ct = ct + 1;
        subject_name = subjects{ss};
        symFile = fullfile(tmpDir,[hemi '.' template '.' subject_name '.Zcorr.sym.nii.gz']);
        tmp = load_nifti(symFile);
        s(ct,:) = tmp.vol;
    end
end
allmean = squeeze(mean(s))';
ecc = load_nifti('/data/jet/abock/data/2014-10-29.eccen-template.nii.gz');
thresh = ecc.vol<=100;
surface_plot('rbco',allmean,'fsaverage_sym','lh','inflated',thresh);
cd(figDir);
savefigs('png','mh_mean_Zcorr');
close all
%% Save areas, ecc, and polar for Zcorr figure
cd(figDir);
ecc = load_nifti(fullfile(templateDir,'2014-10-29.eccen-template.nii.gz'));
pol = load_nifti(fullfile(templateDir,'2014-10-29.angle-template.RADS.nii.gz'));
areas = load_nifti(fullfile(templateDir,'2014-10-29.areas-template.nii.gz'));
thresh = (ecc.vol <= 6.75 & ecc.vol >= 5.75 | ecc.vol <= 12.75 & ecc.vol >= 11.75);
surface_plot('ecc',ecc.vol,'fsaverage_sym','lh','inflated',thresh);
savefigs('png','ecc_6_12');
close all
surface_plot('blueareas',areas.vol,'fsaverage_sym','lh','inflated');
savefigs('png','areas');
close all
%% Compute variance inside and outside stimulated region
in = ecc.vol<=12;
out = ecc.vol>12;
disp(mean(allmean(in&thresh)));
disp(std(allmean(in&thresh)));
disp(mean(allmean(out&thresh)));
disp(std(allmean(out&thresh)));
[H,P,CI,STATS] = ttest2(allmean(in&thresh),allmean(out&thresh));