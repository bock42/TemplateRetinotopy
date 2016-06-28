%% Retinotopic maps
%   Have example pRF data on the cortical surface (pol and ecc) for an
%       example subject / hemisphere
%   Show the Noah template, deformed, on the cortical surface for V1-V3
%       for the example subject (for pol and ecc)
%   Make these images for all subjects / hemisphere as they will be needed
%       for the supplementary data figures.

%% pRF data on cortical surface
%  Load in data
sessions = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
threshs = {...
    [0.648624 0.619279] ...
    [0.916082 0.871358] ...
    [0.903127 0.863659] ...
    }; % 0.9 quantile for pRF correlation
subjectNames = {'AEK' 'ASB' 'GKA'};
hemis = {'lh' 'rh'};
cd('/Users/abock/Dropbox (Aguirre-Brainard Lab)/Bock/Bock_Aguirre_manuscripts/Retinotopy_Templates/Figures/Figure4_retinotopic_maps/raw');
%%
for i = 1:length(sessions)
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        session_dir = sessions{i};
        subject_name = subjects{i};
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.coecc.avg.prfs.nii.gz']));
        ecc = tmp.vol;
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.co.avg.prfs.nii.gz']));
        co = tmp.vol;
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.copol.avg.prfs.nii.gz']));
        pol = tmp.vol;
        % flip rh to be like lh
        if strcmp(hemi,'rh')
            upper = pol>0;
            lower = pol<0;
            pol(upper) = -(pol(upper) - pi);
            pol(lower) = -(pol(lower) + pi);
        end
        % Convert polar
        pol(pol<-pi/2) = -pi/2;
        pol(pol>pi/2) = pi/2;
        % Set R-squared threshold
        thresh = co>threshs{i}(hh);
        % Plot images
        surface_plot('ecc',ecc,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_ecc_pRF'];
        savefigs('png',savename);
        close all;
        surface_plot('pol',pol,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_pol_pRF'];
        savefigs('png',savename);
        close all;
    end
end
%% Low threshold GKA lh
GKAlhthresh = 0.410236; % 0.7 quantile for GKA lh
for i = 3
    for hh = 1
        hemi = hemis{hh};
        session_dir = sessions{i};
        subject_name = subjects{i};
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.coecc.avg.prfs.nii.gz']));
        ecc = tmp.vol;
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.co.avg.prfs.nii.gz']));
        co = tmp.vol;
        tmp = load_nifti(fullfile(session_dir,'pRFs',...
            [hemi '.s5.wdrf.tf.cortex.copol.avg.prfs.nii.gz']));
        pol = tmp.vol;
        % flip rh to be like lh
        if strcmp(hemi,'rh')
            upper = pol>0;
            lower = pol<0;
            pol(upper) = -(pol(upper) - pi);
            pol(lower) = -(pol(lower) + pi);
        end
        % Convert polar
        pol(pol<-pi/2) = -pi/2;
        pol(pol>pi/2) = pi/2;
        % Set R-squared threshold
        thresh = co>GKAlhthresh;
        % Plot images
        surface_plot('ecc',ecc,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_ecc_pRF_0.7thresh'];
        savefigs('png',savename);
        close all;
        surface_plot('pol',pol,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_pol_pRF_0.7thresh'];
        savefigs('png',savename);
        close all;
    end
end
%% Deformed template on cortical surface
%  Load in data
subjectNames = {'AEK' 'ASB' 'GKA'};
hemis = {'lh' 'rh'};
cd('/Users/abock/Dropbox (Aguirre-Brainard Lab)/Bock/Bock_Aguirre_manuscripts/Retinotopy_Templates/Figures/Figure4_retinotopic_maps/raw');
for i = 1:length(sessions)
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        session_dir = sessions{i};
        subject_name = subjects{i};
        outDir = fullfile(session_dir,'pRFs','pRF_templates');
        tmp = load_nifti(fullfile(outDir,[hemi '.ecc.pRF.nii.gz']));
        ecc = tmp.vol;
        tmp = load_nifti(fullfile(outDir,[hemi '.pol.pRF.nii.gz']));
        pol = tmp.vol;
        % flip rh to be like lh
        if strcmp(hemi,'rh')
            upper = pol>0;
            lower = pol<0;
            pol(upper) = -(pol(upper) - pi);
            pol(lower) = -(pol(lower) + pi);
        end
        tmp = load_nifti(fullfile(outDir,[hemi '.areas.pRF.nii.gz']));
        areas = tmp.vol;
        % Set areas threshold (V1-V3)
        thresh = abs(areas)<=3;
        surface_plot('ecc',ecc,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_ecc_deformed'];
        savefigs('png',savename);
        close all;
        surface_plot('pol',pol,subject_name,hemi,'inflated',thresh,0.85,'hemi');
        savename = [subjectNames{i} '_surface_' hemi '_pol_deformed'];
        savefigs('png',savename);
        close all;
    end
end