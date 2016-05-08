%% TemplateRetinotopy MASTER script
%
%	This script will preprocess and analyze data for the TemplateRetinotopy project
%
%   MRklar - v1.0.0

%% Set session and subject name
logDir = '/data/jet/abock/LOGS';
templateDir = '/data/jet/abock/data/Retinotopy_Templates/fsaverage_sym';
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
jobNames = {'A100114K' 'A102714B' 'G101514A'};
pRFruns = {[1,2,5] [1,3,5] [1,3,5]};
movieRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
surfFunc = 's5.wdrf.tf.surf';
volFunc = 's5.wdrf.tf';
%% Run preprocessing
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    job_name = jobNames{ss};
    
    outDir = fullfile(session_dir,'preprocessing_scripts');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    numRuns = 6; % number of bold runs
    reconall = 0;
    slicetiming = 1; % correct slice timings
    refvol = 1; % motion correct to 1st TR
    filtType = 'detrend';
    lowHz = 0.01;
    highHz = 0.10;
    physio = 1;
    motion = 1;
    task = 0;
    localWM = 1;
    anat = 1;
    amem = 20;
    fmem = 50;
    create_preprocessing_scripts(session_dir,subject_name,outDir,logDir,...
        job_name,numRuns,reconall,slicetiming,refvol,filtType,lowHz,highHz,...
        physio,motion,task,localWM,anat,amem,fmem);
end
%% Run the above shell scripts in terminal
%   for example:
%   sh /data/jet/abock/data/Retinotopy_Template/AEK/10012014/preprocessing_scripts/submit_A102714B_all.sh

%% Project retinotopic templates to subject space
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    project_template(session_dir,subject_name)
end
%% Decimate the anatomical template files
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    tDir = fullfile(session_dir,'anat_templates');
    decimate_templates(subject_name,tDir);
end
%% pRF analysis
% Generates a population receptive field (pRF) estimate using data obtained
%   while subjects viewed retinotopy stimuli (e.g. drifting bars). The
%   resulting pRF maps are then averaged across runs.
% If ROI = 'occipital', the averaged maps are plotted on the fsaverage_sym
%   surface, averaged across hemispheres, and converted to a format for
%   template fitting using Mathematica.
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = pRFruns{ss};
    
    hemis = {'lh' 'rh'};
    srcROI = 'cortex';
    for runNum = runs
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            run_pRF(session_dir,subject_name,runNum,hemi,srcROI,surfFunc);
        end
    end
end
%% Average pRF
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = pRFruns{ss};
    srcROI = 'cortex';
    average_pRF(session_dir,subject_name,runs,srcROI);
end
%% Create the fsaverage_sym templates in Mathematica
% output is found in :
%    fullfile(templateDir,'pRFs','coarse_model_templates');
%
%% Copy an anatomical template, to be used later
mkdir(fullfile(templateDir,'anat_templates'));
system(['cp ~/data/2014-10-29.areas-template.nii.gz ' ...
    fullfile(templateDir,'anat_templates','lh.areas.anat.nii.gz')]);
%% Convert coarse_model_templates, decimate
%   Takes the .mgz output from Mathematica (above), convertes to nii.gz and
%   separates out the pol, ecc, and areas maps into individual volumes.
%
%   Assumes this has to be done in terminal: (IN LINUX!)
%       cd $SUBJECTS_DIR/<subject_name>/surf
%       mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%       mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated
convert_Mathematica_templates(templateDir);
%% Project templates to subject space
% This will create scripts to project the templates from the fsaverag_sym 
% surface to the individual subject surfaces
tDir = fullfile(templateDir,'pRFs','coarse_model_templates');
tFiles = listdir(fullfile(tDir,'*.nii.gz'),'files');
hemi = 'lh';
srcsubject = 'fsaverage_sym';
mem = 5;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    submit_name = jobNames{ss};
    outDir = fullfile(session_dir,'project_templates_scripts');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    sDir = fullfile(session_dir,'pRFs','coarse_model_templates');
    if ~exist(sDir,'dir')
        mkdir(sDir);
    end
    job_string = cell(1,length(tFiles));
    for tt = 1:length(tFiles);
        sval = fullfile(tDir,tFiles{tt});
        tval = fullfile(sDir,tFiles{tt});
        job_name = tFiles{tt};
        matlab_string = ['mri_surf2surf(''' srcsubject ''',''' subject_name ''',''' ...
            sval ''',''' tval ''',''' hemi ''');'];
        create_job_shell(outDir,job_name,matlab_string);
        job_string{tt} = fullfile(outDir,[job_name '.sh']);
    end
    create_submit_shell(outDir,logDir,submit_name,job_string,mem)
end
%% Run the above shell scripts


%% Decimate the templates
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    tDir = fullfile(session_dir,'pRFs','coarse_model_templates');
    decimate_templates(subject_name,tDir);
end
%% Decimate the bold runs
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    decimate_bold(session_dir,subject_name,volFunc);
end
%% Decimate the bold runs
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    decimate_bold(session_dir,subject_name,volFunc);
end
%% Create cluster shell scripts (anat)
templateType = 'anat';
tcPart = 'full';
leaveOut = '0';
V2V3 = '0'; % '0' = V1-V2, V1-V3; '1' = V1-V2, V1-V3, V2-V3
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    outDir = fullfile(session_dir,'fit_template_scripts',templateType);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    saveDir = fullfile(session_dir,'pRFs',templateType,volFunc,'Movie','V1');
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    runs = movieRuns{ss};
    create_regress_template_scripts(session_dir,templateType,outDir,runs,...
        volFunc,saveDir,tcPart,leaveOut,V2V3,logDir);
end
%%


















%% Deprecated










%% Prepare for Mathematica
prepare_pRF_Mathematica(session_dir,subject_name)

%% Run template fitting in Mathematica
% In a notebook in Mathematica, run the template fitting developed by Noah
%   Benson.
%% Create the pRF template .nii.gz files, using the output .mgz from Mathematica (above)
% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    create_pRF_template(session_dir,subject_name);
    tDir = fullfile(session_dir,'pRFs','pRF_templates');
    decimate_templates(session_dir,subject_name,tDir);
end
%% If doing correlation template fitting
%
%   see below
%
%%%
%% Create cluster shell scripts (pRF)
templateType = 'pRF';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Create cluster shell scripts (anat)
templateType = 'anat';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Convert coarse_model_templates
%Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
%   note: this can also be run on the cluster
convert_Mathematica_templates(session_dir);

%% Decimate surfaces
% This has to be done in terminal IN LINUX!
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated

%% Decimate the pRF templates and bold runs
%   This can also be run on the cluster
tDir = fullfile(session_dir,'pRFs','coarse_model_templates');
decimate_templates(session_dir,subject_name,tDir);
decimate_bold(session_dir,subject_name,func);

%% Create cluster shell scripts (coarse)
templateType = 'coarse';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Run 'coarse' template fits on the cluster ('regress_template')
% i.e. run the shell scripts created above

%% Find the best template
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
hemis = {'lh' 'rh'};
templateType = 'coarse';
func = 's5.dbrf.tf';
fitType = 'V2V3'; % 'V1 = V1<->V2, V1<->V3; 'V2V3' =  V1<->V2, V1<->V3, AND V2<->V3
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    disp(session_dir);
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        disp(hemi);
        tDir = fullfile(session_dir,'pRFs',templateType,func,'Movie',fitType);
        [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,[],[],[],fitType);
        disp(params(1));
    end
end
%% Create fine templates, centered on the best template (above)
% In a notebook in Mathematica, create the fine templates

%% Convert fine_model_templates
%Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
% note: this can also be run on the cluster
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            tDir = fullfile(session_dir,'pRFs','fine_model_templates','V2V3');
        else
            tDir = fullfile(session_dir,'pRFs','fine_model_templates','V1');
        end
        convert_Mathematica_fine_templates(session_dir,tDir);
        decimate_templates(session_dir,subject_name,tDir);
    end
end
%% Create cluster shell scripts (fine)
templateType = 'fine';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Run the 'fine' template fits on the cluster ('regress_template')
% i.e. run the shell scripts created above

%% Overwrite feat dir


