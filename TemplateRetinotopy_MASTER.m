%% TemplateRetinotopy MASTER script
%
%	This script will preprocess and analyze data for the TemplateRetinotopy project
%
%   MRklar - v1.0.0

%% Set session and subject name
logDir = '/data/jet/abock/LOGS';
dataDir = '/data/jet/abock/data';
templateDir = fullfile(dataDir,'Retinotopy_Templates','fsaverage_sym');
sessions = {...
    fullfile(dataDir,'Retinotopy_Templates','AEK','10012014') ...
    fullfile(dataDir,'Retinotopy_Templates','ASB','10272014') ...
    fullfile(dataDir,'Retinotopy_Templates','GKA','10152014') ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
jobNames = {'A100114K' 'A102714B' 'G101514A'};
outNames = {'AEK' 'ASB' 'GKA'};
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
%   sh /data/jet/abock/data/Retinotopy_Template/AEK/10012014/preprocessing_scripts/submit_A100114K_all.sh

%% Temporally filter the rf.nii.gz
func = 'rf';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    boldDirs = find_bold(session_dir);
    for i = 1:length(boldDirs)
        runNum = i;
        temporal_filter(session_dir,runNum,func,filtType);
    end
end
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
hemis = {'lh' 'rh'};
srcROI = 'cortex';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = pRFruns{ss};
    for runNum = runs
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            run_pRF(session_dir,subject_name,runNum,hemi,srcROI,surfFunc);
        end
    end
end
%% Average pRF
srcROI = 'cortex';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = pRFruns{ss};
    average_pRF(session_dir,subject_name,runs,srcROI);
end
%% Prepare pRF template for fitting in Mathematica
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    outName = outNames{ss};
    prepare_pRF_Mathematica(session_dir,subject_name,outName);
end
%% Create pRF templates following Mathematica template fitting
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    create_pRF_template(session_dir,subject_name);
    tDir = fullfile(session_dir,'pRFs','pRF_templates');
    decimate_templates(session_dir,subject_name,tDir);
end
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
hemis = {'lh' 'rh'};
srcsubject = 'fsaverage_sym';
mem = 5;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    submit_name = jobNames{ss};
    outDir = fullfile(session_dir,'project_templates_scripts');
    system(['rm -rf ' outDir]);
    mkdir(outDir);
    sDir = fullfile(session_dir,'pRFs','coarse_model_templates');
    system(['rm -rf ' sDir]);
    mkdir(sDir);
    job_string = cell(1,length(tFiles)*2);
    ct = 0;
    for tt = 1:length(tFiles);
        sval = fullfile(tDir,tFiles{tt});
        for hh = 1:length(hemis)
            ct = ct +1;
            hemi = hemis{hh};
            job_name = [hemi tFiles{tt}(3:end)];
            tval = fullfile(sDir,job_name);
            if strcmp(hemi,'lh')
            matlab_string = ['mri_surf2surf(''' srcsubject ''',''' subject_name ''',''' ...
                sval ''',''' tval ''',''lh'');'];
            else
                matlab_string = ['mri_surf2surf(''' srcsubject ''',''' subject_name '/xhemi'',''' ...
                sval ''',''' tval ''',''lh'');']; 
            end
            create_job_shell(outDir,job_name,matlab_string);
            job_string{ct} = fullfile(outDir,[job_name '.sh']);
        end
    end
    create_submit_shell(outDir,logDir,submit_name,job_string,mem)
end
%% Run the above shell scripts
%   for example:
%   sh /data/jet/abock/data/Retinotopy_Templates/AEK/10012014/project_templates_scripts/A100114K.sh

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
%% Create cluster shell scripts (pRF,anat,coarse)
tTypes = {'anat' 'pRF' 'coarse'};
for tt = 1:length(tTypes)
    templateType = tTypes{tt};
    tcPart = 'full';
    leaveOut = '0';
    V2V3 = '1'; % '0' = V1-V2, V1-V3; '1' = V1-V2, V1-V3, V2-V3
    if strcmp(V2V3,'1')
        fitType = 'V2V3';
    else
        fitType = 'V1';
    end
    for ss = 1:length(sessions)
        session_dir = sessions{ss};
        outDir = fullfile(session_dir,'fit_template_scripts',templateType);
        system(['rm -rf ' outDir]);
        mkdir(outDir);
        saveDir = fullfile(session_dir,'pRFs',templateType,volFunc,'Movie',fitType);
        system(['rm -rf ' saveDir]);
        mkdir(saveDir);
        runs = movieRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,...
            volFunc,saveDir,tcPart,leaveOut,V2V3,logDir);
    end
end
%% Submit the regress scripts (must be run from chead on UPenn cluster)
script_dirs = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014/fit_template_scripts/' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014/fit_template_scripts/' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014/fit_template_scripts/' ...
    };
logDir = '/data/jet/abock/LOGS';
hemis = {'lh' 'rh'};
for ss = 1:length(script_dirs)
    script_dir = script_dirs{ss};
    cDirs = listdir(script_dir,'dirs');
    for i = 1:2%length(cDirs);
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            system(['rm ' fullfile(logDir,'*')]);
            pause(5);
             scriptDir = fullfile(script_dir,cDirs{i});
             cd(scriptDir);
             inScript=['submit_' hemi '_regress.sh'];
            system(['sh ' fullfile(scriptDir,inScript)]);
            disp(['running ' fullfile(scriptDir,inScript)]);
            pause(5);
            system(['qstat > ' fullfile(logDir,'tmp.txt')]);
            pause(5);
            fid = fopen(fullfile(logDir,'tmp.txt'));
            tmp = fread(fid);
            fclose(fid);
            while ~isempty(tmp)
                system(['qstat > ' fullfile(logDir,'tmp.txt')]);
                pause(5);
                fid = fopen(fullfile(logDir,'tmp.txt'));
                tmp = fread(fid);
                fclose(fid);
            end
            fclose('all');
            rerun_regress(inScript,scriptDir,logDir);
            pause(5);
            system(['qstat > ' fullfile(logDir,'tmp.txt')]);
            pause(5);
            fid = fopen(fullfile(logDir,'tmp.txt'));
            tmp = fread(fid);
            fclose(fid);
            while ~isempty(tmp)
                system(['qstat > ' fullfile(logDir,'tmp.txt')]);
                pause(5);
                fid = fopen(fullfile(logDir,'tmp.txt'));
                tmp = fread(fid);
                fclose(fid);
            end
            fclose('all');
        end
    end
end
%%
hemis = {'lh' 'rh'};
templateType = 'coarse';
func = 's5.wdrf.tf';
fitType = 'V2V3'; % 'V1 = V1<->V2, V1<->V3; 'V2V3' =  V1<->V2, V1<->V3, AND V2<->V3
for ss = 1%:length(sessions)
    session_dir = sessions{ss};
    disp(session_dir);
    for hh = 1%:length(hemis)
        hemi = hemis{hh};
        disp(hemi);
        tDir = fullfile(session_dir,'pRFs',templateType,func,'Movie',fitType);
        if strcmp(templateType,'anat')
            vals = load(fullfile(tDir,[hemi '.' templateType '.varexp.txt']));
            if strcmp(fitType,'V1')
                varexp = nansum(vals(1:2));
            elseif strcmp(fitType,'V2V3')
                varexp = nansum(vals(1:3));
            end
            disp(varexp);
        else
            [varexp,params,sorted_templates] = find_best_template(templateType,tDir,hemi,[],[],[],fitType);
            disp(params(1));
            disp(sorted_templates(1));
            disp(varexp(1));
        end
    end
end
%%

%%% finalized code up to here %%%

%%

%%
session_dir = '/Volumes/ExFat5TB/data/jet/abock/data/Retinotopy_Templates/AEK/10012014';
rawsession_dir = '/Volumes/ExFat5TB/data/jet/abock/data/Template_Retinotopy/AEK/10012014';
template = 'anat';
func = 's5.wdrf.tf';
fitType = 'V1';
cond = 'Movie';
templates = {'pRF' 'coarse' 'anat' 'raw'};
hemi = 'lh';
hh = 1;
ss = 1;
for tt = 1:length(templates)
    template = templates{tt};
    vdir = fullfile(session_dir,'pRFs',template,func,cond,fitType);
    if strcmp(template,'coarse')
        tdir = fullfile(session_dir,'pRFs',[template '_model_templates']);
        [~,~,sorted_templates] = find_best_template(template,vdir,hemi,[],[],[],fitType);
        dotinds = strfind(sorted_templates{1},'.');
        % load the 'best' pol template
        pol = load_nifti(fullfile(tdir,...
            [hemi '.pol.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']));
        % load the 'best' ecc template
        ecc = load_nifti(fullfile(tdir,...
            [hemi '.ecc.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']));
        % load the 'best' ecc template
        areas = load_nifti(fullfile(tdir,...
            [hemi '.areas.' sorted_templates{1}(dotinds(1)+1:dotinds(4)-1) '.nii.gz']));
        % pull out vals
        allVals{tt,ss,hh,1} = squeeze(pol.vol(:));
        allVals{tt,ss,hh,2} = squeeze(ecc.vol(:));
        allVals{tt,ss,hh,3} = squeeze(areas.vol(:));
    elseif strcmp(template,'pRF') || strcmp(template,'anat')
        tdir = fullfile(rawsession_dir,'pRFs',[template '_templates']);
        % load the 'best' pol template
        pol = load_nifti(fullfile(tdir,[hemi '.pol.' template '.nii.gz']));
        % load the 'best' ecc template
        ecc = load_nifti(fullfile(tdir,[hemi '.ecc.' template '.nii.gz']));
        % load the 'best' ecc template
        areas = load_nifti(fullfile(tdir,[hemi '.areas.' template '.nii.gz']));
        % pull out vals
        allVals{tt,ss,hh,1} = squeeze(pol.vol(:));
        allVals{tt,ss,hh,2} = squeeze(ecc.vol(:));
        allVals{tt,ss,hh,3} = squeeze(areas.vol(:));
    elseif strcmp(template,'raw')
        tdir = fullfile(session_dir,'pRFs');
        % load the 'best' pol template
        pol = load_nifti(fullfile(tdir,[hemi '.cortex.avg.copol.prfs.nii.gz']));
        % load the 'best' ecc template
        ecc = load_nifti(fullfile(tdir,[hemi '.cortex.avg.coecc.prfs.nii.gz']));
        % load the 'best' co template
        co = load_nifti(fullfile(tdir,[hemi '.cortex.avg.co.prfs.nii.gz']));
        % pull out vals
        allVals{tt,ss,hh,1} = squeeze(pol.vol(:));
        allVals{tt,ss,hh,2} = squeeze(ecc.vol(:));
        allVals{tt,ss,hh,3} = squeeze(co.vol(:));
    end
end

goodind = ones(size(allVals(1,1,hh,1)));
% loop to find nans
for tt = 1:length(templates)-1
    goodind(abs(allVals{tt,ss,hh,3})>3) = 0; % outside V1-V3
    goodind(isnan(allVals{tt,ss,hh,1})) = 0; % nan value
    goodind(allVals{tt,ss,hh,2}<1.25) = 0; % ecc thresh
    goodind(allVals{tt,ss,hh,2}>8.75) = 0; % ecc thresh
end
goodind = logical(goodind);
ThreshInds{ss,hh} = zeros(size(goodind));
%ThreshInds{ss,hh}(goodind & allVals{end,ss,hh,3}>0.2236) = 1;

ThreshInds{ss,hh}(goodind) = 1;

for tt = 1:length(templates)
    % raw
    % Loop through templates
    tmppol = abs(circ_dist(allVals{4,ss,hh,1},allVals{tt,ss,hh,1}));
    tmpecc = abs(diff([allVals{4,ss,hh,2},allVals{tt,ss,hh,2}]'));
    tmppol(~ThreshInds{ss,hh}) = nan;
    tmpecc(~ThreshInds{ss,hh}) = nan;
    PolDiff(ss,hh,tt) = nanmedian(tmppol);
    EccDiff(ss,hh,tt) = nanmedian(tmpecc);
end
for tt = 1:length(templates)
    tmppol = PolDiff(:,:,tt);
    tmpecc = EccDiff(:,:,tt);
    PolDiffmean(tt) = rad2deg(circ_mean(tmppol(:)));
    PolDiffSEM(tt) = rad2deg(circ_std(tmppol(:)) / 1);
    EccDiffmean(tt) = mean(tmpecc(:));
    EccDiffSEM(tt) = std(tmpecc(:)) / 1;
end
templates
PolDiffmean



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


