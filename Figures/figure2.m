%% Schira model search space
% movie cross-correlation matrix from example subject / hemisphere. Can
%   use the set of vertices that are identified by the "gold standard" pRF,
%   template deformed V1-V3 model
% Shira models from the corners and center of the coarse parameter search
%       space (6,6,6), (7,7,7), (8,8,8)
% Shira models from the corners and center of the fine parameter search
%       space
% For each model, save an image of the cortical surface with the visual
%       areas indicated in shades of blue for V1, V2, and V3
% Need figure of the fit of the predictive model of the cross-correlation
%       matrix to each of the Shira models listed above.
% For Panel C, need to integrate the coarse and fine scale variance
%       explained plots into single plots. Have a set of 3, 2D plots
%       corresponding to each search parameter pairing.

%% Setup savedir
session_dir = '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014';
savedir = ['/Users/abock/Dropbox (Aguirre-Brainard Lab)/Bock/Bock_Aguirre_manuscripts/' ...
    'Retinotopy_Templates/Figures/Figure2_schira_model_search_space/raw'];
cd(savedir);
%% Example movie cross-correlation matrix
hemi = 'lh';
func = 's5.wdrf.tf';
cond = 'Movie';
runs = [2 4 6];
templateSize = 'full';%'V2_V3'; % 'full';
V2V3 = 1; %[0 1];
for j = V2V3
    save_template_correlation_matrices(session_dir,hemi,func,cond,runs,templateSize,j)
end
%% Manually rename above images to 'best' and 'worst'
% for GKA, '4.6.4' = 'best' and '5.2.7' = 'worst'

%% Schira models - coarse
corners = {...
    '1.1.1' '1.1.7' '1.7.1' '1.7.7' ...
    '7.1.1' '7.1.7' '7.7.1' '7.7.7' ...
    };
templateType = 'coarse';
func = 's5.wdrf.tf';
cond = 'Movie';
hemi = 'lh';
session_dir = '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014';
tdir = fullfile(session_dir,'pRFs',templateType,func,cond,'V2V3');
subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% Get the variance explained by each template
[varexp,~,sorted_templates] = find_best_template(templateType,tdir,hemi);
% Find the variance explained by the 'corners'
cVarexp = nan(length(corners),1);
for cc = 1:length(corners)
    for tt = 1:length(sorted_templates)
        templateName = sorted_templates{tt}(4:(strfind(sorted_templates{tt},'varexp')-2));
        switch corners{cc}
            case templateName
                cVarexp(cc) = varexp(tt);
        end
    end
end
%% View/Save the templates
templatedir = '/data/jet/abock/data/Retinotopy_Templates/fsaverage_sym/pRFs/coarse_model_templates';
show_curv = 0; % don't show curvature (i.e. gray brain)
cd(savedir);
for j = 1:length(corners)
    tmp = load_nifti(fullfile(templatedir,['lh.areas.' corners{j} '.nii.gz']));
    areas = tmp.vol;
    thresh = abs(areas)<=3;
    surface_plot('blueareas',areas,'fsaverage_sym','lh','inflated',thresh,0.85,'hemi',show_curv);
    title(['Variance Explained = ' num2str(cVarexp(j))],'FontSize',20);
    savename = ['fsaverage_sym_' templateType '_' corners{j} '_varexp_' num2str(cVarexp(j))];
    savefigs('png',savename);
    close all;
end
%% Schira models - fine
% corners = {...
%     '1.1.1' '1.1.5' '1.5.1' '1.5.5' ...
%     '5.1.1' '5.1.5' '5.5.1' '5.5.5' ...
%     };
% templateType = 'fine';
% func = 's5.dbrf.tf';
% cond = 'Movie';
% hemi = 'lh';
% session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014';
% tdir = fullfile(session_dir,'pRFs',templateType,func,cond,'V1');
% subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% templatedir = fullfile(session_dir,'pRFs','fine_model_templates','V1');
% % Get the variance explained by each template
% [varexp,~,sorted_templates] = find_best_template(templateType,tdir,hemi);
% % Find the variance explained by the 'corners'
% cVarexp = nan(length(corners),1);
% for cc = 1:length(corners)
%     for tt = 1:length(sorted_templates)
%         templateName = sorted_templates{tt}(4:(strfind(sorted_templates{tt},'varexp')-2));
%         switch corners{cc}
%             case templateName
%                 cVarexp(cc) = varexp(tt);
%         end
%     end
% end
%% View/Save the templates
% show_curv = 0; % don't show curvature (i.e. gray brain)
% cd(savedir);
% for j = 1:length(corners)
%     tmp = load_nifti(fullfile(templatedir,['lh.areas.' corners{j} '.nii.gz']));
%     areas = tmp.vol;
%     thresh = abs(areas)<=3;
%     surface_plot('blueareas',areas,subject_name,'lh','inflated',thresh,0.85,'hemi',show_curv);
%     title(['Variance Explained = ' num2str(cVarexp(j))],'FontSize',20);
%     savename = ['GKA_' templateType '_' corners{j} '_varexp_' num2str(cVarexp(j))];
%     savefigs('png',savename);
%     close all;
% end
%% Variance explained plots
% Make the color scale 'parula' (good for colorblind, b/w, etc), square
% pixels
session_dirs = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014' ...
    };
templates = {'coarse'};
func = 's5.wdrf.tf';
cond = 'Movie';
plot_error_bars = 0;
figsave = 1; % save figures
cd(savedir);
V2V3 = 1;%[0 1];
for j = V2V3
    for jj = 1:length(templates)
        template = templates{jj};
        plot_template_mesh(session_dirs,template,func,cond,j);
    end
end
%% Calculate distance from best template
cd(savedir);
session_dirs = {...
    '/data/jet/abock/data/Retinotopy_Templates/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy_Templates/ASB/10272014' ...
    '/data/jet/abock/data/Retinotopy_Templates/GKA/10152014' ...
    };
hemis = {'lh' 'rh'};
func = 's5.wdrf.tf';
cond = 'Movie';
templateScale = 30; % for fine scale (0.1 / 30)
vals = -0.8:1/templateScale:0.8;
[sumX,sumY,sumZ] = meshgrid(vals,vals,vals);
sumX = sumX(:);
sumY = sumY(:);
sumZ = sumZ(:);
numComps = length(session_dirs)*length(hemis);
% Calculate distance
ct = 0;
for j = 1:length(vals)
    for k = 1:length(vals)
        for k = 1:length(vals)
            ct = ct + 1;
            % Minkowski sum (city-block geometry)
            dists(ct) = abs(sumX(ct)) + abs(sumY(ct)) + abs(sumZ(ct));
            % Euclidean distance
            %dists(ct) = sqrt(sumX(ct)^2 + sumY(ct)^2 + sumZ(ct)^2);
        end
    end
end
dists = dists * 10; % 0.1 = 1 step
%% Plot error by steps from best template
V2V3 = 1;%[0 1];
for i = V2V3
    clear varAll meanVals stdVals cvarexp cx cy cz fvarexp params fx fy fz
    if i
        saveName = 'V2V3';
    else
        saveName = 'V1';
    end
    ct = 0;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for ss = 1:length(session_dirs)
            ct = ct + 1;
            session_dir = session_dirs{ss};
            % coarse
            ctdir = fullfile(session_dir,'pRFs','coarse',func,cond,saveName);
            [cvarexp,~,~,cx,cy,cz] = plot_template_fits(ctdir,'coarse',hemi,saveName);
            varAll(ct,:) = cvarexp;
            X(ct,:) = cx;
            Y(ct,:) = cy;
            Z(ct,:) = cz;
%            % fine
%             ftdir = fullfile(session_dir,'pRFs','fine',func,cond,saveName);
%             [fvarexp,params,~,fx,fy,fz] = plot_template_fits(ftdir,'fine',hemi,saveName);
%             fx = fx + params(1).FCx; % adjust fine center
%             fy = fy + params(1).FCy; % adjust fine center
%             fz = fz + params(1).psi; % adjust fine center
%            % combine 'coarse' and 'fine'
%             varAll(ct,:) = [cvarexp',fvarexp'];
%             X(ct,:) = [cx,fx];
%             Y(ct,:) = [cy,fy];
%             Z(ct,:) = [cz,fz];
        end
    end
    % Fill in subject values
    sumAll = nan(numComps,length(vals),length(vals),length(vals));
    sumCt = zeros(length(vals),length(vals),length(vals));
    for j = 1:length(hemis)*length(session_dirs)
        for k = 1:length(X)
            tmpCoords = [X(j,k),Y(j,k),Z(j,k)];
            tmpCoords = round(tmpCoords*templateScale) + round(length(vals)/2);
            tmpVarexp = varAll(j,k);
            sumAll(j,tmpCoords(1),tmpCoords(2),tmpCoords(3)) = tmpVarexp;
            sumCt(tmpCoords(1),tmpCoords(2),tmpCoords(3)) = ...
                sumCt(tmpCoords(1),tmpCoords(2),tmpCoords(3)) + 1;
        end
    end
    allVar = reshape(sumAll,size(sumAll,1),size(sumAll,2)*size(sumAll,3)*size(sumAll,4));
    % Get max dists
    allMean = nanmean(allVar);
    distind = find(~isnan(allMean));
    maxDist = max(dists(distind))+0.001;
    % Bin using linspace
    bins = 0:1:maxDist;
    %bins = linspace(min(dists),floor(maxDist),floor(maxDist)+1);
    [N,bin] = histc(dists,bins);
    for b = 1:length(N)
        tmp = b == bin;
        binVals = allVar(:,tmp);
        binVals = binVals(:);
        meanVals(b) = nanmean(binVals);
        stdVals(b) = nanstd(binVals);
    end
    % Plot Distance vs Variance
    badind = isnan(meanVals);
    meanVals(badind) = [];
    stdVals(badind) = [];
    bins(badind) = [];
    xx = 0:max(bins)/10000:max(bins);
    fm = fit(bins',meanVals','smoothingspline');
    fs = fit(bins',stdVals','smoothingspline');
    mm = fm(xx)';
    ss = fs(xx)';
    % Plot
    figure('units','normalized','position',[0 0 1 1]);
    h1 = fill([(xx),fliplr(xx)],[(mm+ss),fliplr(mm-ss)],...
        [0.75 0.75 0.75],'LineWidth',1,'EdgeColor','none'); hold on;
    h2 = plot(xx,mm,'k','LineWidth',2);
    h3 = plot(bins,meanVals,'.','MarkerSize',20);
    legend([h2 h1 h3],{'mean','std','data'});
    %xlim([0 max(bins)]);
    xlim([0 15])
    ylim([0 0.25]);
    %set(gca,'XTick',bins);
    set(gca,'XTick',0:1:15);
    set(gca,'YTick',[0 0.05 0.1 0.15 0.2 0.25]);
    axis square;
    xlabel('Steps from best template','FontSize',20);
    ylabel('Variance explained','FontSize',20);
    set(gca,'FontSize',15);
    savefigs('pdf',['Distance_vs_variance_explained_' saveName]);
    close all;
end