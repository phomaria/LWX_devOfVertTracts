% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
% addpath(genpath([rootDir 'proj-5a74d1e26ed91402ce400cca/']));

% diary log_LWXdevofVerticalWM_2_analyzeData_step3_regression_hv_ind_fu3
wm_measure = 'fa'; %fa, ad, rd, md, od, icvf, isovf
beh_measures = {'age_mo', 'c_lit', 'c_vm', 'c_fm'};

%% Tractography

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load([rootDir 'supportFiles/LWX_data_' wm_measure '_age_tractz.mat'])
clearvars -except wm_z wm_childrenOnly wm_childrenOnly_z list_tract wm_measure rootDir covariates sub sub_childrenOnly beh_measure beh_measures vertical

% diary on

% Get index matrices for hypothesis-driven grouping of WM tracts. LH only right now.
for k = 1:size(list_tract, 1)
    
    % Indices of horizontal tracts.
    h_idx(k) = strcmp(list_tract{k}, 'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
        || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
        || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
        || strcmp(list_tract{k}, 'leftArc') || strcmp(list_tract{k}, 'rightArc') ...
        || strcmp(list_tract{k}, 'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3');
    
    % Indices of vertical tracts.
    v_idx(k) = strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
        || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
        || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
        || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
        || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
        || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang');
    
end

% Select the measurements of the tracts that I care about and categorize them into h or v. Convert all zeros to NaN.
ht = wm_childrenOnly(:, h_idx); ht(ht==0) = NaN;
vt = wm_childrenOnly(:, v_idx); vt(vt==0) = NaN;

% tract group mean z-score: performed within-category, across-subjects to control for subject-level
% differences in 'reactivity' while keeping subject-level differences in the WM measurement of interest
for c = 1:size(ht, 2)
    ht(isnan(ht(:, c)), c)=nanmean(ht(:, c)); % set NaNs to the mean for now
end
% z_ht = zscore(ht, 0, 1);
z_ht = (nanmean(ht, 1) - ht)./nanstd(ht, [], 1);

for c = 1:size(vt, 2)
    vt(isnan(vt(:, c)), c)=nanmean(vt(:, c)); % set NaNs to the mean for now
end
% z_vt = zscore(vt, 0, 1);
z_vt = (nanmean(vt, 1) - vt)./nanstd(vt, [], 1);

% Subset list_tracts so that we can call the correct tract names later.
list_tract_ht = list_tract(h_idx, :);
list_tract_vt = list_tract(v_idx, :);

%% Behavior.

% Read in behavioral data.
load([rootDir 'supportFiles/LWX_all_groupings.mat']);

% Header for LWX Behavioral Data.
header = {'subID', 'age_mo', 'vmi', 'vp', 'mc', 'pegs_dom', 'pegs_ndom', 'lwi', 'spell', 'wa', 'sos', 'c_vm', 'c_fm', 'c_lit', 'gp_age', 'gp_lit', 'gp_vm', 'gp_fm'};

% Select only those subjects for whom we have both behavioral and tractography data.
subID_list = data_lwx(find(ismember(sub_childrenOnly, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'subID')));

% Subselect only data for children and for the covariates of interest. Remove unnecessary columns (e.g., day).
for b = 1:length(beh_measures)
    beh(:, b) = data_lwx(find(ismember(sub_childrenOnly, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, beh_measures{1, b})));
end

% Get measure-specific z-scores, based on within measure means/std.
z_beh = (beh - nanmean(beh, 1))./nanstd(beh, 1);

% Set significance thresholds.
p_val = 0.05; p_bonf = p_val/4;

fcount = 0;
for b = 1:length(beh_measures)
    
    % Concatenate each tract category into a panel and construct BIG matrix.
    if strcmp(beh_measures{b}, 'age_mo')
        z_mat = (z_beh(:, 1));
    elseif strcmp(beh_measures{b}, 'c_lit')
        z_mat = (z_beh(:, 2));
    elseif strcmp(beh_measures{b}, 'c_vm')
        z_mat = (z_beh(:, 3));
    elseif strcmp(beh_measures{b}, 'c_fm')
        z_mat = (z_beh(:, 4));
    end
    
    for v = 1:length(list_tract_vt)
        
        for h = 1:length(list_tract_ht)
            
            % Display tract name.
            disp(['----------- ' beh_measures{b} ': ' list_tract_ht{h} ' and ' list_tract_vt{v} ' -----------']);
            
            % Handle Missing Not At Random (MNAR) data -- right now replacing with mean (not best option).
            
            %% Determine if Random Effect is needed for 'subject' and Perform GLM.
            
            % BEH: Design matrix, (i.e., behavior), fixed effects
            X = blkdiag(z_mat, z_mat); 
            
            % CONSTANT: Design matrix, random effects
            Z = ones(size(X, 1), 1);
            
            % GROUPING Variable:
            G = cat(1, subID_list, subID_list);
            
            % TRACTS: Define the outcome data, (i.e., [WMh, WMv]).
            y = cat(1, z_ht(:, h), z_vt(:, v));
            
            % Perform fitting procedure with fitlm: 'Recog ~ WMa + WMm + Wmz + 1|subj';
            mdl = fitlmematrix([ones(size(X, 1), 1) X], y, Z, G, 'FixedEffectPredictors', [{'intercept'}; list_tract_ht(h); list_tract_vt(v)], ...
                'RandomEffectPredictors', {'intercept'}, 'RandomEffectGroups', {'Subject'});
            tag = 1;
            
            format short
            
            % Display betas.
            disp(mdl.Coefficients)
            
            % Display R2.
            disp(['Rsquared is: ' num2str(mdl.Rsquared.Ordinary) '.']);
            
            % Posthoc comparison for v > h.
            [pVal, F, df1, df2] = coefTest(mdl, [0 -1 1]);
            disp(['Planned comparison ' list_tract_vt{v} ' > ' list_tract_ht{h} ' for ' beh_measures{b} ' : F(' ...
                num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])
            
            % If the beta for the vertical tract is significantly *greater* than the beta for the horizontal tract, 
            coef_check = dataset2cell(mdl.Coefficients);
            if pVal <= .05 && coef_check{strcmp(coef_check(:, 1), list_tract_vt{v}), 2} >= coef_check{strcmp(coef_check(:, 1), list_tract_ht{h}), 2}
                
                % Code that location for the vertical tract.
                keep(v, h) = 1;
                
            %If the beta for the vertical tract is significantly *less* than the beta for the horizontal tract, 
            elseif pVal <= .05 && coef_check{strcmp(coef_check(:, 1), list_tract_vt{v}), 2} < coef_check{strcmp(coef_check(:, 1), list_tract_ht{h}), 2}
                
                % Code that location for the horizontal tract.
                keep(v, h) = -1;
             
                % If there is no difference,
            else
                
                % Code that location for 'none'.
                keep(v, h) = 0;
                
            end
            
            % Add annotation for passing bonferonni-corrected threshold.
            if pVal <= .05/120
                
                % Code that location for the bonf sig.
                keep_bf(v, h) = 1;
    
                % If there is no difference,
            else
                
                % Code that location for bonf ns.
                keep_bf(v, h) = 0;
                
            end
            
        end
        
    end
    
    fcount = fcount + 1;
    
    figure(fcount)
    spy(ones(size(keep)), 'k', 1)
    hold on;
    spy(keep == 1, 'r', 50) 
    spy(keep == -1, 'b', 50) 
    spy(keep_bf == 1, '*k', 8)
    spy(keep_bf == 1, 'ok', 15)

    xticks(1:length(list_tract_ht));
    xticklabels(list_tract_ht);
    xlabel('Horizontal Tracts');
    xtickangle(45)
    
    yticks(1:length(list_tract_vt));
    yticklabels(list_tract_vt);
    ylabel('Vertical Tracts');
    
    % Title.
    if strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Fine-Motor Skill.');
        
    elseif strcmp(beh_measures{b}, 'age_mo') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Age (months).')
    elseif strcmp(beh_measures{b}, 'c_lit') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Literacy.')
    elseif strcmp(beh_measures{b}, 'c_vm') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Visual-Motor Skill.')
    elseif strcmp(beh_measures{b}, 'c_fm') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Fine-Motor Skill.');
        
    end
    
    print([rootDir 'plots/plot_spy_' beh_measures{b} '_' wm_measure '_hv'], '-dpng')
    print([rootDir 'plots/eps/plot_spy_' beh_measures{b} '_' wm_measure '_hv'], '-depsc')
    
    hold off;
    
end


% diary off