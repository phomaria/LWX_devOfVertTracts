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
beh_measure = {'age_mo', 'c_lit', 'c_vm', 'c_fm'};

%% Tractography

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure '_age_raw.mat']))

% Convert into array and header for ease.
data_all_in = table2array(data_tbl);
data_all_in_header = data_tbl.Properties.VariableNames;

% Get grouping variable. NOTE: need to add lit, vm, and fm.
if strcmp(beh_measure{1}, 'age_mo')
    
    group = data_tbl.group_age;
    
end

% Get index matrices for hypothesis-driven grouping of WM tracts.
for k = 1:length(data_all_in_header)
    
    % Indices of horizontal tracts.
    h_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
        || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF') ...
        || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
        || strcmp(data_all_in_header{k}, 'leftArc') || strcmp(data_all_in_header{k}, 'rightArc') ...
        || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3');
    
    % Indices of vertical tracts.
    v_idx(k) = strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
        || strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
        || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
        || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
        || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
        || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang');
    
end

% SELECT the measurements of the tracts that I care (h and v) and
% the subjects that I care about (non-adults).
% Categorize into h or v. Convert all zeros to NaN.
ht = data_all_in(group ~= 3, h_idx); ht(ht==0) = NaN;
vt = data_all_in(group ~= 3, v_idx); vt(vt==0) = NaN;

% tract group mean z-score: performed within-category, across-subjects to control for subject-level
% differences in 'reactivity' while keeping subject-level differences in the WM measurement of interest
% for c = 1:size(ht, 2)
%     ht(isnan(ht(:, c)), c)=nanmean(ht(:, c)); % set NaNs to the mean for now
% end
% z_ht = zscore(ht, 0, 1);
z_ht = (nanmean(ht, 1) - ht)./nanstd(ht, [], 1);

% for c = 1:size(vt, 2)
%     vt(isnan(vt(:, c)), c)=nanmean(vt(:, c)); % set NaNs to the mean for now
% end
% z_vt = zscore(vt, 0, 1);
z_vt = (nanmean(vt, 1) - vt)./nanstd(vt, [], 1);

% Subset list_tracts so that we can call the correct tract names later.
list_tract_ht = data_all_in_header(h_idx);
list_tract_vt = data_all_in_header(v_idx);

%% Behavior.

% Subselect only data for children and for the covariates of interest.
beh = cat(2, data_tbl.age(group~=3), data_tbl.c_lit(group~=3), data_tbl.c_vm(group~=3), data_tbl.c_fm(group~=3));

% Get measure-specific z-scores, based on within measure means/std.
z_beh = (beh - nanmean(beh, 1))./nanstd(beh, 1);

% Set significance thresholds.
p_val = 0.05; p_bonf = p_val/4;

fcount = 0;
for b = 1:length(beh_measure)
    
    % Concatenate each tract category into a panel and construct BIG matrix.
    if strcmp(beh_measure{b}, 'age_mo')
        z_mat = (z_beh(:, 1));
    elseif strcmp(beh_measure{b}, 'c_lit')
        z_mat = (z_beh(:, 2));
    elseif strcmp(beh_measure{b}, 'c_vm')
        z_mat = (z_beh(:, 3));
    elseif strcmp(beh_measure{b}, 'c_fm')
        z_mat = (z_beh(:, 4));
    end
    
    for v = 1:length(list_tract_vt)
        
        for h = 1:length(list_tract_ht)
            
            % Display tract name.
            disp(['----------- ' beh_measure{b} ': ' list_tract_ht{h} ' and ' list_tract_vt{v} ' -----------']);
            
            % Handle Missing Not At Random (MNAR) data -- right now replacing with mean (not best option).
            
            %% Determine if Random Effect is needed for 'subject' and Perform GLM.
            
            % BEH: Design matrix, (i.e., behavior), fixed effects
            X = blkdiag(z_mat, z_mat); 
            
            % CONSTANT: Design matrix, random effects
            Z = ones(size(X, 1), 1);
            
            % GROUPING Variable:
            G = cat(1, data_tbl.subID(group ~= 3), data_tbl.subID(group ~= 3));
            
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
            disp(['Planned comparison ' list_tract_vt{v} ' > ' list_tract_ht{h} ' for ' beh_measure{b} ' : F(' ...
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
    if strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'fa')
        title('Fractional Anisotropy predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'ad')
        title('Axial Diffusion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'rd')
        title('Radial Diffusion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'md')
        title('Mean Diffusivity predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'od')
        title('Orientation Dispersion predicts Fine-Motor Skill.')
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'icvf')
        title('Neurite Density (ICVF) predicts Fine-Motor Skill.');
        
    elseif strcmp(beh_measure{b}, 'age_mo') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Age (months).')
    elseif strcmp(beh_measure{b}, 'c_lit') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Literacy.')
    elseif strcmp(beh_measure{b}, 'c_vm') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Visual-Motor Skill.')
    elseif strcmp(beh_measure{b}, 'c_fm') && strcmp(wm_measure, 'isovf')
        title('Isotropic Volume Fraction predicts Fine-Motor Skill.');
        
    end
    
    print([rootDir 'plots/plot_spy_' beh_measure{b} '_' wm_measure '_hv'], '-dpng')
    print([rootDir 'plots/eps/plot_spy_' beh_measure{b} '_' wm_measure '_hv'], '-depsc')
    
    hold off;
    
end


% diary off