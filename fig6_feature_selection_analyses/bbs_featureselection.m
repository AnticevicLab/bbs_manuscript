%matlab -nodisplay
% -- 1) ACROSS-SUBJECTS MODEL PERFORMANCE

% -- Read the data
pc1 = dlmread('../fig3_univariate_mapping_analyses/BSNIP_n436_dPC1.csv'); % matrix of subj by pc score
pc1map = dlmread('../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC1_ztstat.pscalar.txt');

pc2 = dlmread('../fig3_univariate_mapping_analyses/BSNIP_n436_dPC2.csv');
pc2map = dlmread('../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC2_ztstat.pscalar.txt');

pc3 = dlmread('../fig3_univariate_mapping_analyses/BSNIP_n436_dPC3.csv'); 
pc3map = dlmread('../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC3_ztstat.pscalar.txt');

pc4 = dlmread('../fig3_univariate_mapping_analyses/BSNIP_n436_dPC4.csv');
pc4map = dlmread('../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC4_ztstat.pscalar.txt');

pc5 = dlmread('../fig3_univariate_mapping_analyses/BSNIP_n436_dPC5.csv');
pc5map = dlmread('../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC5_ztstat.pscalar.txt');

% -- Read behavioral data (numeric data only)
beh = dlmread('../data_symptoms/BSNIP_PSD_N436_Behavior.tsv');
pcloadings =  dlmread('../fig1_symptom_pca_analyses/BSNIP_N436_BACSPANSS_PCArotations.tsv');

% -- Read neural data
gbcraw=dlmread('../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz.txt')';

% -- Calculate ΔGBC (demean by group mean for each parcel)
deltagbc = gbcraw - repmat(mean(gbcraw,1), 436,1); % mean for every parcel
%dlmwrite('BSNIP_PSD_N436_deltaGBC.txt', deltagbc, 'delimiter','\t') % save the ΔGBC values
gbcdata = deltagbc;

% -- Run the selection process for just PC3 here, as example
betaPC_abs = abs(pc3map); % -- Compute absolute values of the β_PC map
pc=pc3;
betaPCdat = pc3map;
pc_n=3;

nsub = length(pc);
sub_vector = 1:nsub;
betaPC_abs_sorted = sort(betaPC_abs);

% -- Compute projected PC score for holdout subject (repeat for all subjects)
pcproj = zeros(1,436);
for k = 1:nsub; %select subject to leave out
	% project behavior
	behminsub = normalize(beh(1:end ~= k,:));
	pcacoeff = pca(behminsub, 'Centered',false);
	
	%pcaall = pca(normalize(beh), 'Centered',false);
	
	pcacoeffcorrected = pcacoeff;
	for pcno = 1:36
		if corr(pcacoeff(:,pcno), pcloadings(:,pcno)) < 0
			pcacoeffcorrected(:,pcno) = -1 * pcacoeff(:,pcno);
		else
			pcacoeffcorrected(:,pcno) = 1 * pcacoeff(:,pcno);
		end
	end
	
	behminsubm  = mean(beh(1:end ~= k,:));
	behminsubsd = std(beh(1:end ~= k,:));
	behsubnorm = (beh(k,:)-behminsubm)./behminsubsd;
	pcprojsub = behsubnorm * pcacoeffcorrected; % this is the mtx of projected PC scores for the holdout subject
	pcproj(1,k) = pcprojsub(pc_n);
end

% -- Compute β_PC GBC maps for N-i subjects, iterate across all N subjects
betaPC_holdoutsub_full = zeros(718,436);
% Select subject to leave out
for k = 1:nsub;

	% Compute β_PC GBC map without the selected subject
	pcminsub = pc(1:end ~= k);
	gbc_holdoutsub_full = gbcdata(1:end ~= k,:);
	
	betaPC_holdoutsub_est = zeros(1,718);
	betaPC_holdoutsub_tStat = zeros(1,718);

	% Compute the β for each parcel
	for parc = 1:718;
		mdl = fitlm(pcminsub,gbc_holdoutsub_full(:,parc)); 
		betaPC_holdoutsub_est(1,parc) = mdl.Coefficients.Estimate(2);
		betaPC_holdoutsub_tStat(1,parc) = mdl.Coefficients.tStat(2);
	end
	betaPC_holdoutsub_full(:,k) = normalize(betaPC_holdoutsub_tStat); % Full map (normalized t stats)
end

% -- Initiate variables.
% This will store the correlation between the projected PC score and observed (ΔGBC ⋅ β_PC3 GBC) for 436 subjects, across all models (with P_select from 718 to 1)
dotprodkcorr_proj_obsv = zeros(1,717);
dotprodkcorr_predobsv  = zeros(1,717);
nrange = 0:717;

% -- Select neural parcels as features
for n = nrange;
	% First, choose which parcels are in Pselect
	Pselect = zeros(1,718); 
	% Select parcels NOT in the n parcels with minimum β values
	for parcel = 1:718
		for minpar = betaPC_abs_sorted(1:n+1)
			if abs(betaPC_abs(parcel)) ~= minpar
				Pselect(1,parcel) = 1;
			else
				Pselect(1,parcel) = 0;
			end
		end
	end
	
	% -- Select the ΔGBC for all Pselect parcels, for all subjects
	gbc_Pselect = repmat(Pselect, nsub, 1) .* gbcdata;
	gbc_Pselect( :, ~any(gbc_Pselect,1) ) = [];  % remove unselected columns
	gbc_Pselect_size = size(gbc_Pselect);

	% -- Select the β_PC for all Pselect parcels
	betaPC_Pselect =  Pselect .* betaPCdat' ;
	betaPC_Pselect( :, ~any(betaPC_Pselect,1) ) = [];

	% -- These will store the (ΔGBC ⋅ β_PC3 GBC)
	dpGBC_Pselect= zeros(nsub,1);

	% -- This will store predicted and observed dpGBC
	dpGBC_pred_all=zeros(nsub,1);
	dpGBC_obsv_all=zeros(nsub,1);

	% Select holdout subject k
	for k = 1:nsub 
		% Remove the parcels not in Pselect
		betaPC_holdoutsub = Pselect .* (betaPC_holdoutsub_full(:,k)'); % select the β_PC map that was computed without the holdout subject
		betaPC_holdoutsub( :, ~any(betaPC_holdoutsub,1) ) = [];  % select only Pselect for that β_PC map

		% Compute observed dpGBC for subject k, using Pselect parcels from the β_PC map that was computed without subject k
		dpGBC_Pselect(k,:) = dot(gbc_Pselect(k,:), betaPC_holdoutsub); %436 by 1 matrix of (ΔGBC ⋅ β_PC3 GBC) per subject
	end	

	for k = 1:nsub 
		trainind = sub_vector(1:end ~= k); % Indices for training subjects (all subjects except subject k)
		dpGBCtrain = dpGBC_Pselect(trainind); % Select the dpGBC for training subjects
		
		pctrain=pc(trainind,:); % Select the PC symptom scores for training subjects
		pctest=pcproj(k); % Select the projected PC score for subject k

		% -- Run linear regression of (ΔGBC ⋅ β_PC3 GBC) on PC score in N-1 training subjects
		y=dpGBCtrain;
		X = [pctrain];		
		b = X\y;
		
		% -- Predict ΔGBC for holdout subject k, via regression
		dpGBCpred = b(1,1) * pctest; 
		
		% -- Save the predicted and observed ΔGBC for k
		dpGBC_pred_all(k,1) = dpGBCpred;

		% dpGBC normalize by no. of parcels in Pselect
		dpGBC_obsv_all(k,1)=dot(gbc_Pselect(k,:)', betaPC_holdoutsub')/gbc_Pselect_size(2);
	end

	dotprodkcorr_predobsv(1,n+1)  = corr(dpGBC_obsv_all, dpGBC_pred_all); % This is the correlation across 436 subjects between observed dpGBC and predicted dpGBC, for Pselect = n	
	dotprodkcorr_proj_obsv(1,n+1) = corr(dpGBC_obsv_all, pcproj'); % This is the correlation across 436 subjects between observed dpGBC and predicted PC3 symptom score, for Pselect = n	
end

% This is the correlation across 436 subjects between observed dpGBC and predicted dpGBC, for all Pselect models
%dlmwrite('BSNIP_PC3_dpGBCPred_dpGBCObsv_r.txt',dotprodkcorr_predobsv','delimiter','\t')  

% This is the correlation across 436 subjects between observed dpGBC and predicted PC3 symptom score, for all Pselect models
%dlmwrite('BSNIP_PC3_dpGBCObs_PC3ScorePred_r.txt',dotprodkcorr_proj_obsv','delimiter','\t') 