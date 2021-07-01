% Compute correlation between PC maps from 1k runs of split-half replication
% matlab -nodisplay

% -- Initiate matrices
H1_obspred_r=zeros(5000,1);
H2_obspred_r=zeros(5000,1);

% -- Run for all in H1
half=1 
for pc = 1:5
  for run=1:1000
    % The inputs are Z-scored beta coefficient files (output from PALM, see associated bash script)
    obsv=gmrimage(sprintf('BSNIP_StratifiedPRB_H%dScores_Run%d_PC%d_dat_ztstat.pscalar.nii', half, run, pc));
    pred=gmrimage(sprintf('BSNIP_StratifiedPRB_H%dPredScores_Run%d_PC%d_dat_ztstat.pscalar.nii', half, run, pc));

    % -- save these in a matrix
    H1_obspred_r(1000 * (pc-1) + run)=corr(obsv.data, pred.data); % all runs, then all PCs for one Half
    corr(obsv.data, pred.data)
  end
end

% Rerun for H2 to get a separate data matrix
half=2 
for pc = 1:5
  for run=1:1000
    obsv=gmrimage(sprintf('BSNIP_StratifiedPRB_H%dScores_Run%d_PC%d_dat_ztstat.pscalar.nii', half, run, pc));
    pred=gmrimage(sprintf('BSNIP_StratifiedPRB_H%dPredScores_Run%d_PC%d_dat_ztstat.pscalar.nii', half, run, pc));

  % -- save these in a matrix
    H2_obspred_r(1000 * (pc-1) + run)=corr(obsv.data, pred.data); % all runs, then all PCs for one Half
    corr(obsv.data, pred.data)
  end
end


dlmwrite('H1_ObsPred_betaPCGBC_Correlations.txt', H1_obspred_r)
dlmwrite('H2_ObsPred_betaPCGBC_Correlations.txt'', H2_obspred_r)