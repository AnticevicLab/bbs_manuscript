# -- gmri is a tool from QuNex (www.qunex.yale.edu), which offers scheduler support for PALM (a Matlab/Octave based tool).
# This code was originally run with an earlier version of QuNex circa 2019. 
# QuNex is now available for download and will be publically released in 2021. 
# Please refer to the latest version of the documentation for updated usage options

pcs="1 2 3 4 5"
runs=`seq 1000`
halves="H1 H2"

# -- This loop selects the subjects assigned to each half in each run, and saves it as a separate file
for half in ${halves}; do
    for run in ${runs}; do
        a="wb_command -cifti-merge BSNIP_StratifiedPRB_${half}_Run${run}_GSR.udvarsme.CAB-NP-718_gbc_mFz_.ptseries.nii -cifti ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz.ptseries.nii"
        c=""
        while IFS= read -r line; do c="$c -column $line "; done < BSNIP_StratifiedPRB_${half}Indices_Run${run}.txt
        eval "$a $c"
    done
done

# -- This loop will run the regression of PC scores to GBC for each half, for each run of the split-half CCA
for half in $halves; do
    for run in $runs; do
        for pc in $pcs; do
            # -- For Observed scores
            gmri run_palm image="BSNIP_StratifiedPRB_${half}_Run${run}_GSR.udvarsme.CAB-NP-718_gbc_mFz_.ptseries.nii" design="name:BSNIP_StratifiedPRB_${half}Scores|d:dPC${pc}Run${run}|t:t1" overwrite=yes args="n:3|accel:tail|twotail|zstat|demean|transposedata" root="BSNIP_StratifiedPRB_${half}Scores_Run${run}_PC${pc}" 
            # -- For Predicted scores
            gmri run_palm image="BSNIP_StratifiedPRB_${half}_Run${run}_GSR.udvarsme.CAB-NP-718_gbc_mFz_.ptseries.nii" design="name:BSNIP_StratifiedPRB_${half}PredScores|d:dPC${pc}Run${run}|t:t1" overwrite=yes args="n:3|accel:tail|twotail|zstat|demean|transposedata" root="BSNIP_StratifiedPRB_${half}PredScores_Run${run}_PC${pc}" 
        done
    done
done