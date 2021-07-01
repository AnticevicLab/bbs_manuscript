% These were performed using PALM in Matlab
%matlab -nodisplay

%%% -- PC Scores
# -- PSD Subjects
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPC1.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PC1
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPC2.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PC2
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPC3.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PC3
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPC4.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PC4
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPC5.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PC5

% -- CON Subjects
palm -i ../data_neural_gbc/BSNIP_CON_N202_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n202_dPC1.csv -t BSNIP_n202_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o CON_GBC_PC1
palm -i ../data_neural_gbc/BSNIP_CON_N202_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n202_dPC2.csv -t BSNIP_n202_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o CON_GBC_PC2
palm -i ../data_neural_gbc/BSNIP_CON_N202_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n202_dPC3.csv -t BSNIP_n202_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o CON_GBC_PC3
palm -i ../data_neural_gbc/BSNIP_CON_N202_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n202_dPC4.csv -t BSNIP_n202_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o CON_GBC_PC4
palm -i ../data_neural_gbc/BSNIP_CON_N202_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n202_dPC5.csv -t BSNIP_n202_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o CON_GBC_PC5

%%% -- BACS/PANSS Factor Scores
% -- PSD Subjects
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPANSSNeg.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PANSS_Neg
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPANSSPos.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PANSS_Pos
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dPANSSGen.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_PANSS_Gen
palm -i ../data_neural_gbc/BSNIP_PSD_N436_GSR.udvarsme.surface_gbc_mFz.dtseries.nii -d BSNIP_n436_dBACSComp.csv -t BSNIP_n436_t1.csv -twotail -zstat -accel tail -n 5000 -transposedata -demean -o GBC_BACS_Cog