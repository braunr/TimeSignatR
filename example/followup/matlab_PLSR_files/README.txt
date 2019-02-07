This directory contains input and output files for use with Laing &al's PLSR code, for the purpose of comparing TimeSignature to the PLSR and differential PLSR predictions.

Laing &al's PLSR code may be obtained from
   http://sleep-sysbio.fhms.surrey.ac.uk/PLSR_16/Code/
Relevant MatLab scripts by Laing &al are in 
   PLSR_16/Code/PhasePrediction/PLSR_MatlabScripts
of the uncompressed zipfile.

=== INPUT DATA ===

The following are input files comprising the training and validation cohorts used in Braun &al 2018, Braun &al 2019:

* CPtrain.mat = training data (use in place of train.mat in Laing &al's scripts)
* CPtrainPairs.mat = training subject pair data (use in place of trainPairs.mat in Laing &al's scripts)
* CPall.mat = validation data (use in place of all.mat in Laing &al's scripts)
* CPallPairs.mat = validation subject pair data (use in place of allPairs.mat in Laing &al's scripts)

=== OUTPUT DATA ===

The following are output files generated using Laing &al's matlab
scripts with the training and validation datasets from Braun &al 2018,
Braun &al 2019.  Note that these can be reproduced using Laing &al's
code (from http://sleep-sysbio.fhms.surrey.ac.uk/PLSR_16/Code/) with
the input files listed above:

* CP_Prediction_PLSR_1sample_5_100.csv = output from running:
    1) GenerateModel_PLSR.m using CPtrain.mat 
    2) Validation_PLSR.m using CPall.mat

* CP_Prediction_PLSR_2samples_5_100.csv = output from running:
    1) GenerateModel_PLSR_2samples.m using CPtrain.mat and CPtrainPairs.mat
    2) Validation_PLSR_2samples.m using CPall.mat and CPallPairs.mat


