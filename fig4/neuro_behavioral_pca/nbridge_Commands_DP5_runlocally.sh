#!/bin/sh
# set -x



# ------------------------------------------------------------------------------
# -- Running the N-BRIDGE Container
# ------------------------------------------------------------------------------


nbridgeversion="0_9_93"
nbridgedatadir="/Users/floramoujaes/n-bridge/nbridgeExample/DP5_delta_delta_CAN-NP_V2"
# mkdir "${nbridgedatadir}"
docker container run -it -v "${nbridgedatadir}":/nbridgeExampleOutput docker.io/nbridge/n-bridge:${nbridgeversion} bash 


# ------------------------------------------------------------------------------
# -- N-BRIDGE Prep Example Composite Score run 30.03.2020
# ------------------------------------------------------------------------------
# -- Define general prep variables
StudyName="DP5_delta_delta_CAN-NP_V2"
OutPath="/nbridgeExampleOutput"
unset OutDataPrep
# The possible outputs of Prep are documented in the inline help and the wiki
OutDataPrep="bcovmtx,bcormtx,bcorchart,bgroupfxbar,bgroupfxridge,bpcaload,bpcascores,bpcascorescon,bpcatriplot,bpcascree,bpcabar,bpcapie,bpcaradar,bpcaridge,bpcacormtx,bicascores,bicaload,bicacormtx,bicapie,bicascree"
PermPrep="5000" #number of permutations to run for determining null distributions

# -- Define variables for example with BSNIP Data:
#
# -- Define variables for example with BSNIP Data:
# input is DataDx = Ket, Data CON = Placebo
InputBehavior="/nbridgeExampleOutput/data/DP5_Behaviour_Delta_ordered.tsv"
ControlBehavior="/nbridgeExampleOutput/data/DP5_BehavioralDataCON_ordered.tsv"
OutPrefix="DP5_delta_delta_CAN-NP_V2"
DropZeroSD="FALSE" # This option keeps all of the behavioral items in the geometry which may be desirable if other datasets are to be projected

# only using 1 composite cognition score so change cog to 1
bash nbridge.sh \
--nbridge_function="prep" \
--behaviorcat="Cognitive:1,Positive:7,Negative:7,General:16" \
--outdata="${OutDataPrep}" \
--outpath="${OutPath}" \
--inputbehavior="${InputBehavior}" \
--controlbehavior="${ControlBehavior}" \
--permutations="${PermPrep}" \
--outprefix="${OutPrefix}" \
--dropzerosd="${DropZeroSD}" \
--createreport="YES"




# ------------------------------------------------------------------------------
# -- N-BRIDGE initialize ran on grace with separate matlab script see /SAY/standard8/Grace-CC1096-MEDPSYPsychDivisions_CNRU/Studies/Anticevic.DP5/analysis/results_PALM_GBC_40/Parcellation_NBridge

# ------------------------------------------------------------------------------
# -- Define general initialize variables

 OutPath="/nbridgeExampleOutput"
 OutDataInit="univarmap"
 #OutDataInit="ccascree,ccapie,ccaridge,ccabweights,ccabradar,ccaloadings,ccaweights,ccascores,ccaloadshuffle,ccaloopred,ccascorecor"
 PermInit="5000"
 OutPrefix="DP5_delta_delta_CAN-NP_V2"
 ### Note that OutPrefix (the input for the --outprefix parameter contains <dataset>_<feature_set>_<unique_run_name> per the specification)

 # Neural data in PTSERIES format:
 ControlNeural="/nbridgeExampleOutput/data/Placebo_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz.ptseries.nii"
 SubjectNeural="/nbridgeExampleOutput/data/Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz.ptseries.nii"
 
 BehPCALoadings="/nbridgeExampleOutput/analysis/results/NBRIDGE_DP5_delta_delta_CAN-NP_V2_BehaviorPCALoadings.tsv"
 InputScores="/nbridgeExampleOutput/analysis/results/NBRIDGE_DP5_delta_delta_CAN-NP_V2_BehaviorPCAScores.tsv"
 ControlScores="/nbridgeExampleOutput/analysis/results/NBRIDGE_DP5_delta_delta_CAN-NP_V2_BehaviorPCAScoresCON.tsv"
 
 
# If the input behavioral data has been dimensionality-reduced via nbridgePrep, specify --dimreduct=TRUE and set --behpcaloadings 
# to the output of --bpcaload from nbridgePrep.
 
# # Univariate mapping -- note that "--symmetrize" needs to be set to "no" for subcortical output

unset OutDataInit
OutDataInit="univarmap"
OutputSuffix="PC"
nbridge.sh \
--nbridge_function="initialize" \
--outdata="${OutDataInit}" \
--outpath="${OutPath}" \
--inputbehavior="${InputScores}" \
--subjectneural="${SubjectNeural}" \
--symmetrize="no" \
--permutations="${PermInit}" \
--outprefix="${OutPrefix}" \
--outsuffix="${OutputSuffix}" \
--createreport="YES"
 
# -> review in figures BSNIP_BACSPANSS_2BACS_NBRIDGEinitialize_AllFigures_2019-12-02_12.19.05.pdf does not have all panss -> now it does (for univarmap this does not apply)







