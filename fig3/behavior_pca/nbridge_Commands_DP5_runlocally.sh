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






