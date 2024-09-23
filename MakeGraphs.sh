#!/usr/bin/env bash

#set -x
set -e

# User input
OutDir="${OutDir-Graphs}"
RealPath="$(which grealpath)"
Results="${Results-Results}"

# Validation
if [ -z "$RealPath" ]; then
  echo "This script requires grealpath (GNU utilities realpath)"
  exit 1
fi

# Make sure destination exists
mkdir -p "$OutDir"

# Make a link in the destination directory - careful, this overwrites!
function MakeLink()
{
  local Source="$Results/$1"
  local DestName="$2"
  local LinkSource="$("$RealPath" --relative-to="$OutDir" "$Source")"
  ln -sf "$LinkSource" "$OutDir/$DestName"
}

# Make a series of links
function MakeSeries()
{
  local Prefix="$1"
  local DestPrefix="$2"
  local Suffix Dest
  for Suffix in F2DTvsQ2 F2DTvsX TheoryExpt CorrelTheory q_up q_down q_gluon q_ubar q_dbar
  do
    [ -n "$DestPrefix" ] && Dest="${DestPrefix}_${Suffix}.pdf"
    MakeLink "${Prefix}_${Suffix}.pdf" "$Dest"
  done
}

# Main loop - Make an entry for each of the graphs in the dissertation (in diss'n order)
# Overall dataset stats - S3.2
MakeLink NNPDF/PAN_DataSet.pdf DataSet_NNPDF_All.pdf
MakeLink NNPDF/PHN_DataSet.pdf DataSet_NNPDF_High.pdf
MakeLink NNPDF/PLN_DataSet.pdf DataSet_NNPDF_Low.pdf
MakeLink reanalyzed/PAr_DataSet.pdf DataSet_reanalyzed_All.pdf
MakeLink reanalyzed/PHr_DataSet.pdf DataSet_reanalyzed_High.pdf
MakeLink reanalyzed/PLr_DataSet.pdf DataSet_reanalyzed_Low.pdf
# NNPDF SLAC data - S4.1 & S4.2
MakeLink NNPDF/Proton_Low/PLN.chi.pdf
MakeLink NNPDF/Deuteron_Low/DLN.chi.pdf
MakeLink NNPDF/Proton_All/PAN.chi.pdf
MakeLink NNPDF/Deuteron_All/DAN.chi.pdf
MakeSeries NNPDF/Proton_All/PAN_5
MakeSeries NNPDF/Deuteron_All/DAN_5
# Reanalysed SLAC data - S4.2.2
MakeLink reanalyzed/Proton_Low/PLr.chi.pdf # Not referred to in MSc
MakeLink reanalyzed/Deuteron_Low/DLr.chi.pdf # Not referred to in MSc
MakeLink reanalyzed/Proton_All/PAr.chi.pdf
MakeLink reanalyzed/Deuteron_All/DAr.chi.pdf
MakeSeries reanalyzed/Proton_All/PAr_6
MakeSeries reanalyzed/Deuteron_All/DAr_6
# JLAB higher twist model with reanalysed SLAC data - S4.3
MakeLink 0/Proton_All/PA0.chi.pdf PAJ.chi.pdf
MakeLink 0/Deuteron_All/DA0.chi.pdf DAJ.chi.pdf
MakeLink reanalyzed/Proton_All_MJ/PArMJ.chi.pdf PA0.chi.pdf
MakeLink reanalyzed/Deuteron_All_MJ/DArMJ.chi.pdf DA0.chi.pdf
MakeSeries reanalyzed/Proton_All_MJ/PArMJ_6 PA0_6
MakeSeries reanalyzed/Deuteron_All_MJ/DArMJ_6 DA0_6
# NNPDF and reanalysed higher kinematic cuts - S4.4
MakeLink NNPDF/Proton_High/PHN.chi.pdf
MakeLink reanalyzed/Proton_High/PHr.chi.pdf
MakeSeries NNPDF/Proton_High/PHN_5
MakeSeries reanalyzed/Proton_High/PHr_9
# ditto deuteron - not referred to in MSc
MakeLink NNPDF/Deuteron_High/DHN.chi.pdf
MakeLink reanalyzed/Deuteron_High/DHr.chi.pdf
MakeSeries NNPDF/Deuteron_High/DHN_5
MakeSeries reanalyzed/Deuteron_High/DHr_9
# Appendix
MakeLink NNPDF/Proton_All/PAN_5_CorrelExpt.pdf
MakeLink NNPDF/Proton_All/PAN_5_CorrelExptTheory.pdf
