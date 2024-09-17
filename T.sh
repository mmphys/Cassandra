#!/bin/bash

DATAROOT="Pass2"
RESULTSROOT="Results"

#DATASET can be Low, High or All (which processes all data)
#TARGET can be Proton, Deuteron or All (which processes all data)
#SLAC_TYPE selects whether/(which version of) the SLAC data are included:
# i.e. NNPDF, reanalyzed or 0 (Don't include any SLAC data)

function Process()
{
#Validate Parameters (DATASET, TARGET and SLAC_TYPE)
SetTarget=0
if [ "$DATASET" == High ] || [ "$DATASET" == Low ]
then
    ((SetTarget++))
else
    DATASET=All
fi
if [ "$TARGET" == Proton ] || [ "$TARGET" == Deuteron ]
then
    ((SetTarget++))
else
    TARGET=All
fi
if [ "$SLAC_TYPE" != reanalyzed ] && [ "$SLAC_TYPE" != NNPDF ]
then
    SLAC_TYPE=0
fi

local PARAMS="-l4 -bx0.05 -bq2 -o3 -x0 -p"

# Select an appropriate cutoff
local CUTOFF
if [ "$DATASET" == High ]
then
    CUTOFF="-wl12.5 -wu17 -n9"
else
    if [ "$SLAC_TYPE" == 0 ]
    then
	CUTOFF="-wl2.5 -wu3.7 -n12"
    elif [ "$SLAC_TYPE" == reanalyzed ]
    then
	CUTOFF="-wl3 -wu12.5 -n19"
    else
	CUTOFF="-wl5 -wu12.5 -n15"
	#CUTOFF="-wl5 -wu12.5 -n30"
	#CUTOFF="-wl1 -wu12.5 -n46"
    fi
fi

# Select an appropriate model
local MODEL
if [ "$DATASET" == Low ]
then
    if [ "$TARGET" == Proton ]
    then
	# I'm judging by cut #7 because chi low
	#MODEL="-ad0.,0.,1"
	#MODEL="-aa0.75,1.0 -ab3.5,5.5 -ac-1.4,0. -ad0.,0.,1"
	# The following line is best proton fit in September
	#MODEL="-tm0 -aa0.87,1.0 -ab5.15,5.3 -ac-1.3,-0.45 -ad0.,0.,1"

	#Testing 1-Aug-2018
	#MODEL="-tm1 -ad-2.,-0.2"
	#MODEL="-tm1 -aa0.5,1 -ab0.,0.5 -ac-1.2,-0.6 -ad-2.,-0.5"
	#MODEL="-tm1 -aa0.65,1.15 -ab0.,0.5 -ac-1.2,-0.6 -ad-2.4,-1"
	#MODEL="-tm1 -aa0.75,1.3 -ab0.,1 -ac-1.6,-0.7 -ad-2.6,-1.2"
	#Best fit as of 3 Aug
	#MODEL="-tm1 -aa0.9,1.4 -ab0.,1 -ac-2,-0.9 -ad-2.8,-1.2"

	#Testing 4-Aug-2018
	#MODEL="-tm0 -aa0.87,1.0 -ab5.15,5.3 -ac-1.3,-0.45 -ad0.,0.,1"

	#MODEL="-ta1 -ac-0.8,0.8,21 -ad-2.,-0.2"
	#MODEL="-ta1 -aa0.5,1 -ab0,6 -ac-0.75,-0.1 -ad-2.,-0.5"
	#MODEL="-ta1 -aa0.7,1 -ab0.4,4.4 -ac-0.7,-0.15 -ad-2.,-1.3"
	#MODEL="-ta1 -aa0.7,1 -ab0.5,3 -ac-0.6,-0.2 -ad-2.,-1.45"
	#Best model 4-Aug-2018
	MODEL="-ta1 -aa0.7,1 -ab0.6,2.3 -ac-0.6,-0.2 -ad-2.,-1.45"

	#MODEL="-ta2 -ab-2,0,11 -ac-0.8,0.8,11 -ad0,0"
    elif [ "$TARGET" == Deuteron ]
    then
	#MODEL="-ad0.,0.,1"
	#MODEL="-aa0.8,1 -ab5,10 -ac-1.5,-0.5 -ad0.,0.,1"
	# The following line is best deuteron fit to date
	#MODEL="-aa0.9,1 -ab4.5,7 -ac-1.15,-0.75 -ad0.,0.,1"
	#11/08/2018
	#MODEL="-ta1 -aa0.7,1 -ab0.6,2.3 -ac-0.6,-0.2 -ad-2.,-1.45"
	MODEL="-ta1 -aa0.85,1 -ab0.8,1.8 -ac-0.5,-0.25 -ad-2.,-1.7"
    else # All
	MODEL="-ad0.,0.,1"
    fi
else
    # Use pre-built HT model
    #MODEL_NUM="_5"
    MODEL="-ta1 -m${RESULTSROOT}/${SLAC_TYPE}/${TARGET}_Low/${TARGET:0:1}L${SLAC_TYPE:0:1}${MODEL_NUM}"
    # This stops the graphs blowing up
    #CUTOFF="${CUTOFF} -xc0.11"
fi

# Translate parameters to source files
local DATA_FILE_NAME
if [ $SetTarget -eq 2 ]
then
   DATA_FILE_NAME="${TARGET}_${DATASET}_*"
elif [ $SetTarget -eq 0 ]
then
     DATA_FILE_NAME="*"
elif [ "$DATASET" == All ]
then
     DATA_FILE_NAME="${TARGET}*"
else
     DATA_FILE_NAME="*_${DATASET}*"
fi

local DATA="$DATAROOT/${DATA_FILE_NAME}"
[ $SLAC_TYPE != 0 ] && DATA+=" $DATAROOT/${SLAC_TYPE}/${DATA_FILE_NAME}"

# Make sure the results directory exists
local Dir="$SLAC_TYPE/${TARGET}_${DATASET}"
mkdir -p "$RESULTSROOT/$Dir"

# Now run
local Combo="${TARGET:0:1}${DATASET:0:1}${SLAC_TYPE:0:1}"
local File="$RESULTSROOT/$Dir/$Combo"
echo "Performing Combined $SLAC_TYPE $TARGET $DATASET run"
local Cmd="./Twiggy '-l$File' $DATA $PARAMS $MODEL $CUTOFF"
echo " $Cmd"
eval "$Cmd"
(
  cd "${RESULTSROOT}"
  tar -czf "$Combo.tar" "$Dir/"
)
}

# Main

for DATASET in Low High All
do
for TARGET in Proton Deuteron #All
do
for SLAC_TYPE in NNPDF reanalyzed 0
do
  Process &
done
done
done
wait
