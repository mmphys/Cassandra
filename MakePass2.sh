#!/usr/bin/env bash

OutDir=Pass2

# Read InFile extracting Nucleon records

function ProcessFile()
{
(
  # Get name of OutFile
  local OutFile="$OutDir/"
  if [ -v Suffix ]; then
    OutFile+="$Suffix/"
  fi
  mkdir -p "$OutFile"
  OutFile+="${Nucleon}_${HighLow}_${Type}"
  [ -v Suffix ] && OutFile+="_$Suffix"
  OutFile+="_F2"
  # Make sure it doesn't exist
  if [ -e "$OutFile" ]; then
    echo " $OutFile exists - skipping"
    return
  fi
  # How many records are there
  Count=$(awk "/$Nucleon/{++n}; END {print n}" "$InFile")
  echo " $Count $OutFile"
  # Write the header
  awk "/Count/{print \"             Count : $Count\";next};{print};!NF{exit}" \
      "$InFile" > "$OutFile"
  # Now write selected records
  awk "/$Nucleon/ {print}" "$InFile" >> "$OutFile"
)
}

# Main loop - process each file, breaking into parts

SourceSuffix="_0_Cassandra.F2"
for InFile in Pass1/*"$SourceSuffix"
do
  Warn=1
  File="${InFile##*/}"
  FileParts="${File%$SourceSuffix}"
  Parts=(${FileParts//_/ })
  if((${#Parts[@]}>=2)); then
    Prefix="${Parts[0]}"
    Type="${Parts[-1]}"
    case "$Prefix$Type" in
      BBCDMS | JJLAB | SSLAC) Warn=0;;
    esac
    if((!Warn)); then
      unset Parts[0]
      unset Parts[-1]
      if (( ${#Parts[@]} )); then
        HighLow="${Parts[-1]}"
        unset Parts[-1]
        case "$HighLow" in
          High | Low) : ;;
          *) Warn=1;;
        esac
      else
        HighLow=Low
      fi
      if((!Warn)); then
        if (( ${#Parts[@]} )); then
          Suffix="${Parts[-1]}"
          unset Parts[-1]
        else
          unset Suffix
        fi
        echo "$File"
        for Nucleon in Proton Deuteron; do ProcessFile; done
      fi
    fi
  fi
  ((Warn)) && echo "$File - unknown type"
done
