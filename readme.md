Cassandra
=========

There are two project executables:

1. **Cassandra** (pass 1) Add F2 predictions to DIS data
2. **Twiggy** (pass 2) Model creation / validation

## **Cassandra** (pass 1)

Run the scripts B, J and S which outputs to Pass1

## Manually create input for pass 2

Manually create the Pass2 directory. For each .F2 file created in pass 1, make two copies of the file in Pass2 directory, one for Deuteron data only and one for Proton data only. Update the file names to match the convention for pass 2, i.e.:

Bullet list:

* Deuteron_High_BCDMS_F2
* Deuteron_Low_BCDMS_F2
* Deuteron_Low_JLAB_F2
* Proton_High_BCDMS_F2
* Proton_Low_BCDMS_F2
* Proton_Low_JLAB_F2
* NNPDF/Deuteron_High_SLAC_NNPDF_F2
* NNPDF/Deuteron_Low_SLAC_NNPDF_F2
* NNPDF/Proton_High_SLAC_NNPDF_F2
* NNPDF/Proton_Low_SLAC_NNPDF_F2
* reanalyzed/Deuteron_High_SLAC_reanalyzed_F2
* reanalyzed/Deuteron_Low_SLAC_reanalyzed_F2
* reanalyzed/Proton_High_SLAC_reanalyzed_F2
* reanalyzed/Proton_Low_SLAC_reanalyzed_F2

Manually edit the files:

Numbered list:

1. Delete the data you don't want (Deuteron or Proton data)
2. Append the target name to the end of the first line in the header, e.g. `Cassandra Prediction BCDMS F2 Proton`
3. Update the Count on line 5 to the number of rows in the file

## **Twiggy** (pass 2)

Edit the script T to choose SLAC_TYPE, DATASET and TARGET, then run the script T

There are plenty of other parameters in the script you can change - once you are up-to-speed, e.g. you can specify a specific model (created in Pass1 at one specific cut) to run against every cut in pass 2 by uncommenting this line:

`#MODEL_NUM="_5"`

NB: There are more options in the scripts, and far more options supported by the executables.
