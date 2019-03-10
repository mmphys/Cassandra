Cassandra

To run the project:

Run the scripts B, J and S which outputs to Pass1

Delete junk from Pass1 directory
•  i.e. everything EXCEPT .F2 and .log files
•  I started writing this as monolithic program, but moved model creation to pass 2 without tidying up (deleting) model creation is pass 1

Create Pass2 input directory
•  For each .F2 file created in pass 1:
•  Make two copies of the file in Pass2 directory, one for Deuteron data only and one for Proton data only
•  Yes, you need to manually edit the files (sorry!):
•  1. Delete the data you don't want (Deuteron or Proton data)
•  2. Append the target name to the end of the first line in the header, e.g. Cassandra Prediction BCDMS F2 Proton
•  3. Update the Count on line 5 to the number of rows in the file
•  Update the file names to match the convention Pass 2 expects, i.e.:
Deuteron_High_BCDMS_F2
Deuteron_Low_BCDMS_F2
Deuteron_Low_JLAB_F2
Proton_High_BCDMS_F2
Proton_Low_BCDMS_F2
Proton_Low_JLAB_F2
NNPDF/Deuteron_High_SLAC_NNPDF_F2
NNPDF/Deuteron_Low_SLAC_NNPDF_F2
NNPDF/Proton_High_SLAC_NNPDF_F2
NNPDF/Proton_Low_SLAC_NNPDF_F2
reanalyzed/Deuteron_High_SLAC_reanalyzed_F2
reanalyzed/Deuteron_Low_SLAC_reanalyzed_F2
reanalyzed/Proton_High_SLAC_reanalyzed_F2
reanalyzed/Proton_Low_SLAC_reanalyzed_F2

Run Pass2

•  Edit the script T to choose SLAC_TYPE, DATASET and TARGET, then run the script T
•  There are plenty of other parameters in the script you can change - once you are up-to-speed
•  E.g. you can specify a specific model (created in Pass1 at one specific cut) to run against every cut in pass 2 by uncommenting this line:
#MODEL_NUM="_5"

NB: There are more options in the scripts, and far more options supported by the programs
