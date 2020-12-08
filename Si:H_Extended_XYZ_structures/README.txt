"Training_Data_Revised_Stress" contains all of the base structures used to create "GAP_glue_soap_v5", the first successful GAP for Si:H
After this, iterative training was implemented. Each round of iterative training has two folders:
Iterative_Training_round*: contains the results of using each new version of the potential to conduct MD
selectedStructures_round*: contains those structures which GAP fit poorly enough that it was deemed important to include them in the next training set. Typically 5-10 structures
