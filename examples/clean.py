#===============================
# This is a file to clean up the output info from GENOA
# Author: Zhizhao Wang
# Time: 
#===============================

# Please put 1 for the part need to be cleaned
options={
          'SSH':1; # 1 for remove the generated IDchem_ref and IDchem_rdc
          'SOA Conc.':1; # remove the generated SOA concentrations (except for the final reduction and the one with the named preserved in the variable IDkeep)
          'Chems':1; # remove the mechanism and the record files generated for each reduction step (except for the final reduction and those names kept in IDkeep)

IDkeep = []

        }

# clean 
