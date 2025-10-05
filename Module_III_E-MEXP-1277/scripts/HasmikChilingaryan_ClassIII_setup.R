# HasmikChilingaryan_ClassIII_setup.R
# Setup script: inspect raw files and SDRF

# Setup Project folders
dir.create("Module_III_E-MEXP-1277", showWarnings = FALSE)
dir.create("Module_III_E-MEXP-1277/raw_data", recursive = TRUE, showWarnings = FALSE)
dir.create("Module_III_E-MEXP-1277/scripts", showWarnings = FALSE)
dir.create("Module_III_E-MEXP-1277/results", showWarnings = FALSE)
dir.create("Module_III_E-MEXP-1277/figures", showWarnings = FALSE)

list.dirs("Module_III_E-MEXP-1277", recursive = FALSE) #confirm

list.files("AI_BiotechBioinf_Internship2025/Module_III_E-MEXP-1277/raw_data") #confirm file transfer
