## Process this file with automake to produce Makefile.in

# Background rate files to distribute
bkgs = bkg_kb_A_50h_v3.txt \
       bkg_kb_B_50h_v3.txt \
       bkg_kb_C_50h_v3.txt \
       bkg_kb_D_50h_v3.txt \
       bkg_kb_E_50h_v3.txt \
       bkg_kb_F_50h_v3.txt \
       bkg_kb_G_50h_v3.txt \
       bkg_kb_H_50h_v3.txt \
       bkg_kb_I_50h_v3.txt \
       bkg_kb_J_50h_v3.txt \
       bkg_kb_K_50h_v3.txt \
       bkg_kb_NA_50h_v3.txt \
       bkg_kb_NB_50h_v3.txt
bkgsdir   = $(datadir)/models
bkgs_DATA = $(bkgs)

# Models to distribute
models      = crab.xml \
              disk.xml \
              gauss.xml \
              shell.xml
modelsdir   = $(datadir)/models
models_DATA = $(models)

# Extra files to be distributed in the tar file
EXTRA_DIST = $(bkgs) $(models)
