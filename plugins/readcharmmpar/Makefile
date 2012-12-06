.SILENT:

VMFILES = readcharmmpar.tcl par_all27_prot_lipid_na.inp pkgIndex.tcl par_amber2charmm.inp par_all36_lipid.prm
VMVERSION = 1.2
DIR = $(PLUGINDIR)/noarch/tcl/readcharmmpar$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying readcharmmpar $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp $(VMFILES) $(DIR) 

	
