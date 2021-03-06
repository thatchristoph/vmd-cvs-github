#
# vmd documentation Makefile
#
.SILENT:

############################################################################
##################  no more configuration beyond this point  ###############
############################################################################

# makefile configuration
#VPATH                   = .:../src
SHELL 			= /bin/sh
.SUFFIXES: .C .c ..c .i .o

# directories for object files, scripts, and different parts of code.
# these are all relative to being in any of the subdirs of the main directory
# such as the src, lib, or bin dirs
MAINDIR			= ../.

#####
# specific names of programs used in make commands
#####
GLINCDIR	= 
GLLIBDIR	= 
GLLIBNAME	= 
CD		= cd
COPY		= cp
COPYDIR		= cp -r
DOXYGEN		= doxygen
MAKEDIR		= mkdir
MOVE		= mv -f
DELETE		= rm -f
DELETEDIR	= rm -rf
ECHO		= echo
TAR		= tar
TOUCH		= touch
GZIP		= gzip
COMPRESS	= compress
GNUCOMPRESS	= /usr/local/bin/gzip
AR		= ar rsc
RANLIB		= /bin/true
COCMD		= co -q
COLOCKCMD	= co -l
CICMD		= ci -m"No message"
MAKEINDEX	= makeindex
LATEX		= latex
LATEX2HTML	= latex2html

# PDF generating version of Latex.
PDFLATEX	= pdflatex

# Adobe distiller (Postscript to PDF converter) for Unix, version 3.x
DISTILL  	= distill

# Free Postscript to PDF converter
PS2PDF		= ps2pdf

# Ghostscript, used to convert PS files to PDF, need version 5.x or better
GHOSTSCRIPT	= gs

###############################
# Documentation files 
###############################
DOC_STYLE_FILES	= \
		graphic.sty \
		picinpar.sty

# files needed for the installation guide
MAIN_INST_DOC	= ig.tex
DOC_INST_FILES	= \
		ig_chapters.tex \
		ig_custom.tex \
		ig_install.tex \
		ig_intro.tex \
		ig_macros.tex 
ALL_INST_FILES	= $(MAIN_INST_DOC) $(DOC_INST_FILES)

# files needed for the user's guide
MAIN_USER_DOC	= ug.tex
DOC_USER_FILES	= \
		ug_basics.tex \
		ug_chapters.tex \
		ug_exec_env.tex \
		ug_file_types.tex \
		ug_form_animate.tex \
		ug_form_color.tex \
		ug_form_display.tex \
		ug_form_edit.tex \
		ug_form_file.tex \
		ug_form_graphics.tex \
		ug_form_labels.tex \
		ug_form_main.tex \
		ug_form_mol.tex \
		ug_form_remote.tex \
		ug_form_render.tex \
		ug_form_tracker.tex \
		ug_forms.tex \
		ug_future.tex \
		ug_graphics.tex \
		ug_intro.tex \
		ug_mdcomm.tex \
		ug_misc.tex \
		ug_mol_colors.tex \
		ug_mol_reps.tex \
		ug_mol_selection.tex \
		ug_molecule_info.tex \
		ug_mouse_display.tex \
		ug_mouse_forms.tex \
		ug_other_ui.tex \
		ug_rendering.tex \
		ug_script_writing.tex \
		ug_stereo.tex \
		ug_text_ui.tex \
		ug_tutorial.tex \
		ug_vectors.tex

ALL_USER_FILES	= $(MAIN_USER_DOC) $(DOC_USER_FILES)

# files needs by all the different guides, and for the book format
MAIN_ALL_DOC	= vmd.tex
DOC_ALL_FILES	= \
		vmd_authors.tex \
		vmd_begindoc.tex \
		vmd_configure.tex \
		vmd_macros.tex \
		vmd_moreinfo.tex \
		vmd_paramcomp.tex \
		vmd_paraminst.tex \
		vmd_ref.tex \
		vmd_style.tex \
		vmd_title.tex

ALL_MAIN_FILES	= $(MAIN_ALL_DOC) $(DOC_ALL_FILES)

DOC_FILES = \
		$(ALL_INST_FILES) \
		$(ALL_USER_FILES) \
		$(ALL_MAIN_FILES) \
		$(DOC_STYLE_FILES)

DOC_DRIVER_FILES	= $(MAIN_INST_DOC) $(MAIN_USER_DOC) $(MAIN_ALL_DOC)
DOC_DVI_FILES		= $(DOC_DRIVER_FILES:.tex=.dvi)
DOC_PS_FILES		= $(DOC_DRIVER_FILES:.tex=.ps)

###############################
### latex the documentation
###############################
default: 
	$(ECHO) "Please choose one of the following forms of documentation:"
	$(ECHO) "----------------------------------------------------------"
	$(ECHO) "User, Installation, and old Programmer Guides:"
	$(ECHO) "  all      -- Build all ps, pdf, and html documentation   "
	$(ECHO) "  ps       -- Adobe Postscript                            "
	$(ECHO) "  pdf      -- Adobe PDF (pdflatex)                        "
	$(ECHO) "  distill  -- Adobe PDF (distiller, translated from PS)   "
	$(ECHO) "  ps2pdf   -- Adobe PDF (ghostscript, translated from PS) "
	$(ECHO) "  html     -- HTML                                        "
	$(ECHO) ""
	$(ECHO) "Programmer Documentation:"
	$(ECHO) "  doxygen  -- Autogenerated source documentation (doxygen)"
	$(ECHO) ""
	$(ECHO) "clean      -- clean out generated files for fresh start  "

all: 
	$(MAKE) ps;
	$(MAKE) html
	$(MAKE) distill

#	$(MAKE) pdf
#	$(MAKE) ps2pdf

ps: doc.ig.ps doc.ug.ps

pdf: doc.ig.pdf doc.ug.pdf

distill: doc.ig.ps.pdf doc.ug.ps.pdf

ps2pdf: doc.ig.ps2pdf doc.ug.ps2pdf

html: doc.ig.html doc.ug.html

#
# Auto-generated source code documentation using "doxygen"
#
doxygen: 
	$(DELETEDIR) dox
	$(DOXYGEN) Doxyfile

#
# Installation guide
#
doc.ig.ps:
	$(ECHO) "Creating the Postscript installation guide ..."
	$(DELETE) ig.ps
	$(LATEX) $(MAIN_INST_DOC) && \
	$(LATEX) $(MAIN_INST_DOC) && \
	$(LATEX) $(MAIN_INST_DOC) && \
		dvips -Pcmz -o -Z $(MAIN_INST_DOC:.tex=.dvi)

doc.ig.pdf:
	$(ECHO) "Creating the PDF (pdflatex) installation guide ..." 
	$(PDFLATEX) $(MAIN_INST_DOC)

doc.ig.ps.pdf: ig.ps
	$(ECHO) "Creating the PDF (distiller) installation guide ..." 
	$(DELETE) ig.pdf
	$(PS2PDF) ig.ps

doc.ig.ps2pdf: ig.ps
	$(ECHO) "Creating the PDF (ghostview) installation guide ..." 
	$(DELETE) ig.pdf
	$(GHOSTSCRIPT) -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="ig.pdf" -c save pop -f "ig.ps" -c quit

doc.ig.html:
	$(ECHO) "Creating the HTML installation guide ..."
	$(DELETEDIR) ./ig
	$(LATEX2HTML) \
		-local_icons \
		-split 5 \
		-toc_depth 4 \
		-toc_stars \
		-t "VMD Installation Guide" \
		-address "vmd@ks.uiuc.edu" \
		-up_url "http://www.ks.uiuc.edu/Research/vmd/" \
		-up_title "Main VMD Page" \
		 $(MAIN_INST_DOC)
	$(TAR) -cvf - ig | $(GZIP) > ig.html.tar.gz 

#
# Users guide
#
doc.ug.ps:
	$(ECHO) "Creating the Postscript users guide ..."
	$(DELETE) ug.ps
	$(LATEX) $(MAIN_USER_DOC) && \
	$(MAKEINDEX) $(MAIN_USER_DOC:.tex=) && \
	$(LATEX) $(MAIN_USER_DOC) && \
	$(MAKEINDEX) $(MAIN_USER_DOC:.tex=) && \
	$(LATEX) $(MAIN_USER_DOC) && \
	$(MAKEINDEX) $(MAIN_USER_DOC:.tex=) && \
	$(LATEX) $(MAIN_USER_DOC) && \
		dvips -Pcmz -o -Z $(MAIN_USER_DOC:.tex=.dvi)

doc.ug.pdf:
	$(ECHO) "Creating the PDF users guide ..."
	$(PDFLATEX) $(MAIN_USER_DOC)

doc.ug.ps.pdf: ug.ps
	$(ECHO) "Creating the PDF (distiller) users guide ..."
	$(DELETE) ug.pdf
	$(PS2PDF) ug.ps

doc.ug.ps2pdf: ug.ps
	$(ECHO) "Creating the PDF (ghostview) users guide ..." 
	$(DELETE) ug.pdf
	$(GHOSTSCRIPT) -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="ug.pdf" -c save pop -f "ug.ps" -c quit

doc.ug.html:
	$(ECHO) "Creating the HTML users guide ..."
	$(DELETEDIR) ./ug
	$(LATEX2HTML) \
		-local_icons \
		-split 5 \
		-toc_depth 4 \
		-toc_stars \
		-t "VMD User's Guide" \
		-address "vmd@ks.uiuc.edu" \
		-up_url "http://www.ks.uiuc.edu/Research/vmd/" \
		-up_title "Main VMD Page" \
		$(MAIN_USER_DOC)
	$(TAR) -cvf - ug | $(GZIP) > ug.html.tar.gz

#
# Complete book of all documentation
#
doc.book.ps:
	$(ECHO) "Creating the complete Postscript documentation book ..."
	$(LATEX) $(MAIN_ALL_DOC) && \
		dvips -Pcmz -o -Z $(MAIN_ALL_DOC:.tex=.dvi)

doc.book.pdf:
	$(ECHO) "Creating the complete PDF documentation book ..."
	$(PDFLATEX) $(MAIN_ALL_DOC)

doc.book.html:
	$(ECHO) "Creating the complete HTML documentation book ..."
	$(LATEX2HTML) \
		-local_icons \
		-split 5 \
		-toc_depth 4 \
		-toc_stars \
		-t "VMD Documentation" \
		-address "vmd@ks.uiuc.edu" \
		-up_url "http://www.ks.uiuc.edu/Research/vmd/" \
		-up_title "Main VMD Page" \
		$(MAIN_ALL_DOC)

config.h vmd_version.tex:
	true


clean:
	$(ECHO) "Erasing all generated files and directories..."
	$(ECHO) " + Erasing Latex generated files."
	$(DELETE) *.aux *.dvi *.idx *.lof *.log *.lot *.pdf *.ps *.toc ?g.out
	$(DELETE) *.ilg *.ind
	$(TOUCH)  ug.ind 
	$(ECHO) " + Erasing HTML tar files and directories."
	$(DELETE) *.html.tar.gz *.html.tar
	$(ECHO) "   - Erasing ig HTML directory."
	$(DELETEDIR) ig 
	$(ECHO) "   - Erasing pg HTML directory."
	$(DELETEDIR) pg
	$(ECHO) "   - Erasing ug HTML directory."
	$(DELETEDIR) ug
	$(ECHO) " + Erasing doxygen auto-generated source documentation."
	$(ECHO) "   - Erasing dox HTML/Latex directories."
	$(DELETEDIR) dox

