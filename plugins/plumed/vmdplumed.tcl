# VMD Plumed tool  - a GUI to compute collective variables
# over a trajectory
#
#     Copyright (C) 2012  Universitat Pompeu Fabra and
#                         National Research Council of Italy
#
#     Author              Toni Giorgino  (toni.giorgino@gmail.com)
#
#     This program is available under either the 3-clause BSD license,
#     (e.g. see http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html)
#
# $Id: vmdplumed.tcl,v 1.4 2012/07/19 15:24:55 johns Exp $ 


package provide plumed 1.0

# vmd_install_extension plumed plumed_tk "Analysis/Collective variable analysis (PLUMED)"

namespace eval Plumed:: {
    namespace export plumed
    variable debug 0
    variable plugin_version "1.0"
    variable w                                          ;# handle to main window
    variable textfile "META_INP"
    variable plugin_name "PLUMED collective variable analysis tool"
    variable driver_path "(Executable not found in path. Please see manual to install, or click 'Browse...' to locate it.)"
    variable text_instructions \
"Enter collective variable definitions below, in PLUMED syntax.  
Click 'Plot' to evaluate them on the 'top' trajectory.  
VMD atom selections in square brackets expand automatically.
For example:

    protein-> \[chain A and name CA\] protein<-
    ligand->  \[chain B and noh\]     ligand<-

    DISTANCE LIST <protein> <ligand>
    COORD LIST <protein> <ligand>  NN 6 MM 12 D_0 5.0 R_0 0.5 
"
    variable empty_meta_inp "\nDISTANCE LIST 1 200      ! Just an example\n"
}

proc plumed_tk {} {
    Plumed::plumed
    return $Plumed::w
}


proc Plumed::plumed {} { 
    variable w
    variable textfile
    variable plugin_name
    variable driver_path
    variable cv_n 1
    variable pbc_type 1
    variable pbc_boxx
    variable pbc_boxy
    variable pbc_boxz
    variable plot_wall 0
    variable plot_points 0
    variable empty_meta_inp
    variable text_instructions

    # If already initialized, just turn on
    if { [winfo exists .textview] } {
	wm deiconify $w
	return
    }

    variable file_types {
	{"Plumed Files" { .plumed .metainp .meta_inp } }
	{"Text Files" { .txt .TXT} }
	{"All Files" * }
    }

    set w [toplevel ".plumed"]
    wm title $w "$plugin_name"
#    wm resizable $w 0 0

    # If driver is in the path, make it default
    set dr [ auto_execok driver ]
    if { $dr != "" } { set driver_path $dr }

    # If PBC exist, use them
    catch { molinfo top get {a b c} } vmdcell
    if { [llength $vmdcell] == 3 } {
	lassign $vmdcell a b c 
	if { [expr $a * $b * $c ] > 1.0 } {
	    lassign $vmdcell pbc_boxx pbc_boxy pbc_boxz
	}
    }


    frame $w.menubar -relief raised -bd 2 ;# frame for menubar
    pack $w.menubar -padx 1 -fill x

    ## file menu
    menubutton $w.menubar.file -text File -underline 0 -menu $w.menubar.file.menu
    menu $w.menubar.file.menu -tearoff no
    $w.menubar.file.menu add command -label "New" -command  Plumed::file_new
    $w.menubar.file.menu add command -label "Open..." -command Plumed::file_open
    $w.menubar.file.menu add command -label "Save" -command  Plumed::file_save -acce Ctrl-S
    $w.menubar.file.menu add command -label "Save as..." -command  Plumed::file_saveas
    $w.menubar.file.menu add command -label "Export..." -command  Plumed::file_export
    $w.menubar.file.menu add separator
    $w.menubar.file.menu add command -label "Batch analysis..." -command  Plumed::batch_gui
    $w.menubar.file.menu add separator
    $w.menubar.file.menu add command -label "Quit" -command  Plumed::file_quit
    $w.menubar.file config -width 5
    bind $w <Control-s> Plumed::file_save

    ## edit
    menubutton $w.menubar.edit -text Edit -underline 0 -menu $w.menubar.edit.menu
    menu $w.menubar.edit.menu -tearoff no
    $w.menubar.edit.menu add command -label "Undo" -command  "$::Plumed::w.txt.text edit undo" -acce Ctrl-Z
    $w.menubar.edit.menu add command -label "Redo" -command  "$::Plumed::w.txt.text edit redo"
    $w.menubar.edit.menu add separator
    $w.menubar.edit.menu add command -label "Cut" -command  "tk_textCut $::Plumed::w.txt.text" -acce Ctrl-X
    $w.menubar.edit.menu add command -label "Copy" -command  "tk_textCopy $::Plumed::w.txt.text" -acce Ctrl-C
    $w.menubar.edit.menu add command -label "Paste" -command  "tk_textPaste $::Plumed::w.txt.text" -acce Ctrl-V
    $w.menubar.edit config -width 5

    ## Templates
    set templates {  
	"Group definition" {groupname-> \[chain A\] groupname<-}
	"-" "-"
	"Absolute position" "POSITION LIST <xx>   DIR XYZ"
	"Distance" "DISTANCE LIST <xx> <yy>    DIR XYZ"
	"Minimum distance" "MINDIST LIST <xx> <yy>    BETA 500"
	"Angle" "ANGLE LIST <xx> <oo> <yy>"
	"Dihedral" "TORSION LIST <xx> <yy> <zz> <ww>"
	"Contacts" "COORD LIST <xx> <yy>  NN 6 MM 12 D_0 2.5 R_0 0.5 "
	"Hydrogen bonds" "HBONDS LIST <xx> <yy>  TYPE nn"
	"Interfacial water" "WATERBRIDGE LIST <type1> <type2> <solvent> NN 8 MM 12 R_0 4.0"
	"Radius of gyration" "RGYR LIST <xx>"
	"Trace of the inertia tensor" "INERTIA LIST <xx>"
	"Electrostatic potential" "ELSTPOT LIST <xx> <yy>  R_0 4.0   CUT 12.0"
	"Electric dipole" "DIPOLE LIST <xx>"
	"-" "-"
	"RMSD from a reference" "TARGETED TYPE RMSD FRAMESET /tmp/reference.pdb SQRT  ! use absolute pathname"
	"S-path variable" "S_PATH TYPE RMSD FRAMESET /tmp/frame_ NFRAMES 2 LAMBDA 9.0  ! use absolute pathnames"
	"Z-path variable" "Z_PATH TYPE RMSD FRAMESET /tmp/frame_ NFRAMES 2 LAMBDA 9.0  ! use absolute pathnames"
	"-" "-"
	"Dihedral correlation" "DIHCOR NDIH 3\nA1 B1 C1 D1\nA2 B2 C2 D2\nA3 B3 C3 D3\n"
	"Alpha-beta similarity" "ALPHABETA NDIH 3\nA1 B1 C1 D1 ref1\nA2 B2 C2 D2 ref2\nA3 B3 C3 D3 ref3\n"
	"Alpha RMSD" "ALPHARMSD LIST <xx>  R_0 0.8  NN 8  MM 12  NOPBC"
	"Antiparallel beta RMSD" "ANTIBETARMSD LIST <xx>  R_0 0.8  NN 8  MM 12  STRANDS_CUTOFF 10  NOPBC"
	"Parallel beta RMSD" "PARABETARMSD LIST <xx>  R_0 0.8  NN 8  MM 12  STRANDS_CUTOFF 10  NOPBC"
	"Torsional RMSD" "RMSDTOR NDIH 3\nA1 B1 C1 D1 ref1\nA2 B2 C2 D2 ref2\nA3 B3 C3 D3 ref3\n"
	"Puckering coordinates" "PUCKERING LIST <xx>  TYPE  PHI|THETA|Q"
	"-" "-"
	"Upper wall (allow lower)" "UWALL CV nn LIMIT xx   KAPPA 100.0 EXP 4.0 EPS 1.0 OFF 0.0"
	"Lower wall (allow higher)" "LWALL CV nn LIMIT xx   KAPPA 100.0 EXP 4.0 EPS 1.0 OFF 0.0"
    }

    menubutton $w.menubar.insert -text "Templates" -underline 0 -menu $w.menubar.insert.menu
    menu $w.menubar.insert.menu -tearoff yes
    foreach { disp insr } $templates {
	if {$disp == "-" } {
	    $w.menubar.insert.menu add separator
	} else {
	    $w.menubar.insert.menu add command -label $disp -command \
		"$::Plumed::w.txt.text insert insert \"$insr\n\""
	}
    }
    # In preparation for V2
    set [namespace current]::lists_use_comma 0
#    $w.menubar.insert.menu add separator
#    $w.menubar.insert.menu add checkbutton -label {Use comma to expand lists} \
#	-variable [namespace current]::lists_use_comma 
    bind $w <Control-g> "$::Plumed::w.menubar.insert.menu invoke 1"
    $w.menubar.insert.menu entryconfigure 1 -accelerator Ctrl-G
    $w.menubar.insert config -width 10

    ## Structural
    menubutton $w.menubar.structure -text Structure -underline 0 -menu $w.menubar.structure.menu
    menu $w.menubar.structure.menu -tearoff no
    $w.menubar.structure.menu add command -label "Build reference structure..." -command Plumed::reference_gui
    $w.menubar.structure.menu add command -label "Insert native contacts CV..." -command Plumed::nc_gui
    $w.menubar.structure.menu add command -label "Insert Ramachandran \u03c6/\u03c8/\u03c9 CVs..." -command Plumed::rama_gui
    $w.menubar.structure config -width 8



    ## help menu
    menubutton $w.menubar.help -text Help -underline 0 -menu $w.menubar.help.menu
    menu $w.menubar.help.menu -tearoff no
    $w.menubar.help.menu add command -label "On the $plugin_name" -command \
	"vmd_open_url http://www.multiscalelab.org/toni/PlumedCVTool"
    $w.menubar.help.menu add command -label "How to install the 'driver' binary" -command \
	"vmd_open_url http://www.multiscalelab.org/toni/PlumedCVTool"
    $w.menubar.help.menu add command -label "On PLUMED" -command \
	"vmd_open_url http://www.plumed-code.org"
    $w.menubar.help.menu add command -label "PLUMED user's guide and CV syntax" -command \
	"vmd_open_url http://www.plumed-code.org/documentation/manual_1-2-2.pdf?attredirects=0&d=1"
#	"vmd_open_url http://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxwbHVtZWR3ZWJ8Z3g6MmU2NDQ0Mzc0OTAxNDMzNg"
    $w.menubar.help.menu add command -label "About..."  -command [namespace current]::help_about
    # XXX - set menubutton width to avoid truncation in OS X
    $w.menubar.help config -width 5


    pack $w.menubar.file -side left
    pack $w.menubar.edit -side left
    pack $w.menubar.insert -side left
    pack $w.menubar.structure -side left
    pack $w.menubar.help -side right


    ##
    ## main window area
    ## 
    frame $w.txt
    label $w.txt.label  -textvariable Plumed::textfile
    text $w.txt.text -wrap none -undo 1 -autoseparators 1 -bg #ffffff -bd 2 \
	-yscrollcommand "$::Plumed::w.txt.vscr set" -font {courier 12}
    label $w.explanation -text $text_instructions -justify left -bd 2 -relief solid -padx 10 -pady 10
    file_new
    $w.txt.text window create 1.0 -window $w.explanation -padx 100 -pady 10
    scrollbar $w.txt.vscr -command "$::Plumed::w.txt.text yview"
    pack $w.txt.label -fill x 
    pack $w.txt.text  -side left -fill both -expand 1
    pack $w.txt.vscr  -side left -fill y -expand 0
    pack $w.txt  -fill both -expand 1


    ##
    ## main window area
    ## 
    pack [ labelframe $w.options -relief ridge -bd 2 -text "Options" -padx 1m -pady 1m ] \
	-side top -fill x

    pack [ frame $w.options.pbc   ]  -side top -fill x
    pack [ label $w.options.pbc.cvtext -text "Number of collective variables: " ] -side left
    pack [ entry $w.options.pbc.cv -width 6 -textvariable  [namespace current]::cv_n ] -side left
    pack [ label $w.options.pbc.spacer1 -text " " ] -side left -expand true -fill x
    pack [ radiobutton $w.options.pbc.pbcno -value 1 -text "No PBC" -variable [namespace current]::pbc_type ] -side left
    pack [ radiobutton $w.options.pbc.pbcdcd -value 2 -text "From DCD" -variable [namespace current]::pbc_type ] -side left
    pack [ radiobutton $w.options.pbc.pbcbox -value 3 -text "Box:" -variable [namespace current]::pbc_type ] -side left
    pack [ entry $w.options.pbc.boxx -width 6 -textvariable  [namespace current]::pbc_boxx ] -side left
    pack [ entry $w.options.pbc.boxy -width 6 -textvariable  [namespace current]::pbc_boxy ] -side left
    pack [ entry $w.options.pbc.boxz -width 6 -textvariable  [namespace current]::pbc_boxz ] -side left
    pack [ label $w.options.pbc.spacer2 -text " " ] -side left -expand true -fill x
    pack [ checkbutton $w.options.pbc.plotwall -text "Plot wall potential" -variable  [namespace current]::plot_wall ] -side left
    pack [ checkbutton $w.options.pbc.inspector -text "Points" -variable  [namespace current]::plot_points ] -side left


    pack [ frame $w.options.location ] -side top -fill x
    pack [ label $w.options.location.text -text "Location of PLUMED's driver binary: " ] -side left -expand 0
    pack [ entry $w.options.location.path -width 40 -textvariable \
	       [namespace current]::driver_path ] -side left -expand 1 -fill x
    pack [ button $w.options.location.browse -text "Browse..." -relief raised \
	   -command [namespace current]::location_browse   ] -side left -expand 0

    pack [ frame $w.plot ] -side top -fill x
    pack [ button $w.plot.plot -text "Plot" -relief raised -command [namespace current]::do_compute \
	       -pady 2 -bd 2 ] -side left -fill x -expand 1
}


# ==================================================


proc Plumed::file_new { } {
    variable w
    variable textfile
    variable empty_meta_inp

    $w.txt.text delete 1.0 {end - 1c}
    set textfile "untitled.plumed"
    $w.txt.text insert end $empty_meta_inp
}


proc Plumed::file_open { } {
    variable w
    variable textfile
    variable file_types

    set textfile [tk_getOpenFile -filetypes $file_types \
		      -initialfile "$Plumed::textfile"    ]

    if { $textfile == "" } { return }
    set rc [ catch { set fd [open $textfile "r"] } ]
    if { $rc == 1} { return }
    close $fd

    $w.txt.text delete 1.0 {end - 1c}
    $w.txt.text insert end [read_file $textfile ]
}

proc Plumed::file_save { } {
    variable w
    variable textfile

    if { $textfile == "untitled.plumed" } {
	Plumed::file_saveas
	return
    }

    set rc [ catch { set fd [open $textfile "w"] } ]
    if { $rc == 1} {
	tk_messageBox -title "Error" -parent .plumed -message "Failed to open file $textfile"
	return
    } else {
	puts $fd  [$w.txt.text get 1.0 {end -1c}] 
	close $fd
    }
}

proc Plumed::file_saveas { } {
    variable w
    variable file_types
    variable textfile

    set textfile [tk_getSaveFile -filetypes $file_types \
		       -initialfile $Plumed::textfile ]
    set rc [ catch { set fd [open $textfile "w"] } ]
    if { $rc == 1} {
	puts "failed to open file $textfile"
	return
    }
    puts $fd [$w.txt.text get 1.0 {end -1c}]
    close $fd
}

proc Plumed::file_export { } {
    variable w
    set file_types {
	{"All Files" * }
    }
    set textfile [tk_getSaveFile -filetypes $file_types \
		       -initialfile "META_INP"       ]
    set rc [ catch { set fd [open $textfile "w"] } ]
    if { $rc == 1} {
	puts "failed to open file $textfile"
	return
    }
    puts $fd  [ Plumed::replace_serials [$w.txt.text get 1.0 {end -1c}] ]
    puts $fd  "ENDMETA"
    close $fd
}


# Well, not really quit
proc Plumed::file_quit { } {
    variable w
    wm withdraw $w
}

# Browse for executable
proc Plumed::location_browse { } {
    variable driver_path
    set tmp [ tk_getOpenFile  ]
    if { $tmp != "" } {
	set driver_path $tmp
    }
}


proc Plumed::help_about { {parent .plumed} } {
    variable plugin_name
    variable plugin_version

    tk_messageBox -title "About" -parent $parent -message \
"
$plugin_name 
Version $plugin_version

Toni Giorgino <toni.giorgino@gmail.com>

Copyright (c) 2010-2012:

Institute for Biomedical Engineering (ISIB-CNR),
National Research Council of Italy

and

Computational Biophysics Group (GRIB-IMIM-UPF),
Universitat Pompeu Fabra, Barcelona


\$Id: vmdplumed.tcl,v 1.4 2012/07/19 15:24:55 johns Exp $
"
}



# ==================================================                                                 

# http://wiki.tcl.tk/772
proc Plumed::tmpdir { } {
    global tcl_platform
    switch $tcl_platform(platform) {
	unix {
	    set tmpdir /tmp   ;# or even $::env(TMPDIR), at times.
	} macintosh {
	    set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
	} default {
	    set tmpdir [pwd]
	    catch {set tmpdir $::env(TMP)}
	    catch {set tmpdir $::env(TEMP)}
	}
    }
    return $tmpdir
}

# Write a PDB file using charges and masses in the topology
# (required by PLUMED's "driver" utility). 
proc Plumed::writePlumed { sel filename } {
    set old [ $sel get { x y z resid occupancy beta } ]

    $sel set x 0;		# workaround for xyz>1000 in PDB
    $sel set y 0
    $sel set z 0
    $sel set resid 1;	# workaround for PLUMED not reading ascii RESID
    $sel set occupancy [ $sel get mass ]	
    $sel set beta [ $sel get charge ]
    # FIXME: plumed does not read serial > 100k
    $sel writepdb $filename

    $sel set { x y z resid occupancy beta } $old
}


proc Plumed::replace_serials { intxt }  {
    variable lists_use_comma
    set re {\[(.+?)\]}
    set lorig {}
    set lnew {}
    set lcount {}
    while { [ regexp $re $intxt junk orig ] } {
	lappend lorig $orig
	set as [ atomselect top $orig ]
	set new [ $as get serial ]
	$as delete
	if $lists_use_comma {
	    set new [string map { " " , } $new ]
	}
	lappend lnew $new
	lappend lcount [llength $new]
	regsub $re $intxt $new intxt
    }
    set out $intxt
    set out "$out
!
! The above script includes the following replacements, based on 
! a structure named [molinfo top get name] with [molinfo top get numatoms] atoms.\n!\n"
    foreach orig $lorig new $lnew cnt $lcount {
	set out "${out}! $orig -> (list of $cnt atoms)\n"
    }
    return $out
}


proc Plumed::write_meta_inp { meta } {
    variable w
    set text [$w.txt.text get 1.0 {end -1c}]
    set fd [open $meta w]
    puts $fd [ Plumed::replace_serials $text ]
    puts $fd "PRINT W_STRIDE 1  ! line added by vmdplumed for visualization"
    puts $fd "ENDMETA"
    close $fd
}

proc Plumed::get_pbc { } {
    variable pbc_type
    variable pbc_boxx
    variable pbc_boxy
    variable pbc_boxz
    set pbc [ switch $pbc_type {
	1 {format "-nopbc"}
	2 {format "" }
	3 {format "-cell $pbc_boxx $pbc_boxy $pbc_boxz" } } ]
    return $pbc
}

proc Plumed::transpose matrix {
    set cmd list
    set i -1
    foreach col [lindex $matrix 0] {append cmd " \$[incr i]"}
    foreach row $matrix {
        set i -1
        foreach col $row {lappend [incr i] $col}
    }
    eval $cmd
}


# Return the contents of a file
proc Plumed::read_file { fname } {
    set fd [open $fname r]
    set dtext [read $fd]
    close $fd
    return $dtext
}

# from rmsdtt
proc Plumed::index2rgb {i} {
  set len 2
  lassign [colorinfo rgb $i] r g b
  set r [expr int($r*255)]
  set g [expr int($g*255)]
  set b [expr int($b*255)]
  #puts "$i      $r $g $b"
  return [format "#%.${len}X%.${len}X%.${len}X" $r $g $b]
}

proc Plumed::dputs { text } {
    variable debug
    if {$debug} {
	puts "DEBUG $text"
    }
}


# ==================================================

# TONI TOOLTIPS  From http://wiki.tcl.tk/1954
proc Plumed::setTooltip {widget text} {
        if { $text != "" } {
                # 2) Adjusted timings and added key and button bindings. These seem to
                # make artifacts tolerably rare.
                bind $widget <Any-Enter>    [list after 500 [list Plumed::showTooltip %W $text]]
                bind $widget <Any-Leave>    [list after 500 [list destroy %W.tooltip]]
                bind $widget <Any-KeyPress> [list after 500 [list destroy %W.tooltip]]
                bind $widget <Any-Button>   [list after 500 [list destroy %W.tooltip]]
        }
 }

proc Plumed::showTooltip {widget text} {
        global tcl_platform
        if { [string match $widget* [winfo containing  [winfo pointerx .] [winfo pointery .]] ] == 0  } {
                return
        }

        catch { destroy $widget.tooltip }

        set scrh [winfo screenheight $widget]    ; # 1) flashing window fix
        set scrw [winfo screenwidth $widget]     ; # 1) flashing window fix
        set tooltip [toplevel $widget.tooltip -bd 1 -bg black]
        wm geometry $tooltip +$scrh+$scrw        ; # 1) flashing window fix
        wm overrideredirect $tooltip 1

        if {$tcl_platform(platform) == {windows}} { ; # 3) wm attributes...
                wm attributes $tooltip -topmost 1   ; # 3) assumes...
        }                                           ; # 3) Windows
        pack [label $tooltip.label -bg lightyellow -fg black -text $text -justify left]

        set width [winfo reqwidth $tooltip.label]
        set height [winfo reqheight $tooltip.label]

        set pointer_below_midline [expr [winfo pointery .] > [expr [winfo screenheight .] / 2.0]]                ; # b.) Is the pointer in the bottom half of the screen?

        set positionX [expr [winfo pointerx .] - round($width / 2.0)]    ; # c.) Tooltip is centred horizontally on pointer.
        set positionY [expr [winfo pointery .] + 35 * ($pointer_below_midline * -2 + 1) - round($height / 2.0)]  ; # b.) Tooltip is displayed above or below depending on pointer Y position.

        # a.) Ad-hockery: Set positionX so the entire tooltip widget will be displayed.
        # c.) Simplified slightly and modified to handle horizontally-centred tooltips and the left screen edge.
        if  {[expr $positionX + $width] > [winfo screenwidth .]} {
                set positionX [expr [winfo screenwidth .] - $width]
        } elseif {$positionX < 0} {
                set positionX 0
        }

        wm geometry $tooltip [join  "$width x $height + $positionX + $positionY" {}]
        raise $tooltip

        # 2) Kludge: defeat rare artifact by passing mouse over a tooltip to destroy it.
        bind $widget.tooltip <Any-Enter> {destroy %W}
        bind $widget.tooltip <Any-Leave> {destroy %W}
 }


# ==================================================                                                 

proc Plumed::rama_gui { } {
    variable rama_sel "protein"
    variable rama_phi 1
    variable rama_psi 1
    variable rama_omega 0

    toplevel .plumedrama -bd 4
    wm title .plumedrama "Insert Ramachandran CVs"
    pack [ label .plumedrama.head1 -text "Insert CVs for the Ramachandran angles of the matched residues.
N-CA-C atom naming is assumed for backbone atoms.
Dihedrals involving atoms outside the selection are skipped.
" ] -side top -fill x 

    pack [ frame .plumedrama.sel ] -side top -fill x
    pack [ label .plumedrama.sel.txt -text "Selection: backbone and " ] -side left -fill x
    pack [ entry .plumedrama.sel.in -width 20 -textvariable [namespace current]::rama_sel ] -side left -expand 1 -fill x

    pack [ frame .plumedrama.cv ] -side top -fill x
    pack [ label .plumedrama.cv.txt -text "Dihedral angles: " ] -side left -fill x
    pack [ checkbutton .plumedrama.cv.phi -text Phi -variable  [namespace current]::rama_phi ] -side left
    pack [ checkbutton .plumedrama.cv.psi -text Psi -variable  [namespace current]::rama_psi ] -side left
    pack [ checkbutton .plumedrama.cv.omega -text Omega -variable  [namespace current]::rama_omega ] -side left

    pack [ frame .plumedrama.act ] -side top -fill x
    pack [ button .plumedrama.act.ok -text "Insert" -relief raised -command \
	       { Plumed::rama_insert } ] -side left -fill x -expand 1
    pack [ button .plumedrama.act.close -text "Close" -relief raised \
	       -command {  destroy .plumedrama }   ] -side left -fill x -expand 1

}

# Insert Ramachandran angles. Uses the following atoms
#       Cm   N CA C   Np CAp
# phi    +   +  + +
# psi        +  + +   +
# omega         + +   +  +
# where "m/p" means previous/next residue Computation is done on the
# base of "residue" (unique), but they will be printed as "resid"
# (human-readable)

proc Plumed::rama_insert {} {
    variable rama_sel
    variable rama_phi
    variable rama_psi
    variable rama_omega
    variable w

    set nnew 0

    set sel [atomselect top "($rama_sel) and name CA and backbone"]
    set rlist [lsort -integer -uniq [$sel get residue]]
    $sel delete

    if {[llength $rlist] == 0} {
	tk_messageBox -title "Error" -parent .plumedrama -message "Selection is empty."
	return	
    }

    foreach r $rlist {
	if {$r == 0} {
	    set rm1 1000000;	# non-existent residue, kludge not to atomselect a negative one
	} else {
	    set rm1 [expr $r-1]
	}
	set rp1 [expr $r+1]
	set Cm  [atomselect top "($rama_sel) and backbone and residue $rm1 and name C"]
	set N   [atomselect top "($rama_sel) and backbone and residue $r and name N"]
	set CA  [atomselect top "($rama_sel) and backbone and residue $r and name CA"]
	set C   [atomselect top "($rama_sel) and backbone and residue $r and name C"]
	set Np  [atomselect top "($rama_sel) and backbone and residue $rp1 and name N"]
	set CAp [atomselect top "($rama_sel) and backbone and residue $rp1 and name CA"]
	set rid [format "%s%d" [$CA get resname] [$CA get resid]]; # human readable
	if {$rama_phi && [rama_insert_cv_maybe $Cm $N $CA $C             PHI $rid] } { 	incr nnew }
	if {$rama_psi && [rama_insert_cv_maybe     $N $CA $C  $Np        PSI $rid] } {	incr nnew }
	if {$rama_omega && [rama_insert_cv_maybe      $CA $C  $Np $CAp  OMEGA $rid] } { incr nnew }
	$Cm delete; 
	$N delete; $CA delete; $C delete
	$Np delete; $CAp delete
    }
    $w.txt.text insert insert "! The above list contains $nnew Ramachandan CVs\n" 
}


# Return the line computing a torsion CV defined by the arguments iff all of them are valid
proc Plumed::rama_insert_cv_maybe {A B C D angle rid} {
    variable w
    set cv_line "TORSION LIST %d %d %d %d  ! %s_%s\n"
    set oos_msg "! No dihedral %s for residue %s: out of selection\n"
    set topo_msg "! No dihedral %s for residue %s: chain break\n"
    if { [$A num]==0 || [$B num]==0 || [$C num]==0 || [$D num]==0 } {
	set r [format $oos_msg $angle $rid]
	set ok 0
    } elseif { [llength [lsort -uniq -integer [list \
  	         [$A get fragment] [$B get fragment] [$C get fragment] [$D get fragment] ]]] != 1 } {
	# above; check that all fragment IDs are equal (check enforced on topology)
	# could also be segid or chain (check enforced on "logical" structure)
	set r [format $topo_msg $angle $rid]
	set ok 0
    } else {
	set r [format $cv_line  \
		   [$A get serial] [$B get serial] [$C get serial] [$D get serial]  \
		   $rid $angle ]
	set ok 1
    } 
    $w.txt.text insert insert $r
    return $ok
}





# ==================================================                                                 

proc Plumed::reference_gui { } {
    variable refalign "backbone"
    variable refmeas "name CA"
    variable reffile "reference.pdb"
    variable refmol top
    toplevel .plumedref -bd 4
    wm title .plumedref "Build reference structure"
    pack [ label .plumedref.title -text "Convert top molecule's current frame\ninto a reference file for FRAMESET-type analysis:" ] -side top
    pack [ frame .plumedref.align ] -side top -fill x
    pack [ label .plumedref.align.aligntext -text "Alignment set: " ] -side left
    pack [ entry .plumedref.align.align -width 20 -textvariable [namespace current]::refalign ] -side left -expand 1 -fill x
    pack [ frame .plumedref.meas ] -side top -fill x
    pack [ label .plumedref.meas.meastext -text "Displacement set: " ] -side left
    pack [ entry .plumedref.meas.meas -width 20 -textvariable [namespace current]::refmeas ] -side left -expand 1 -fill x
    pack [ frame .plumedref.mol ] -side top -fill x
    pack [ label .plumedref.mol.moltext -text "Numbering for molecule: " ] -side left
    pack [ entry .plumedref.mol.mol -width 20 -textvariable [namespace current]::refmol ] -side left -expand 1 -fill x
    pack [ frame .plumedref.file ] -side top -fill x
    pack [ label .plumedref.file.filetxt -text "File to write: " ] -side left
    pack [ entry .plumedref.file.file -width 20 -textvariable [namespace current]::reffile ] -side left -expand 1 -fill x
    pack [ button .plumedref.file.filebrowse -text "Browse..." -relief raised \
	       -command { Plumed::reference_set_reffile [ tk_getSaveFile  -initialfile "$::Plumed::reffile" ] }   ] -side left -expand 0
    pack [ frame .plumedref.act ] -side top -fill x
    pack [ button .plumedref.act.ok -text "Write" -relief raised -command \
	       { Plumed::reference_write } ] -side left -fill x -expand 1
    pack [ button .plumedref.act.cancel -text "Close" -relief raised \
	       -command {  destroy .plumedref }   ] -side left -fill x -expand 1
}

proc Plumed::reference_set_reffile { x } { variable reffile; if { $x != "" } {set reffile $x} }; # why??

proc Plumed::reference_write { } {
    variable refalign
    variable refmeas
    variable reffile
    variable refmol

    set tmpf [ file join [ Plumed::tmpdir ] "reftmp.[pid].pdb" ]

    # From where new serials are taken
    set asnew [ atomselect $refmol "($refalign) or ($refmeas)" ]
    set newserial [ $asnew  get serial ]
    $asnew delete

    set asall [ atomselect top all]
    set asref [ atomselect top "($refalign) or ($refmeas)" ]
    set asalign [ atomselect top $refalign ] 
    set asmeas  [ atomselect top $refmeas ] 
    
    set oldserial [ $asref  get serial ]
    if { [llength $oldserial] != [llength $newserial] } {
	tk_messageBox -title "Error" -parent .plumedref -message "Selection ($refalign) or ($refmeas) has different number of atoms in molecule $refmol ([llength $newserial]) versus top ([llength $oldserial])."
	return
    }

    set old [ $asall get {occupancy beta segname} ]; # backup

    $asall set occupancy 0
    $asall set beta 0
    $asall set segname XXXX

    $asalign set occupancy 1
    $asmeas  set beta 1
    $asref   set segname YYYY

    $asall writepdb $tmpf

    $asall set {occupancy beta segname} $old; # restore
    $asall delete
    $asref delete
    $asalign delete
    $asmeas delete


    # i.e. grep YYYY $tmpd/reftmp.pdb > $reffile
    set fdr [ open $tmpf r ]
    set fdw [ open $reffile w ]
    set i 0
    while { [gets $fdr line] != -1 } {
	if { [ regexp {YYYY} $line ] } {
	    # workaround plumed bug in PDB reader
	    set line [string replace $line 21 21 " "]
	    # replace serial
	    set line [string replace $line 6 10 [ format "%5s" [ lindex $newserial $i ] ] ]
	    puts $fdw $line
	    incr i
	}
    }
    close $fdr
    close $fdw
    file delete $tmpf
    puts "Done."
}


# ==================================================

proc Plumed::batch_gui {} {
    variable batch_dir
    variable batch_nthreads 1
    variable batch_extension colvar
    variable merge_results 0
    variable batch_keep 0
    variable batch_progress_width 400.0

    if { ! [ info exists batch_dir ] } { set batch_dir [pwd] }

    global tcl_platform
    if { $tcl_platform(platform) != "unix" } {
	tk_messageBox -title "Error" -parent .plumed -message "Sorry, unix only"; return; }

    toplevel .plumedbatch -bd 4
    wm title .plumedbatch "Batch analysis"
    pack [ label .plumedbatch.info -text "Process all DCD files in directory..." \
	        ] -side top -fill x 
    pack [ frame .plumedbatch.dir ] -side top -fill x
    pack [ entry .plumedbatch.dir.dir -width 60 -textvariable [namespace current]::batch_dir ] -side left -expand 1 -fill x
    pack [ button .plumedbatch.dir.dirbrowse -text "Browse..." -relief raised \
	       -command { Plumed::batch_setdir [ tk_chooseDirectory -mustexist 1 -parent .plumedbatch ] }   ] -side left -expand 0

    pack [ frame .plumedbatch.merge ] -side top -fill x
    pack [ checkbutton .plumedbatch.merge.button -text "Combine results"  -variable \
	       [namespace current]::merge_results -command Plumed::batch_update_gui ] -side left -fill x
    pack [ label .plumedbatch.merge.info  -relief sunken -bd 1 -text " " -pady 3] -side left   -fill x -expand 1

    pack [ frame .plumedbatch.threads ] -side top -fill x
    pack [ checkbutton .plumedbatch.threads.keep -text "Keep existing"  \
	   -variable  [namespace current]::batch_keep ] -side left -fill x 
    pack [ label .plumedbatch.threads.info2 -text "   File extension: " ] -side left -fill x
    pack [ entry .plumedbatch.threads.extension -width 15 -textvariable [namespace current]::batch_extension \
	       -validate key -validatecommand [list Plumed::batch_update_gui %P] ] \
	-side left -expand 1 -fill x
    pack [ label .plumedbatch.threads.info -text "Concurrent processes: " ] -side left -fill x
    pack [ entry .plumedbatch.threads.nthreads -width 5 -textvariable [namespace current]::batch_nthreads ] \
	-side left -expand 1 -fill x
    catch { set batch_nthreads [exec grep processor /proc/cpuinfo | wc -l] }

    pack [ canvas .plumedbatch.progress -relief sunken -bd 1 -width $batch_progress_width -height 20  ]  \
	-side top  -padx 10 -pady 10 -expand 1
    .plumedbatch.progress create rectangle 0 0 0 20 -tags bar -fill lavender
    .plumedbatch.progress create text [ expr $batch_progress_width / 2 ] 10 -tags text 

    pack [ frame .plumedbatch.act ] -side top -fill x
    pack [ button .plumedbatch.act.ok -text "Start" -relief raised -command \
	       { Plumed::batch_start } ] -side left -fill x -expand 1 
    pack [ button .plumedbatch.act.abort -text "Abort" -relief raised -command \
	       { set Plumed::batch_abort 1; Plumed::batch_abort_do } ] -side left -fill x -expand 1 
    pack [ button .plumedbatch.act.cancel -text "Close" -relief raised \
	       -command {  destroy .plumedbatch }   ] -side left -fill x -expand 1
    Plumed::batch_update_gui
}

proc Plumed::batch_setdir { x } { 
    variable batch_dir; if {$x != "" } {
	set batch_dir $x
	.plumedbatch.progress coords bar 0 0 0 20 
	.plumedbatch.progress itemconfigure text -text "Counting files..."
	update
	.plumedbatch.progress itemconfigure text -text [ format "%d files to process" \
	     [ llength [ glob -directory $batch_dir *.dcd ] ] ]
    } 
}

proc Plumed::batch_update_gui { { ext __UNSET__ } } {
    variable merge_results
    variable batch_extension
    if { $ext == "__UNSET__" } { 	# not called from validate
	set ext $batch_extension 
    }
    .plumedbatch.merge.info config -text \
	[ switch -- $merge_results  \
	      0 { format "Results will be stored in separate files: *.dcd -> *.dcd.$ext" } \
	      1 { format "A single '$ext.all' file will be created in the directory" } ]
    return 1
}

proc Plumed::batch_progress { x { y 100 } } { 
    variable batch_progress_width
    .plumedbatch.progress coords bar 0 0 [ expr {int($batch_progress_width*$x/$y) } ] 20 
    .plumedbatch.progress itemconfigure text -text [ format "Processed %d of %d" $x $y ]
}


# somewhat contrived implementation of a process pool
proc Plumed::batch_start { } {
    variable merge_results
    variable batch_nthreads 
    variable batch_ongoing {}
    variable batch_abort 0
    variable batch_ok 0
    variable batch_nfiles
    variable batch_extension
    variable batch_dir;		# must be absolute
    variable batch_keep
    variable driver_path
    variable cv_n

    .plumedbatch.progress itemconfigure text -text "Please wait..."

    set tmpd "[ Plumed::tmpdir ]/vmdplumed.[pid]"
    file mkdir $tmpd
    set owd [ pwd ]
    cd $tmpd

    set pdb temp.pdb
    writePlumed [atomselect top all] $pdb

    set meta META_INP
    write_meta_inp $meta

    set pbc [ get_pbc ]

    set filelist [ lsort -dictionary [ glob -directory $batch_dir *.dcd ] ]
    set batch_nfiles [ llength $filelist ]
    
    set fi 0;			# file index
    while { [llength $filelist] > 0  && !$batch_abort } {
	dputs "Batch ongoing is now $batch_ongoing"
	if { [llength $batch_ongoing] < $batch_nthreads } {
	    # spawn a new process, pop it from filelist and put it in batch_ongoing (kindof)
	    set fn [ lindex $filelist 0 ]
	    set filelist [ lreplace $filelist 0 0 ]

	    # if file exists and user requested not overwrite, skip
	    if { $batch_keep==1 && [file exists $fn.$batch_extension ] } {
		dputs "File  $fn.$batch_extension exists, skipping"
		incr batch_ok
		batch_progress $batch_ok $batch_nfiles
		continue
	    }

	    # else...
	    set dn batch.$fi;	# make a temp sub dir
	    cd $tmpd; file mkdir $dn; cd $tmpd/$dn
	    file link link.dcd $fn 
	    file link $pdb ../$pdb; # workaround: driver does not accept paths on cl
	    file link $meta ../$meta

	    set cmd "$driver_path -dcd link.dcd -pdb $pdb -plumed $meta -ncv $cv_n $pbc"
	    puts "Processing $fn"
	    dputs "Spawning $cmd"
	    set io [ open "| $cmd 2>driver.stderr " r ]
	    lappend batch_ongoing $io
	    fconfigure $io -blocking 0
	    fileevent $io readable [ list Plumed::batch_event $io $fi $fn ]
	    incr fi
	} else {
	    dputs "C batch_ongoing $batch_ongoing"
	    update
	    dputs "D batch_ongoing $batch_ongoing"
	}
    }

    while { [llength $batch_ongoing] > 0 } {
	dputs "Still waiting for processes to finish"
	update
	dputs "F batch_ongoing $batch_ongoing"
    }

    catch { cd $owd }

    # cat, prepending file name
    if { $merge_results } {
	.plumedbatch.progress itemconfigure text -text "Merging results..."
	update
	set filelist [ lsort -dictionary [ glob -directory $batch_dir *.dcd ] ]
	set wch [ open "$batch_dir/$batch_extension.all" w ]
	foreach fn $filelist {
	    set fns [ file tail $fn ]
	    set rch [ open "$fn.$batch_extension" r ]
	    while {[gets $rch line] != -1} {
		if { [ regexp {^#} $line ] } { continue }
		puts $wch "$fns $line"
	    }
	    close $rch
	    if { $batch_keep==0 } {
		file delete "$fn.$batch_extension"
	    }
	}
	close $wch
    }

    .plumedbatch.progress itemconfigure text -text "Complete."
}


proc Plumed::batch_abort_do {  } {
    variable batch_ongoing
    set tmpd "[ Plumed::tmpdir ]/vmdplumed.[pid]"
    puts "Cleaning up [llength $batch_ongoing] ongoing processes."
    foreach fh $batch_ongoing {
	catch { close $fh }
    }
    foreach f [ glob $tmpd/batch.* ] {
	file delete -force $f
    }
    set batch_ongoing {}
    return
}

proc Plumed::batch_event { io findex fname } {
    variable batch_abort
    variable batch_nfiles
    variable batch_ok
    variable batch_ongoing
    variable batch_extension
    
    while { [gets $io line] != -1 } { }
    if { [eof $io] } {
	set tmpd "[ Plumed::tmpdir ]/vmdplumed.[pid]"
	set colvarfile "$tmpd/batch.$findex/COLVAR"
	if { ! [file exists $colvarfile] } {
	    puts "$line"
	    puts "-----------"
	    puts "Error processing file $fname, aborting. Check above output."
	    set batch_abort 1
	    batch_abort_do
	} else {
	    puts "File handle $io, index $findex, $fname completed ok"
	    incr batch_ok
	    Plumed::batch_progress $batch_ok $batch_nfiles
	    file rename -force $colvarfile $fname.$batch_extension
	    catch { close $io }
	    file delete -force $tmpd/batch.$findex
	    dputs "A batch_ongoing $batch_ongoing"
	    set batch_ongoing [ lreplace $batch_ongoing [ lsearch $batch_ongoing $io ] [ lsearch $batch_ongoing $io ] ]
	    dputs "B batch_ongoing $batch_ongoing"
	}
    }
}



# ==================================================

proc Plumed::nc_gui { } { 
    variable nc_selA "protein and name CA"
    variable nc_selB ""
    variable nc_cutoff 7
    variable nc_dresid 0
    variable nc_destmol top
    variable nc_groupname nc

    toplevel .plumednc -bd 4
    wm title .plumednc "Native contacts CV"
    pack [ label .plumednc.head1 -text "Insert a CV and group definitions required to define a native contacts CV.
The current frame of the top molecule is taken as the native state." ] -side top -fill x 

    pack [ frame .plumednc.sel1 ] -side top -fill x
    pack [ label .plumednc.sel1.txt -text "Selection 1: " ] -side left -fill x
    pack [ entry .plumednc.sel1.sel -width 50 -textvariable [namespace current]::nc_selA ] -side left -expand 1 -fill x

    pack [ frame .plumednc.sel2 ] -side top -fill x
    pack [ label .plumednc.sel2.txt -text "Selection 2 (optional): " ] -side left -fill x
    pack [ entry .plumednc.sel2.sel -width 40 -textvariable [namespace current]::nc_selB ] -side left -expand 1 -fill x

    pack [ frame .plumednc.cutoff ] -side top -fill x
    pack [ label .plumednc.cutoff.txt -text "Distance cutoff (A): " ] -side left -fill x
    pack [ entry .plumednc.cutoff.sel -width 10 -textvariable [namespace current]::nc_cutoff ] -side left -expand 1 -fill x
    pack [ label .plumednc.cutoff.txt2 -text "      Single selection: |\u0394 resid| \u2265 " ] -side left -fill x
    pack [ entry .plumednc.cutoff.dresid -width 10 -textvariable [namespace current]::nc_dresid ] -side left -expand 1 -fill x
    Plumed::setTooltip .plumednc.cutoff.dresid "Consider contact pairs only if they span at least N monomers in the sequence (by resid attribute). 
   0 - consider all contact pairs;
   1 - ignore contacts within the same residue;
   2 - also ignore contacts between neighboring monomers; and so on."

    pack [ frame .plumednc.destmol ] -side top -fill x
    pack [ label .plumednc.destmol.txt -text "Target molecule ID: " ] -side left -fill x
    pack [ entry .plumednc.destmol.sel -width 10 -textvariable [namespace current]::nc_destmol ] -side left -expand 1 -fill x

    pack [ frame .plumednc.groupname ] -side top -fill x
    pack [ label .plumednc.groupname.txt -text "Prefix for PLUMED groups: " ] -side left -fill x
    pack [ entry .plumednc.groupname.sel -width 20 -textvariable [namespace current]::nc_groupname ] -side left -expand 1 -fill x

    pack [ label .plumednc.preview -text "Click `Count' to compute the number of contacts." ] -side top -fill x 

    pack [ frame .plumednc.act ] -side top -fill x
    pack [ button .plumednc.act.preview -text "Count" -relief raised -command \
	       { Plumed::nc_preview } ] -side left -fill x -expand 1 
    pack [ button .plumednc.act.insert -text "Insert" -relief raised -command \
	       { Plumed::nc_insert } ] -side left -fill x -expand 1 
    pack [ button .plumednc.act.close -text "Close" -relief raised \
	       -command {  destroy .plumednc }   ] -side left -fill x -expand 1
}


proc Plumed::nc_compute { } {
    variable nc_selA
    variable nc_selB
    variable nc_cutoff
    variable nc_destmol
    variable nc_dresid
    
    # See RMSD trajectory tool enhanced with native contacts
    # http://www.multiscalelab.org/utilities/RMSDTTNC
    set sel1 [ atomselect top $nc_selA ] 
    set sel2 [ atomselect $nc_destmol $nc_selA ]
    if { [$sel1 num] != [$sel2 num] } {
	tk_messageBox -title "Error" -parent .plumednc -message "Selection ($nc_selA) has different number of atoms in molecule top ([$sel1 num]) versus $nc_destmol ([$sel2 num])."
	return
    }

    if { $nc_selB != "" } {
	set sel1B [ atomselect top $nc_selB ] 
	set sel2B [ atomselect $nc_destmol $nc_selB ]
	if { [$sel1B num] != [$sel2B num] } {
	    tk_messageBox -title "Error" -parent .plumednc -message "Selection ($nc_selB) has different number of atoms in molecule top ([$sel1B num]) versus $nc_destmol ([$sel2B num])."
	    return
	}
    } else {
	set sel1B 0
	set sel2B 0 
    }

    # mapping index of top -> serials of $nc_destmol
    # sel1 and sel1B are in top
    # sel2 and sel2B are in $nc_destmol
    array set i1_s2 {}
    foreach idx1 [$sel1 get index]   ser2 [$sel2 get serial] {
	set i1_s2($idx1) $ser2
    }

    # Add atoms in selB if intermolecular
    if { $sel1B != 0 } {
	foreach idx1 [$sel1B get index]  ser2 [$sel2B get serial] {
	    set i1_s2($idx1) $ser2
	}
    }

    # Prepare list of resids in reference
    array set i1resid {}
    if { $sel1B == 0 } {
	foreach i [$sel1 get index] resid [$sel1 get resid] {
	    set i1resid($i) $resid
	}
    }

    # Get native contacts (sel1)
    if { $sel1B != 0 } {
	set nclist [transpose [ measure contacts $nc_cutoff $sel1 $sel1B ] ]
	puts "\nDEBUG: reference is [$sel1 text], sel1B is [$sel1B text], nclist has [llength $nclist]"
    } else {
	set nclist [transpose [ measure contacts $nc_cutoff $sel1 ] ]
    }
    set ncref_full [ llength $nclist ]

    # Convert pair list as atomnos, removing close pairs if needed
    set ncl {}
    foreach pairs $nclist {
	set i1 [lindex $pairs 0]
	set i2 [lindex $pairs 1]
	if { $sel1B == 0 } {
	    if { [ expr abs( $i1resid($i1) - $i1resid($i2) ) ] < $nc_dresid } {
		puts "DEBUG: Removing contact pair $i1-$i2 (resid $i1resid($i1) - $i1resid($i2) )"
		continue
	    }
	}
	set p1 $i1_s2($i1)
	set p2 $i1_s2($i2)
	lappend ncl [ list $p1 $p2 ]
    }

    puts "CONTACTS: [ llength $nclist ] in reference, [ llength $ncl ] after removing close resids"
    return $ncl
}

proc Plumed::nc_preview { } {
    .plumednc.preview configure -text "Counting, please wait..."
    update
    set ncl [ Plumed::nc_compute ]
    set ncn [ llength $ncl ]
    puts "NC: $ncl "
    .plumednc.preview configure -text "There are $ncn native contacts."
}


proc Plumed::nc_insert { } {
    variable nc_groupname 
    variable nc_cutoff
    variable w

    set nc [ Plumed::nc_compute ]
    .plumednc.preview configure -text "There are [llength $nc] native contacts."
    if { [llength $nc ] == 0 } {
	tk_messageBox -title "Error" -parent .plumednc -message "There are no contacts in the currently selected frame."
	return
    }
    set ncl [ transpose $nc  ]

    set txt1 "${nc_groupname}_1-> [lindex $ncl 0] ${nc_groupname}_1<-"
    set txt2 "${nc_groupname}_2-> [lindex $ncl 1] ${nc_groupname}_2<-"
    set txt3 "COORD LIST <${nc_groupname}_1> <${nc_groupname}_2> PAIR NN 6 MM 12 D_0 $nc_cutoff R_0 0.5"

    $w.txt.text insert 1.0 "$txt1\n$txt2\n\n"
    $w.txt.text insert insert "$txt3\n"
}


# ==================================================

proc Plumed::do_compute {} {
    variable driver_path
    variable cv_n

    if {[molinfo top]==-1 || [molinfo top get numframes] < 2} {
	tk_messageBox -title "Error" -parent .plumed -message "A top molecule and at least two frames are required to plot."
	return    }

    if {![file executable $driver_path] || ![file isfile $driver_path]} {
	tk_messageBox -title "Error" -parent .plumed -message "The plumed 'driver' executable is required. See manual for installation instructions."
	return     }

    set tmpd "[ Plumed::tmpdir ]/vmdplumed.[pid]"
    file mkdir $tmpd
    set owd [ pwd ]
    cd $tmpd

    set dcd temp.dcd
    animate write dcd $dcd waitfor all

    set pdb temp.pdb
    Plumed::writePlumed [atomselect top all] $pdb

    set meta META_INP
    Plumed::write_meta_inp $meta

    set out COLVAR
    file delete $out
    set pbc [ get_pbc ]

    puts "Executing: $driver_path -dcd $dcd -pdb $pdb -plumed $meta -ncv $cv_n $pbc"
    if { [ catch { eval exec $driver_path -dcd $dcd -pdb $pdb -plumed $meta -ncv $cv_n $pbc } driver_stdout ] } {
	set dontplot 1
    } else {
	set dontplot 0
    }

    puts $driver_stdout
    puts "-----------"
    puts "Temporary files are in directory $tmpd"
    cd $owd

    if { $dontplot } {
	puts "Something went wrong. Check above messages."
    } else {
	Plumed::do_plot "$tmpd/$out" "$driver_stdout"
    }
}


proc Plumed::do_plot { { out COLVAR } { txt ""  } } {
    variable w
    variable cv_n
    variable plot_points
    variable plot_wall

    set fd [open $out r]
    while { [gets $fd header] >= 0 } {
	if { ! [ regexp {^#} $header ] } {
	    # End of header lines
	    set line $header
	    break
	} 
	# puts "header: $header";

	if { [ regexp {FIELDS} $header ] } {
	    # Find column numbers from header
	    set hlist [ split $header ]
	    set tcol [expr  [ lsearch $header time ] - 2 ]
	    set wcol [expr  [ lsearch $header vwall ] - 2 ]
	    for { set i 1 } { $i <= $cv_n } { incr i } {
		set fld   [ lsearch $header "cv$i" ]
		if { $fld != -1 } {
		    set cvcol($i)  [expr $fld - 2 ]
		    if { [ regexp -lineanchor "^$i-(.+?)$" $txt junk cvname ] } {
			set cvlabel($i) "CV $i: $cvname"
		    } else {
			set cvlabel($i) "CV $i"
		    }
		} else {
		    puts "Can't find CV $i in output - not plotting"
		    set cvcol($i) -1
		}
	    }
	}
    }

    if { ![info exists tcol] } {
	puts "Missing FIELDS header. 'driver' executable of PLUMED version > 1.2.0 is required"
	close $fd
	return
    }
	    
    # Read in data in lists
    set ltime {}
    set lwall {}
    while 1 {
	# puts "data: $line"
	set sline [regexp -inline -all -- {\S+} $line];	# split http://wiki.tcl.tk/989
	lappend ltime [lindex $sline $tcol]
	lappend lwall [lindex $sline $wcol]
	for { set i 1 } { $i <= $cv_n } { incr i } {
	    if { $cvcol($i) == -1 } { 	continue }
	    lappend lcv($i) [lindex $sline $cvcol($i)]
	}
	if { [gets $fd line] <0 } break
    }
    close $fd
    
    if { [ llength $ltime ] == 0 } {
	puts "No output in COLVAR. Please check above messages."
	return
    }

    if { [ llength $ltime ] == 1 } {
	puts "Single frame output. Omitting plot"
	return
    }

#    puts "$ltime $lwall"
#    parray lcv
    if { $plot_points } {
	set pt circle
    } else {
	set pt none
    }
    
    set cvplot [multiplot -title "Collective variables" -xlabel "Frame"  -nostats ]

    for { set i 1 } { $i <= $cv_n } {incr i } {
	if { $cvcol($i) == -1 } { continue }
	set coln [ expr $i - 1 ]
	while {$coln > 15} { set coln [expr $coln - 16] }
	set color [index2rgb $coln]
	$cvplot add $ltime $lcv($i) -legend $cvlabel($i) -lines -marker $pt \
	    -radius 2 -fillcolor $color -color $color -nostats
    }

    if { $plot_wall } {
	$cvplot add $ltime $lwall -legend "Wall" -lines -marker $pt \
	    -radius 2 -nostats -dash ","
    }

    $cvplot replot

}




