############################################################################
#cr
#cr            (C) Copyright 1995-2009 The Board of Trustees of the
#cr                        University of Illinois
#cr                         All Rights Reserved
#cr
############################################################################

############################################################################
# RCS INFORMATION:
#
#       $RCSfile: mdff_setup.tcl,v $
#       $Author: ltrabuco $        $Locker:  $             $State: Exp $
#       $Revision: 1.5 $       $Date: 2010/02/07 05:01:27 $
#
############################################################################


# MDFF package
# Authors: Leonardo Trabuco <ltrabuco@ks.uiuc.edu>
#          Elizabeth Villa <villa@ks.uiuc.edu>
#
# mdff gridpdb -- creates a pdb file docking
# mdff setup   -- writes a NAMD config file for docking
#

package require readcharmmpar
package require exectool
package provide mdff_setup 0.2

namespace eval ::MDFF::Setup:: {

  variable defaultDiel 80
  variable defaultScaling1_4 1.0
  variable defaultGScale 0.3
  variable defaultTemp 300
  variable defaultNumSteps 500000
  variable defaultMinimize 200
  variable defaultConsCol B
  variable defaultFixCol O
  variable defaultMargin 0

  variable defaultK 10.0
  variable defaultGridpdbSel {noh and (protein or nucleic)}

}

proc ::MDFF::Setup::mdff_setup_usage { } {
    
  variable defaultDiel
  variable defaultScaling1_4
  variable defaultGScale
  variable defaultTemp
  variable defaultParFile
  variable defaultNumSteps
  variable defaultMinimize
  variable defaultConsCol 
  variable defaultFixCol
  variable defaultMargin 

  ::MDFF::Setup::init_files

  puts "Usage: mdff setup -o <output prefix> -psf <psf file> -pdb <pdb file> -griddx <griddx file> ?options?"
  puts "Options:" 
  puts "  -gridpdb  -- pdb file for docking (default: -pdb)"
  puts "  -diel     -- dielectric constant (default: $defaultDiel; 1 with -pbc)" 
  puts "  -temp     -- temperature in Kelvin (default: $defaultTemp)" 
  puts "  -ftemp    -- final temperature (default: $defaultTemp)" 
  puts "  -gscale   -- scaling factor for the grid (default: $defaultGScale)" 
  puts "  -extrab   -- extrabonds file (default: none)" 
  puts "  -conspdb  -- pdb file with constrained atoms (default: none)"
  puts "  -conscol  -- force constant column in conspdb (default: beta)"
  puts "  -fixpdb   -- pdb file with fixed atoms (default: none)"
  puts "  -fixcol   -- column in fixpdb (default: occupancy)"
  puts "  -scal14   -- 1-4 scaling (default: $defaultScaling1_4)"
  puts "  -step     -- docking protocol step (default: 1)" 
  #puts "  -parfiles -- parameter file list (default $defaultParFile)"
  puts "  -parfiles -- parameter file list"
  puts "  -minsteps -- number of minimization steps (default $defaultMinimize)"
  puts "  -numsteps -- number of time steps to run (default: $defaultNumSteps)" 
  puts "  -margin   -- extra length in patch dimension during simulation (default: $defaultMargin)"
  puts "  -pbc      -- use periodic boundary conditions (for explicit solvent)"
  
}

proc ::MDFF::Setup::mdff_setup { args } {

  variable defaultDiel 
  variable defaultScaling1_4
  variable defaultGScale 
  variable defaultTemp 
  variable defaultParFile
  variable defaultNumSteps 
  variable defaultMinimize
  variable defaultConsCol 
  variable defaultFixCol
  variable namdTemplateFile
  variable defaultMargin 

  set nargs [llength $args]
  if {$nargs == 0} {
    mdff_setup_usage
    error ""
  }

  # Get NAMD template and parameter files
  ::MDFF::Setup::init_files

  # Periodic simulation?
  set pos [lsearch -exact $args {-pbc}]
  if { $pos != -1 } {
    set pbc 1
    set args [lreplace $args $pos $pos]
  } else {
    set pbc 0
  }

  # parse switches
  foreach {name val} $args {
    switch -- $name {
      -o          { set arg(o)        $val }
      -psf        { set arg(psf)      $val } 
      -pdb        { set arg(pdb)      $val }
      -gridpdb    { set arg(gridpdb)  $val }
      -diel       { set arg(diel)     $val }
      -scal14     { set arg(scal14)   $val }
      -temp       { set arg(temp)     $val }
      -ftemp      { set arg(ftemp)    $val }
      -griddx     { set arg(griddx)   $val }
      -gscale     { set arg(gscale)   $val }
      -extrab     { set arg(extrab)   $val }
      -conspdb    { set arg(conspdb)  $val }
      -conscol    { set arg(conscol)  $val }
      -fixpdb     { set arg(fixpdb)   $val }
      -fixcol     { set arg(fixcol)   $val }
      -step       { set arg(step)     $val }
      -parfiles   { set arg(parfiles) $val }
      -numsteps   { set arg(numsteps) $val }
      -minsteps   { set arg(minsteps) $val }
      -margin     { set arg(margin)   $val }
    }
  }
    
  if { [info exists arg(o)] } {
    set outprefix $arg(o)
  } else {
    mdff_setup_usage
    error "Missing output files prefix."
  }
  
  if { [info exists arg(psf)] } {
    set psf $arg(psf)
  } else {
    mdff_setup_usage
    error "Missing psf file."
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    mdff_setup_usage
    error "Missing pdb file."
  }

  if { [info exists arg(diel)] } {
    set diel $arg(diel)
  } elseif $pbc {
    set diel 1
  } else {
    set diel $defaultDiel
  }

  if { [info exists arg(scal14)] } {
    set scal14 $arg(scal14)
  } else {
    set scal14 $defaultScaling1_4
  }

  if { [info exists arg(temp)] } {
    set itemp $arg(temp)
  } else {
    set itemp $defaultTemp
  }

  if { [info exists arg(ftemp)] } {
    set ftemp $arg(ftemp)
  } else {
    set ftemp $itemp
  }

  if { [info exists arg(parfiles)] } {
    set parfiles $arg(parfiles)
  } else {
    file copy -force $defaultParFile .
    set parfiles [file tail $defaultParFile]
  }

  if { [info exists arg(numsteps)] } {
    set numsteps $arg(numsteps)
  } else {
    set numsteps $defaultNumSteps
  }


  if { [info exists arg(minsteps)] } {
    set minsteps $arg(minsteps)
  } else {
    set minsteps $defaultMinimize
  }

  if { [info exists arg(griddx)] } {
    set grid $arg(griddx) 
  } else {
    mdff_setup_usage 
    error "Missing grid dx file name."
  }
  
  if { [info exists arg(gscale)] } {
    set gscale $arg(gscale)
  } else {
    set gscale $defaultGScale
  }
  
  if { [info exists arg(extrab)] } {
    set extrab $arg(extrab)
  } else {
    set extrab 0
  }
  
  if { [info exists arg(gridpdb)] } {
    set gridpdb $arg(gridpdb)
  } else {
    set gridpdb $pdb
  }

  if { [info exists arg(conspdb)] } {
    set conspdb $arg(conspdb)
  } else {
    set conspdb 0
  }

  if { [info exists arg(conscol)] } {
    set conscol $arg(conscol)
  } else {
    set conscol $defaultConsCol
  }

  if { [info exists arg(fixpdb)] } {
    set fixpdb $arg(fixpdb)
  } else {
    set fixpdb 0
  }

  if { [info exists arg(fixcol)] } {
    set fixcol $arg(fixcol)
  } else {
    set fixcol $defaultFixCol
  }

  if { [info exists arg(minsteps)] } {
    set minsteps $arg(minsteps)
  } else {
    set minsteps $defaultMinimize
  }

  if { [info exists arg(margin)] } {
    set margin $arg(margin)
  } else {
    set margin $defaultMargin
  }

  if { [info exists arg(step)] } {
    set step $arg(step)
  } else {
    # puts "No step number was specified. Assuming step 1.."
    set step 1
  }

  # Copy NAMD template file to current directory
  file copy -force $namdTemplateFile .

  set outname ${outprefix}-step${step}
  puts "mdff) Writing NAMD configuration file ${outname}.namd ..."
  
  set out [open ${outname}.namd w]     

  puts $out "###  Docking -- Step $step" 
  puts $out " "   
  puts $out "set PSFFILE $psf"        
  puts $out "set PDBFILE $pdb"
  puts $out "set GRIDPDB $gridpdb"
  puts $out "set DIEL $diel"        
  puts $out "set SCALING_1_4 $scal14"
  puts $out "set ITEMP $itemp"   
  puts $out "set FTEMP $ftemp"   
  puts $out "set GRIDFILE $grid"   
  puts $out "set GSCALE $gscale"   
  puts $out "set EXTRAB [list $extrab]"   
  puts $out "set CONSPDB $conspdb"
  if {$conspdb != "0" } {
    puts $out "set CONSCOL $conscol"
  }
  puts $out "set FIXPDB  $fixpdb"
  if {$fixpdb != "0" } {
    puts $out "set FIXCOL $fixcol"
  }

  puts $out " " 
  
  if {$step >  1 } {
    set prevstep [expr $step - 1]
    set inputname "${outprefix}-step${prevstep}"
    set prevnamd "${inputname}.namd"
    if { ![file exists $prevnamd] } {
      puts "Warning: Previous NAMD configuration file $prevnamd not found." 
      puts "You may need to manually edit the variable INPUTNAME in the file ${outname}.namd."
    }
    puts $out "set INPUTNAME $inputname"  
  }

  puts $out "set OUTPUTNAME $outname"
  puts $out " "
  puts $out "set TS $numsteps"
  puts $out "set MS $minsteps"
  puts $out " "
  puts $out "set MARGIN $margin"
  puts $out " "
  puts $out "####################################"
  puts $out " "
  puts $out "structure \$PSFFILE"
  puts $out "coordinates \$PDBFILE"
  puts $out " "
  puts $out "paraTypeCharmm on"
  foreach par $parfiles {
    puts $out "parameters $par"
  }
  if $pbc {
    puts $out ""
    puts $out "if {\[info exists INPUTNAME\]} {"
    puts $out "  BinVelocities \$INPUTNAME.restart.vel"
    puts $out "  BinCoordinates \$INPUTNAME.restart.coor"
    puts $out "  ExtendedSystem \$INPUTNAME.restart.xsc"
    puts $out "} else {"
    puts $out "  temperature \$ITEMP"
    ::MDFF::Setup::get_cell $psf $pdb $out
    puts $out "}"

    puts $out "PME yes"
    puts $out "PMEGridSpacing 1.0"
    puts $out "PMEPencils 1"

    puts $out "wrapAll on"

  } else {
    puts $out ""
    puts $out "if {\[info exists INPUTNAME\]} {"
    puts $out "  BinVelocities \$INPUTNAME.restart.vel"
    puts $out "  BinCoordinates \$INPUTNAME.restart.coor"
    puts $out "} else {"
    puts $out "  temperature \$ITEMP"
    puts $out "}"

  }
  puts $out " "
  puts $out "source [file tail $namdTemplateFile]"

  close $out

}

proc ::MDFF::Setup::get_cell {psf pdb out} {
  set molid [mol new $psf type psf waitfor all]
  mol addfile $pdb type pdb waitfor all

  set sel [atomselect $molid {noh water}]

  if { [$sel num] == 0 } {
    $sel delete
    mol delete $molid
    error "Could not determine the periodic cell information. No water molecules were found in the input structure."
  }
  set minmax [measure minmax $sel]
  set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
  puts $out "  cellBasisVector1 [lindex $vec 0] 0 0"
  puts $out "  cellBasisVector2 0 [lindex $vec 1] 0"
  puts $out "  cellBasisVector3 0 0 [lindex $vec 2]"
  set center [measure center $sel]
  puts $out "  cellOrigin $center"
  $sel delete
  
  mol delete $molid

}


proc ::MDFF::Setup::mdff_gridpdb_usage { } {
 
  variable defaultGridpdbSel

  puts "Usage: mdff gridpdb -psf <input psf> -pdb <input pdb> -o <output pdb> ?options?"
  puts "Options:" 
  puts "  -seltext   -- atom selection text  (default: $defaultGridpdbSel)"
}

proc ::MDFF::Setup::mdff_gridpdb { args } {

  variable defaultGridpdbSel

  set nargs [llength $args]
  if {$nargs == 0} {
    mdff_gridpdb_usage
    error ""
  }

  # parse switches
  foreach {name val} $args {
    switch -- $name {
      -psf        { set arg(psf)      $val }
      -pdb        { set arg(pdb)      $val }
      -o          { set arg(o)        $val }
      -seltext    { set arg(seltext)  $val }
    }
  }
    
  if { [info exists arg(o)] } {
    set gridpdb $arg(o)
  } else {
    mdff_gridpdb_usage
    error "Missing output gridpdb file name."
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    mdff_gridpdb_usage
    error "Missing pdb file."
  }

  if { [info exists arg(psf)] } {
    set psf $arg(psf)
  } else {
    mdff_gridpdb_usage
    error "Missing psf file."
  }

  if { [info exists arg(seltext)]} {
    set seltext $arg(seltext)
  } else {
    set seltext $defaultGridpdbSel
  }
  
  set molid [mol new $psf type psf waitfor all]
  mol addfile $pdb type pdb waitfor all
  set all [atomselect $molid all]
  $all set occupancy 0

  if { $seltext == "all" } {
    $all set beta [$all get mass]
    $all set occupancy 1
  } else {
    $all set beta 0
    set sel [atomselect $molid $seltext]
    if {[$sel num] == 0} {
      error "empty atomselection"
    } else {
      $sel set occupancy 1
      $sel set beta [$sel get mass]
    }  
    $sel delete
  }

  $all writepdb $gridpdb
  $all delete
  
  return 

}

proc ::MDFF::Setup::init_files {} {
  global env
  variable defaultParFile
  variable namdTemplateFile

  set defaultParFile [file join $env(CHARMMPARDIR) par_all27_prot_lipid_na.inp]
  set namdTemplateFile [file join $env(MDFFDIR) mdff_template.namd]

}

proc ::MDFF::Setup::mdff_constrain_usage { } {

  variable defaultK
  variable defaultConsCol
 
  puts "Usage: mdff constrain <atomselection> -o <pdb file> ?options?"
  puts "Options:"
  puts "  -col <column> (default: $defaultConsCol)"
  puts "  -k <force constant in kcal/mol/A^2> (default: $defaultK)"
  
}

proc ::MDFF::Setup::mdff_fix_usage { } {

  variable defaultFixCol
 
  puts "Usage: mdff fix <atomselection> -o <pdb file> ?options?"
  puts "Options:"
  puts "  -col <column> (default: $defaultFixCol)"
  
}

proc ::MDFF::Setup::mdff_constrain { args } {

  variable defaultK
  variable defaultConsCol

  set nargs [llength $args]
  if {$nargs == 0} {
    mdff_constrain_usage
    error ""
  }

  set sel [lindex $args 0]
  if { [$sel num] == 0 } {
    error "empty atomselection"
  }

  foreach {name val} [lrange $args 1 end] {
    switch -- $name {
      -o     { set arg(o)     $val }
      -col   { set arg(col)   $val }
      -k     { set arg(k)     $val }
      default { error "unkown argument: $name $val" }
    }
  }

  if { [info exists arg(o)] } {
    set outputFile $arg(o)
  } else {
    mdff_constrain_usage
    error "Missing output pdb file."
  }

  if { [info exists arg(col)] } {
    set col $arg(col)
  } else {
    set col $defaultConsCol
  }

  if { $col == "beta" || $col == "B" } {
    set col "beta"
  } elseif { $col == "occupancy" || $col == "O" } {
    set col "occupancy"
  } elseif { $col == "x" || $col == "X" } {
    set col "x"
  } elseif { $col == "y" || $col == "Y" } {
    set col "y"
  } elseif { $col == "z" || $col == "Z" } {
    set col "z"
  } else {
    error "Unrecognized column."
  }

  if { [info exists arg(k)] } {
    set k $arg(k)
  } else {
    set k $defaultK
  }

  set molid [$sel molid]
  set all [atomselect $molid all]
  set bakCol [$all get $col]
  $all set $col 0
  $sel set $col $k
  $all writepdb $outputFile
  $all set $col $bakCol
  $all delete

  return

}

proc ::MDFF::Setup::mdff_fix { args } {

  variable defaultFixCol

  set nargs [llength $args]
  if {$nargs == 0} {
    mdff_fix_usage
    error ""
  }

  set sel [lindex [lindex $args 0] 0]
  if { [$sel num] == 0 } {
    error "mdff_constrain: empty atomselection."
  }

  foreach {name val} [lrange $args 1 end] {
    switch -- $name {
      -o     { set arg(o)     $val }
      -col   { set arg(col)   $val }
      default { error "unkown argument: $name $val" }
    }
  }

  if { [info exists arg(o)] } {
    set outputFile $arg(o)
  } else {
    mdff_fix_usage
    error "Missing output pdb file."
  }

  if { [info exists arg(col)] } {
    set col $arg(col)
  } else {
    set col $defaultFixCol
  }

  if { $col == "beta" || $col == "B" } {
    set col "beta"
  } elseif { $col == "occupancy" || $col == "O" } {
    set col "occupancy"
  } elseif { $col == "x" || $col == "X" } {
    set col "x"
  } elseif { $col == "y" || $col == "Y" } {
    set col "y"
  } elseif { $col == "z" || $col == "Z" } {
    set col "z"
  } else {
    error "Unrecognized column."
  }

  set molid [$sel molid]
  set all [atomselect $molid all]
  set bakCol [$all get $col]
  $all set $col 0
  $sel set $col 1
  $all writepdb $outputFile
  $all set $col $bakCol
  $all delete

  return

}
