#
# $Id: fftk_GenDihScan.tcl,v 1.8 2012/09/21 19:39:42 mayne Exp $
#

#======================================================
namespace eval ::ForceFieldToolKit::GenDihScan {

    variable psf
    variable pdb
    variable outPath
    variable basename
    variable qmProc
    variable qmCharge
    variable qmMem
    variable qmMult
    variable qmRoute
    
    variable dihData
    
}
#======================================================
proc ::ForceFieldToolKit::GenDihScan::init {} {
    
    # localize variables
    variable psf
    variable pdb
    variable outPath
    variable basename
    variable qmProc
    variable qmCharge
    variable qmMem
    variable qmMult
    variable qmRoute
    
    variable dihData

    # set variables
    set psf {}
    set pdb {}
    set outPath {}
    set basename {}
    ::ForceFieldToolKit::GenDihScan::resetGaussianDefaults
    #set qmProc 1
    #set qmCharge 0
    #set qmMem 1
    #set qmMult 1
    #set qmRoute "# opt=modredundant MP2/6-31g(d)"
    
    set dihData {}
    
}
#======================================================
proc ::ForceFieldToolKit::GenDihScan::sanityCheck {} {
    # checks to see that the appropriate information is set prior to running
    
    # returns 1 if all input is sane
    # returns 0 if there is a problem
    
    # localize relevant variables
    variable psf
    variable pdb
    variable outPath
    variable basename
    variable qmProc
    variable qmCharge
    variable qmMem
    variable qmMult
    variable qmRoute
    variable dihData
    
    # local variables
    set errorList {}
    set errorText ""
    
    # checks
    # psf
    if { $psf eq "" } {
        lappend errorList "No PSF file was specified."
    } else {
        if { ![file exists $psf] } { lappend errorList "Cannot find PSF file." }
    }
    
    # pdb
    if { $pdb eq "" } {
        lappend errorList "No PDB file was specified."
    } else {
        if { ![file exists $pdb] } { lappend errorList "Cannot find PDB file." }
    }

    # make sure that output folder is specified and writable
    if { $outPath eq "" } {
        lappend errorList "No output path was specified."
    } else {
        if { ![file writable $outPath] } { lappend errorList "Cannot write to output path." }
    }
    
    # make sure that basename is not empty
    if { $basename eq "" } { lappend errorList "No basename was specified." }

    # validate dihData
    if { [llength $dihData] == 0 } {
        lappend errorList "No dihedrals were entered for scanning."
    } else {
        foreach ele $dihData {
            # dih indices
            if { [llength [lindex $ele 0]] != 4 } { lappend errorList "Found inappropriate dihedral definition." }
            foreach ind [lindex $ele 0] {
                if { $ind < 0 || ![string is integer $ind] } { lappend errorList "Found inappropriate dihedral index." }
            }
            # plus/minus
            if { ![string is double [lindex $ele 1]] } { lappend errorList "Found inappropriate dihedral +/- value." }
            # step size
            if { [lindex $ele 2] <= 0 || ![string is double [lindex $ele 2]] } { lappend errorList "Found inappropriate dihedral step size." }
        }
    }
    
    # validate gaussian settings (not particularly vigorous validation)
    # qmProc (processors)
    if { $qmProc eq "" } { lappend errorList "No processors were specified." }
    if { $qmProc <= 0 || $qmProc != [expr int($qmProc)] } { lappend errorList "Number of processors must be a positive integer." }
    # qmMem (memory)
    if { $qmMem eq "" } { lappend errorList "No memory was specified." }
    if { $qmMem <= 0 || $qmMem != [expr int($qmMem)]} { lappend errorList "Memory must be a postive integer." }
    # qmCharge (charge)
    if { $qmCharge eq "" } { lappend errorList "No charge was specified." }
    if { $qmCharge != [expr int($qmCharge)] } { lappend errorList "Charge must be an integer." }
    # qmMult (multiplicity)
    if { $qmMult eq "" } { lappend errorList "No multiplicity was specified." }
    if { $qmMult < 0 || $qmMult != [expr int($qmMult)] } { lappend errorList "Multiplicity must be zero or a positive integer." }
    # qmRoute (route card for gaussian; just make sure it isn't empty)
    if { $qmRoute eq "" } { lappend errorList "Route card is empty." }


    # if there is an error, tell the user about it
    # return -1 to tell the calling proc that there is a problem
    if { [llength $errorList] > 0 } {
        foreach ele $errorList {
            set errorText [concat $errorText\n$ele]
        }
        tk_messageBox \
            -type ok \
            -icon warning \
            -message "Application halting due to the following errors:" \
            -detail $errorText
        
        # there are errors, return the error response
        return 0
    }

    # if you've made it this far, there are no errors
    return 1

}
#======================================================
proc ::ForceFieldToolKit::GenDihScan::buildGaussianFiles {} {
    # builds gaussian input files for scanning dihedral angles
    
    # localize variables
    variable psf
    variable pdb
    variable outPath
    variable basename
    variable qmProc
    variable qmCharge
    variable qmMem
    variable qmMult
    variable qmRoute
    
    variable dihData
    
    # run an input sanity check
    if { ![::ForceFieldToolKit::GenDihScan::sanityCheck] } { return }

    # assign Gaussian atom names and gather x,y,z for output com file
    mol new $psf; mol addfile $pdb
    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }
    
    # cycle through each dihedral to scan
    set scanCount 1
    foreach dih $dihData {
        # change 0-based indices to 1-based
        set zeroInds [lindex $dih 0]
        set oneInds {}
        foreach ind $zeroInds {
            lappend oneInds [expr {$ind + 1}]
        }
        
        # negative scan
        # open the output file
        set outfile [open ${outPath}/${basename}.scan${scanCount}.neg.gau w]
        
        # write the header
        puts $outfile "%chk=${basename}.scan${scanCount}.neg.chk"
        puts $outfile "%nproc=$qmProc"
        puts $outfile "%mem=${qmMem}GB"
        puts $outfile "$qmRoute"
        puts $outfile ""
        puts $outfile "$basename Dihedral Scan at MP2/6-31G*"
        puts $outfile ""
        puts $outfile "$qmCharge $qmMult"
        # write coords
       foreach atom_entry $atom_info {
           puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
       }
       # write scan
       puts $outfile ""
       puts $outfile "D $oneInds S [expr int([expr [lindex $dih 1]/[lindex $dih 2]])] [format "%.6f" [expr {-1*[lindex $dih 2]}]]"
       
       close $outfile
       
       # positive scan
        # open the output file
        set outfile [open ${outPath}/${basename}.scan${scanCount}.pos.gau w]
        
        # write the header
        puts $outfile "%chk=${basename}.scan${scanCount}.pos.chk"
        puts $outfile "%nproc=$qmProc"
        puts $outfile "%mem=${qmMem}GB"
        puts $outfile "$qmRoute"
        puts $outfile ""
        puts $outfile "$basename Dihedral Scan at MP2/6-31G*"
        puts $outfile ""
        puts $outfile "$qmCharge $qmMult"
        # write coords
       foreach atom_entry $atom_info {
           puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
       }
       # write scan
       puts $outfile ""
       puts $outfile "D $oneInds S [expr int([expr [lindex $dih 1]/[lindex $dih 2]])] [format "%.6f" [lindex $dih 2]]"
       
       close $outfile    
       
       incr scanCount
        
    }
    
    # clean up
    mol delete top

}
#======================================================
#======================================================
proc ::ForceFieldToolKit::GenDihScan::resetGaussianDefaults {} {
    # resets gaussian settings to the default values

    # localize variables
    variable qmProc
    variable qmCharge
    variable qmMem
    variable qmMult
    variable qmRoute

    # set variables
    set qmProc 1
    set qmCharge 0
    set qmMem 1
    set qmMult 1
    set qmRoute "# opt=modredundant MP2/6-31g(d) Geom=PrintInputOrient"
}
#======================================================
