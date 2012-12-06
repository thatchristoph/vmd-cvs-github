#
# $Id: fftk_BuildPar.tcl,v 1.6 2012/09/19 17:59:14 mayne Exp $
#

#======================================================
namespace eval ::ForceFieldToolKit::BuildPar {

    variable psf
    #variable pdb
    variable RefParList
    variable missingParOutPath
    variable avgParsTemplate
    variable avgParsInput
    variable avgParsOutPath
    variable updateInputParPath
    variable updateLogPath
    variable updateOutParPath

}

#======================================================
proc ::ForceFieldToolKit::BuildPar::init {} {

    # localize variables
    variable psf
    #variable pdb
    variable RefParList
    variable missingParOutPath
    variable avgParsTemplate
    variable avgParsInput
    variable avgParsOutPath
    variable updateInputParPath
    variable updateLogPath
    variable updateOutParPath

    # initialize
    set psf {}
    #set pdb  {}
    set RefParList {}
    set missingParOutPath {}
    set avgParsTemplate {}
    set avgParsInput {}
    set avgParsOutPath {}
    set updateInputParPath {}
    set updateLogPath {}
    set updateOutParPath {}
    
}
#======================================================
proc ::ForceFieldToolKit::BuildPar::sanityCheck { procType } {
    # checks to see that appropriate information is set
    # prior to running any of the BuildPar procedures
    
    # returns 1 if all input is sane
    # returns 0 if there is a problem
    
    # localize buildPar variables
    variable psf
    variable RefParList
    variable missingParOutPath
    variable avgParsTemplate
    variable avgParsInput
    variable avgParsOutPath
    variable updateInputParPath
    variable updateLogPath
    variable updateOutParPath
    
    # local variables
    set errorList {}
    set errorText ""
    
    # build the error list based on what proc is checked
    # and what the existing input is
    switch -exact $procType {
        idMissingPars {
            # make sure that a PSF file was entered and exists
            if { $psf eq "" } { lappend errorList "No PSF file was specified." } \
            else { if { ![file exists $psf] } { lappend errorList "Cannot find PSF file." } }

            # no parameter files is OK

            # make sure that the output path/filename were entered and user has write permissions to the directory
            if { $missingParOutPath eq "" } { lappend errorList "No output path was specificed." }
            if { $missingParOutPath ne "" && ![file writable [file dirname $missingParOutPath]] } { lappend errorList "Cannot write to the specified directory" }
        }

        avgBAPars {
            # make sure that a template parameter file was entered and exists
            if { $avgParsTemplate eq "" } { lappend errorList "No template parameter file was specified." } \
            else { if { ![file exists $avgParsTemplate] } { lappend errorList "Cannot find template parameter file." } }

            # make sure that a paratool parameter file was entered and exists
            if { $avgParsInput eq "" } { lappend errorList "No paratool parameter file was specified." } \
            else { if { ![file exists $avgParsInput] } { lappend errorList "Cannot find paratool parameter file." } }

            # make sure that a savename was entered and that the user has write permissions
            if { $avgParsOutPath eq "" } { lappend errorList "No output path was specified." }
            if { $avgParsOutPath ne "" && ![file writable [file dirname $avgParsOutPath]] } { lappend errorList "Cannot write to output path." }
        }

        updateOptPars {
            # make sure that an input parameter file was entered and exists
            if { $updateInputParPath eq "" } { lappend errorList "No input parameter file was specified." } \
            else { if { ![file exists $updateInputParPath] } { lappend errorList "Cannot find input parameter file." } }

            # make sure that an optimization log file was entered and exists
            if { $updateLogPath eq "" } { lappend errorList "No optimization log file was specified." } \
            else { if { ![file exists $updateLogPath] } { lappend errorList "Cannot find optimization log file." } }

            # make sure a savename was entered and that the user can write to that directory
            if { $updateOutParPath eq "" } { lappend errorList "No output path was specified." }
            if { $updateOutParPath ne "" && ![file writable [file dirname $updateOutParPath]] } { lappend errorList "Cannot write to output path." }

        }
    }

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
proc ::ForceFieldToolKit::BuildPar::getRefPars { prmlist } {
    # reads in a list of parameter files and returns
    # lists of type definitions present in input prm files

    # initialize some lists
    set bonds {}
    set angles {}
    set dihedrals {}
    set impropers {}
    set vdws {}
    
    set partypes {bonds angles dihedrals impropers vdws}
    
    foreach prmfile $prmlist {
        # read in parameters from file
        set tempParList [::ForceFieldToolKit::SharedFcns::readParFile $prmfile]
        # cycle through parameter sections
        for {set i 0} {$i <= [llength $partypes]} {incr i} {
            # cycle through each parameter definition/entry
            foreach parDef [lindex $tempParList $i] {
                # append the typeDef to the appropriate list
                lappend [lindex $partypes $i] [lindex $parDef 0]
            }
        }
    }
    
    # remove duplicates (mostly to reduce list size that is stored/passed)
    set bonds [lsort -unique $bonds]
    set angles [lsort -unique $angles]
    set dihedrals [lsort -unique $dihedrals]
    set impropers [lsort -unique $impropers]
    set vdws [lsort -unique $vdws]

    # return the relevent lists for building a supplemental
    # zero-ed out prm file
    return [list $bonds $angles $dihedrals $vdws]
}
#======================================================
proc ::ForceFieldToolKit::BuildPar::getMolecPars { psfFile } {
    # finds bond, angle, dihedral, vdw parameters that are
    # required for a molecule
    
    # load the psf file
    mol new $psfFile
    
    # retype
    ::ForceFieldToolKit::SharedFcns::reTypeFromPSF $psfFile "top"
    
    # initilize some variables
    set uniq_bond_list {}
    set uniq_ang_list {}
    set uniq_dih_list {}
    set uniq_vdw_list {}

    # bonds
    set bond_list [topo getbondlist]
    foreach bond_entry $bond_list {
        set at1 [[atomselect top "index [lindex $bond_entry 0]"] get type]
        set at2 [[atomselect top "index [lindex $bond_entry 1]"] get type]
        set bond_type [list $at1 $at2]
        # test forward and reverse patterns for duplicate bond types
        if {[lsearch -exact $uniq_bond_list $bond_type] != -1} {
            #puts "bond: $bond_type is a forward duplicate"
        } elseif {[lsearch -exact $uniq_bond_list [lreverse $bond_type]] != -1} {
            #puts "bond: $bond_type is a reverse duplicate"
        } else {
            #puts "bond: $bond_type is unique"
            lappend uniq_bond_list $bond_type
        }
    }

    # angles
    set ang_list [topo getanglelist]
    foreach ang_entry $ang_list {
        set at1 [[atomselect top "index [lindex $ang_entry 1]"] get type]
        set at2 [[atomselect top "index [lindex $ang_entry 2]"] get type]
        set at3 [[atomselect top "index [lindex $ang_entry 3]"] get type]
        set ang_type [list $at1 $at2 $at3]
        # test forward and reverse patterns for duplicate angle types
        if {[lsearch -exact $uniq_ang_list $ang_type] != -1} {
            #puts "angle: $ang_type is a forward duplicate"
        } elseif {[lsearch -exact $uniq_ang_list [lreverse $ang_type]] != -1} {
            #puts "angle: $ang_type is a reverse duplicate"
        } else {
            #puts "angle: $ang_type is unique"
            lappend uniq_ang_list $ang_type
        }
    }

    # dihedrals
    set dih_list [topo getdihedrallist]
    foreach dih_entry $dih_list {
        set at1 [[atomselect top "index [lindex $dih_entry 1]"] get type]
        set at2 [[atomselect top "index [lindex $dih_entry 2]"] get type]
        set at3 [[atomselect top "index [lindex $dih_entry 3]"] get type]
        set at4 [[atomselect top "index [lindex $dih_entry 4]"] get type]
        set dih_type [list $at1 $at2 $at3 $at4]
        # test forward and reverse patterns for duplicate dihedral types
        if {[lsearch -exact $uniq_dih_list $dih_type] != -1} {
            #puts "dihedral: $dih_type is a forward duplicate"
        } elseif {[lsearch -exact $uniq_dih_list [lreverse $dih_type]] != -1} {
            #puts "dihdedral: $dih_type is a reverse duplicate"
        } else {
            #puts "dihedral: $dih_type is unique"
            lappend uniq_dih_list $dih_type
        }
    }

    # vdws
    set uniq_vdw_list [topo atomtypenames]
    
    # clean up
    mol delete top
    
    # return uniq parameters present in molecule
    return [list $uniq_bond_list $uniq_ang_list $uniq_dih_list $uniq_vdw_list]
}
#======================================================
proc ::ForceFieldToolKit::BuildPar::crossCheckPars { molecPars refPars } {
    # cross checks molecule pars against a reference set of pars
    # returns molec pars not found in reference par set
    
    # initialize some lists
    set unpar_bonds {}
    set unpar_angles {}
    set unpar_diheds {}
    set unpar_vdws {}

    for {set i 0} {$i<4} {incr i} {
        foreach prmtype [lindex $molecPars $i] {
            if {[lsearch -exact [lindex $refPars $i] $prmtype] != -1} {
                #puts "\tparameter: $prmtype has forward parameter
                continue
            } elseif {[lsearch -exact [lindex $refPars $i] [lreverse $prmtype]] != -1} {
                #puts "\tparameter: $prmtype has reverse parameter
                continue
            } else {
                switch -exact $i {
                    0 {lappend unpar_bonds $prmtype}
                    1 {lappend unpar_angles $prmtype}
                    2 {lappend unpar_diheds $prmtype}
                    3 {lappend unpar_vdws $prmtype}
                }; # end switch
            }; # end prm search test (else)
        }; # end cycling through actual prm lists (foreach loop)
    }; # end cycling through each prm type (for loop)

    # return the missing parameters
    return [list $unpar_bonds $unpar_angles $unpar_diheds $unpar_vdws]
}
#======================================================
proc ::ForceFieldToolKit::BuildPar::buildInitParFile {} {
    # builds a zeroed-out parameter list for a molecule
    # containing only unparameterized parameters

    # localize some variables
    variable psf
    variable RefParList
    variable missingParOutPath
    
    # run a sanity check
    if { ![::ForceFieldToolKit::BuildPar::sanityCheck idMissingPars ] } { return }

    # read in the type definitions for reference parameter set    
    set refPars [::ForceFieldToolKit::BuildPar::getRefPars $RefParList]

    # read in the type definitions for the molecule parameter set
    set molecPars [::ForceFieldToolKit::BuildPar::getMolecPars $psf]

    # search the reference parameters for molecule parametes
    set missingPars [::ForceFieldToolKit::BuildPar::crossCheckPars $molecPars $refPars]
    # note: returns { bonds angles dihedrals vdws }
    
    # format the missing parameters so that we can pass to a generic par writer proc
    set Bonds {}
    foreach bondDef [lindex $missingPars 0] {
        lappend Bonds [list $bondDef {0.0 0.0} {}]
    }
    set Angles {}
    foreach angleDef [lindex $missingPars 1] {
        lappend Angles [list $angleDef {0.0 0.0} {} {}]
    }
    set Dihedrals {}
    foreach dihDef [lindex $missingPars 2] {
        lappend Dihedrals [list $dihDef {0.0 1 0.0} {}]
    }
    # impropers (not set here)
    set vdws {}
    foreach vdwDef [lindex $missingPars 3] {
        lappend vdws [list $vdwDef {0.0 0.0} {} {! SET BY ANALOGY!!!}]
    }
    
    # build the parlist
    set parList [list $Bonds $Angles $Dihedrals {} $vdws]
    
    # write the zeroed out par file
    ::ForceFieldToolKit::SharedFcns::writeParFile $parList $missingParOutPath
}
#======================================================
# DEPRECIATED WITH ADDITION OF fftk_GenBonded.tcl
#
#proc ::ForceFieldToolKit::BuildPar::buildAvgParFile {} {
#    # reads in a template parameter file
#    # reads in the par file generated by paratool's analysis of Hessian
#    # averages any duplicate parameters (bond, angle, dih) from the paratool prm file
#    # writes out only averaged parameters that are present in template prm file
#    
#    # localize variables
#    variable avgParsTemplate
#    variable avgParsInput
#    variable avgParsOutPath
#    
#    # run a sanity check
#    if { ![::ForceFieldToolKit::BuildPar::sanityCheck avgBAPars ] } { return }
#    
#    
#    # setup the template parameters based on init file
#    # read in the template parameter file
#    set templatePars [::ForceFieldToolKit::SharedFcns::readParFile $avgParsTemplate]
#    
#    # build an array for each parameter type
#    array set templateBonds {}
#    foreach bond [lindex $templatePars 0] {
#        #                    {type def}                {k  b0}          {comment}
#        set templateBonds([lindex $bond 0]) [list [lindex $bond 1] [lindex $bond 2]]
#    }
#    array set templateAngles {}
#    foreach angle [lindex $templatePars 1] {
#        #                       {type def}              {k  theta}       {ksub   s}         {comment}
#        set templateAngles([lindex $angle 0]) [list [lindex $angle 1] [lindex $angle 2] [lindex $angle 3]]
#    }
#    array set templateDihs {}
#    foreach dih [lindex $templatePars 2] {
#        #                  {type def}            {k mult delta}     {comment}
#        set templateDihs([lindex $dih 0]) [list [lindex $dih 1] [lindex $dih 2]]
#    }
#    # skip impropers
#    # skip vdw
#    
#   
#    # read in the paratool-generated parameter file
#    set inputPars [::ForceFieldToolKit::SharedFcns::readParFile $avgParsInput]
#        
#    # BONDS
#    # initialize some variables
#    set b_list {}
#    set b_rts {}
#
#    # parse out bonds section from full parameter set
#    set bonds [lindex $inputPars 0]
#    # cycle through each bond definition
#    foreach bondEntry $bonds {
#        # parse out parameter data
#        # { {bond type def} {k b} {comment} }
#        set typeDef [lindex $bondEntry 0]
#        set k [lindex $bondEntry 1 0]
#        set b [lindex $bondEntry 1 1]
#        
#        # test (forward and reverse)
#        set testfwd [lsearch -exact $b_list $typeDef]
#        set testrev [lsearch -exact $b_list [lreverse $typeDef]]
#        
#        if { $testfwd == -1 && $testrev == -1 } {
#            # new bond type definition, append all values
#            lappend b_list $typeDef
#            lappend b_rts [list $k $b 1]
#        } else {
#            if { $testfwd > -1 } {
#                set ind $testfwd
#            } else {
#                set ind $testrev
#            }
#            # repeat type definition found, add to running totals
#            lset b_rts $ind 0 [expr {[lindex $b_rts $ind 0] + $k}]
#            lset b_rts $ind 1 [expr {[lindex $b_rts $ind 1] + $b}]
#            lset b_rts $ind 2 [expr {[lindex $b_rts $ind 2] + 1}]
#        }
#    }
#
#    # update the bonds array
#    for {set i 0} {$i < [llength $b_list]} {incr i} {
#        # if the paratool bond is present in the template, update the k and b0 values
#        if { [info exists templateBonds([lindex $b_list $i])] } {
#            # calc the avg values from running totals data (b_rts)
#            set avgK [expr {[lindex $b_rts $i 0]/[lindex $b_rts $i 2]}]
#            set avgB0 [expr {[lindex $b_rts $i 1]/[lindex $b_rts $i 2]}]
#            # update the value in the templateBonds array
#            lset templateBonds([lindex $b_list $i]) 0 [list $avgK $avgB0]     
#        }
#        # else, we don't really care
#    }
#    
#    # rebuild the bonds section
#    set newbonds {}
#    foreach key [array names templateBonds] {
#        lappend newbonds [list $key [lindex $templateBonds($key) 0] [lindex $templateBonds($key) 1]]
#    }
#    
#    # replace the bonds section of the template
#    lset templatePars 0 $newbonds
#
#    # DONE with BONDS
#
#
#   
#    # ANGLES
#    # initialize some variables
#    set a_list {}
#    set a_rts {}
#    set a_rtsub {}
#    
#    # parse out angles section from paratool parameter set
#    set angles [lindex $inputPars 1]
#    # cycle through each angle definition
#    foreach angleEntry $angles {
#        # parse out parameter data
#        set typeDef [lindex $angleEntry 0]
#        set k [lindex $angleEntry 1 0]
#        set theta [lindex $angleEntry 1 1]
#        set kub [lindex $angleEntry 2 0]
#        set s [lindex $angleEntry 2 1]
#        
#        # test (forward and reverse)
#        set testfwd [lsearch -exact $a_list $typeDef]
#        set testrev [lsearch -exact $a_list [lreverse $typeDef]]
#        if { $testfwd == -1 && $testrev == -1 } {
#            # new angle definition, append all data
#            lappend a_list $typeDef
#            lappend a_rts [list $k $theta 1]
#            # handle angle UB term, empty or not
#            if { $kub ne "" } {
#                lappend a_rtsub [list $kub $s 1]
#            } else {
#                lappend a_rtsub [list {} {} 0]
#            }
#        } else {
#            # duplicate definition, update running totals and count
#            if { $testfwd > -1 } {
#                set ind $testfwd
#            } else {
#                set ind $testrev
#            }
#            # update angle totals and count
#            lset a_rts $ind 0 [expr {[lindex $a_rts $ind 0] + $k}]
#            lset a_rts $ind 1 [expr {[lindex $a_rts $ind 1] + $theta}]
#            lset a_rts $ind 2 [expr {[lindex $a_rts $ind 2] + 1}]
#            # for UB term, update if not empty string (just ignore empty strings)
#            if { $kub ne "" } {
#                # how the term is updated depends on whether there are any UB terms stored already
#                # i.e., if the count is above 0 we need to update, otherwise, just replace
#                if { [lindex $a_rtsub $ind 2] > 0 } {
#                    lset a_rtsub $ind 0 [expr {[lindex $a_rtsub $ind 0] + $kub}]
#                    lset a_rtsub $ind 1 [expr {[lindex $a_rtsub $ind 1] + $s}]
#                    lset a_rtsub $ind 2 [expr {[lindex $a_rtsub $ind 2] + 1}]
#                } else {
#                    lset a_rtsub $ind 0 $kub
#                    lset a_rtsub $ind 1 $s
#                    lset a_rtsub $ind 2 1
#                }
#            }; # end of UB term if
#        }; # end of angles test
#    }; # end of angles loop
#    
#    # update the angles array
#    for {set i 0} {$i < [llength $a_list]} {incr i} {
#        # if the paratool angle is present in the template, update the values
#        if { [info exists templateAngles([lindex $a_list $i])] } {
#            # calc the avg values from angle running totals
#            set avgK [expr {[lindex $a_rts $i 0]/[lindex $a_rts $i 2]}]
#            set avgTheta [expr {[lindex $a_rts $i 1]/[lindex $a_rts $i 2]}]
#            # update the angles data in the angles array
#            lset templateAngles([lindex $a_list $i]) 0 [list $avgK $avgTheta]
#            
#            # if kub and s are defined (count greater than zero), average them
#            # otherwise set as undefined
#            if { [lindex $a_rtsub $i 2] > 0 } {
#                set avgKub [expr {[lindex $a_rtsub $i 0]/[lindex $a_rtsub $i 2]}]
#                set avgS [expr {[lindex $a_rtsub $i 1]/[lindex $a_rtsub $i 2]}]
#                # update with value
#                lset templateAngles([lindex $a_list $i]) 1 [list $avgKub $avgS]
#            } else {
#                # update as undefined
#                lset templateAngles([lindex $a_list $i]) 1 [list {} {}]
#            }
#        }
#        # else, we don't care about the values
#    }
#    
#    # rebuild the angles section
#    set newangles {}
#    foreach key [array names templateAngles] {
#        lappend newangles [list $key [lindex $templateAngles($key) 0] [lindex $templateAngles($key) 1] [lindex $templateAngles($key) 2]]
#    }
#    
#    # replace the angles section of the template
#    lset templatePars 1 $newangles
#    
#    # DONE with ANGLES
#    
#    
#    # DIHEDRALS
#    # there are several ways to consider how to handle
#    # here we will average all dihedral k and delta values, take only last mult value
#    
#    set d_list {}
#    set d_rts {}
#    
#    # parse out the dihedrals section from full parameter set
#    set dihedrals [lindex $inputPars 2]
#    # cycle through each dihedral definition
#    foreach dihEntry $dihedrals {
#        # parse out parameter data
#        set typeDef [lindex $dihEntry 0]
#        set k [lindex $dihEntry 1 0]
#        set n [lindex $dihEntry 1 1]
#        set delta [lindex $dihEntry 1 2]
#        
#        # test (forward and reverse)
#        set testfwd [lsearch -exact $d_list $typeDef]
#        set testrev [lsearch -exact $d_list [lreverse $typeDef]]
#        if { $testfwd == -1 && $testrev == -1 } {
#            # new dihedral definition, append all data
#            lappend d_list $typeDef
#            lappend d_rts [list $k $n $delta 1]
#        } else {
#            # repeat dih definition, add to running totals
#            if { $testfwd > -1 } {
#                set ind $testfwd
#            } else {
#                set ind $testrev
#            }
#            lset d_rts $ind 0 [expr {[lindex $d_rts $ind 0] + $k}]
#            lset d_rts $ind 1 $n; # note that n is updated, but NOT averaged
#            lset d_rts $ind 2 [expr {[lindex $d_rts $ind 2] + $delta}]
#            lset d_rts $ind 3 [expr {[lindex $d_rts $ind 3] + 1}]
#        }
#    }
#    
#    # update the dihs array
#    for {set i 0} {$i < [llength $d_list]} {incr i} {
#        # if the paratool dihedral is present in the template, update the values
#        if { [info exists templateDihs([lindex $d_list $i])] } {
#            # calc the averages   
#            set avgK [expr {[lindex $d_rts $i 0]/[lindex $d_rts $i 3]}]
#            set avgN [lindex $d_rts $i 1]
#            set avgDelta [expr {[lindex $d_rts $i 2]/[lindex $d_rts $i 3]}]
#            # update the values in the templateDihs array
#            lset templateDihs([lindex $d_list $i]) 0 [list $avgK $avgN $avgDelta]
#        }
#        # else, we don't really care about the values
#    }
#    
#    # rebuild the dih section
#    set newdihs {}
#    foreach key [array names templateDihs] {
#        lappend newdihs [list $key [lindex $templateDihs($key) 0] [lindex $templateDihs($key) 1]]
#    }
#    
#    # replace the dih section of the template
#    lset templatePars 2 $newdihs
#    
#    # DONE with DIHEDRALS
#    
#    # keep all other template data
#    
#    # write the new parameter file
#    ::ForceFieldToolKit::SharedFcns::writeParFile $templatePars $avgParsOutPath
#        
#}
#======================================================
proc ::ForceFieldToolKit::BuildPar::buildUpdatedParFile {} {
    # reads in a parameter file and optimization log (bonds/angles, dihedral)
    # writes a parameter file with updated parameters

    # localize variables
    variable updateInputParPath
    variable updateLogPath
    variable updateOutParPath
    
    # initialize some variables
    set optBonds {}
    set optAngles {}
    set optDihedrals {}

    # run a sanity check
    if { ![::ForceFieldToolKit::BuildPar::sanityCheck updateOptPars ] } { return }

    # read in the template par file
    set parData [::ForceFieldToolKit::SharedFcns::readParFile $updateInputParPath]
    
    # parse the optimization log file
    set inFile [open $updateLogPath r]
    set readState 0
    
    while { ![eof $inFile] } {
        # read a line at a time
        set inLine [gets $inFile]
        
        # determine if we've reached the final parameter data
        switch -exact $inLine {
            "FINAL PARAMETERS" { set readState 1 }
            "END" { set readState 0 }
            default {
                if { $readState } {
                    # parse and append parameters to appropriate list
                    switch -exact [lindex $inLine 0] {
                        "bond" { lappend optBonds [list [lindex $inLine 1] [list [lindex $inLine 2] [lindex $inLine 3]] {}] }
                        "angle" { lappend optAngles [list [lindex $inLine 1] [list [lindex $inLine 2] [lindex $inLine 3]] {} {}] }
                        "dihedral" { lappend optDihedrals [list [lindex $inLine 1] [list [lindex $inLine 2] [lindex $inLine 3] [lindex $inLine 4]] {}] }
                        default { continue }
                    }
                } else {
                    continue
                }
            }
        }; # end outer switch
    }; # end of while loop (parse log file)
    
    # update bonds
    # build a list of bond definitions from the input par data
    set oldBondDefs {}
    foreach bondEntry [lindex $parData 0] {
        lappend oldBondDefs [lindex $bondEntry 0]
    }
    # cycle through each bond with parameters to update
    foreach bond2update $optBonds {
        # search against the input parameter data
        set testfwd [lsearch $oldBondDefs [lindex $bond2update 0]]
        set testrev [lsearch $oldBondDefs [lreverse [lindex $bond2update 0]]]
        if { $testfwd == -1 && $testrev ==-1 } {
            puts "ERROR: Bond definition to update: [lindex $bond2update 0] was not found in input parameter set"
        } elseif { $testfwd > -1 } {
            lset parData 0 $testfwd 1 [lindex $bond2update 1]
        } elseif { $testrev > -1 } {
            lset parData 0 $testrev 1 [lindex $bond2update 1]
        }
    }
    
    # update angles
    # build a list of angle definitions from the input par data
    set oldAngleDefs {}
    foreach angleEntry [lindex $parData 1] {
        lappend oldAngleDefs [lindex $angleEntry 0]
    }
    # cycle through each angle with parameters to update
    foreach angle2update $optAngles {
        # search against the input parameter data
        set testfwd [lsearch $oldAngleDefs [lindex $angle2update 0]]
        set testrev [lsearch $oldAngleDefs [lreverse [lindex $angle2update 0]]]
        if { $testfwd == -1 && $testrev == -1 } {
            puts "ERROR: Angle definition to update: [lindex $angle2update 0] was not found in input parameter set"
        } elseif { $testfwd > -1 } {
            lset parData 1 $testfwd 1 [lindex $angle2update 1]
        } elseif { $testrev > -1 } {
            lset parData 1 $testrev 1 [lindex $angle2update 1]
        }
    }
    
    # update dihedrals
    # due to multiplicities, we have to handle this one differently
    
    # initialize and populate array for old dihedral information
    # oldDihPars(type def) = {  {k1 mult1 delta1} ...{kN multN deltaN}  }
    array set oldDihPars {}
    #array set oldDihCom {}
    foreach dihEntry [lindex $parData 2] {
        lappend oldDihPars([lindex $dihEntry 0]) [lindex $dihEntry 1]
        #lappend oldDihCom([lindex $dihEntry 0]) [lindex $dihEntry 2]
    }
    
    # cycle through each dihedral with parameters to update
    array set newDihPars {}
    foreach dih2update $optDihedrals {
        # search against the input parameter data and determine the type def order (fwd or rev)
        # if found in the old parameter set, blow up the old definition
        # if not found in old paremeter set, check the new one
        # if not found in either, then it's an error of some sort
        if { [info exists oldDihPars([lindex $dih2update 0])] } {
            set currTypeDef [lindex $dih2update 0]
            array unset oldDihPars [lindex $dih2update 0]
        } elseif { [info exists dihPars([lreverse [lindex $dih2update 0]])] } {
            set currTypeDef [lreverse [lindex $dih2update 0]]
            array unset oldDihPars [lreverse [lindex $dih2update 0]]
        } elseif { [info exists newDihPars([lindex $dih2update 0])] } {
            set currTypeDef [lindex $dih2update 0]
        } elseif { [info exists newDihPars([lreverse [lindex $dih2update 0]])] } {
            set currTypeDef [lreverse [lindex $dih2update 0]]
        } else {
            puts "ERROR: Dihedral definition to update: [lindex $dih2update 0] was not found in input parameter set"
            continue
        }
        
        # add to the new array
        lappend newDihPars($currTypeDef) [lindex $dih2update 1]
        #lappend newDihCom($currTypeDef) [lindex $dih2update 2]
    }
    
    # convert dihedral data array back into parData format
    set dihParUpdate {}
    # cycle through any remaining old parameters
    foreach key [array names oldDihPars] {
        # cycle through each mult/comment for a given typeDef
        for {set i 0} {$i < [llength $oldDihPars($key)]} {incr i} {
            #                          {def}      {k mult delta}             {comment}
            #lappend dihParUpdate [list $key [lindex $oldDihPars($key) $i] [lindex $oldDihCom($key) $i]]
            lappend dihParUpdate [list $key [lindex $oldDihPars($key) $i] {}]
        }
    }
    # cycle through new parameters
    foreach key [array names newDihPars] {
        # cycle through each mult/comment for a given typeDef
        for {set i 0} {$i < [llength $newDihPars($key)]} {incr i} {
            #                          {def}      {k mult delta}             {comment}
            #lappend dihParUpdate [list $key [lindex $newDihPars($key) $i] [lindex $newDihCom($key) $i]]
            lappend dihParUpdate [list $key [lindex $newDihPars($key) $i] {}]
        }
    }
    
    # replace the input dihedral parameters with the updated parameters
    lset parData 2 $dihParUpdate    


    
    # write the updated parameter file
    ::ForceFieldToolKit::SharedFcns::writeParFile $parData $updateOutParPath
}
#======================================================
