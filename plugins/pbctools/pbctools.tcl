##
## PBC Tools
##
## A plugin for the handling of periodic boundary conditions.
##
## Authors: 
##   Jerome Henin <Jerome.Henin _at_ edam.uhp-nancy.fr>
##   Olaf Lenz <olenz _at_ icp.uni-stuttgart.de>
##   Cameron Mura <cmura _at_ mccammon.ucsd.edu>
##   Jan Saam <saam _at_ charite.de>
##
## The pbcbox procedure copies a lot of the ideas of Axel Kohlmeyer's
## script vmd_draw_unitcell.
##
## $Id: pbctools.tcl,v 1.14 2011/02/20 20:02:18 akohlmey Exp $
##
package provide pbctools 2.6

###################################################
# Main namespace procedures
###################################################
# Main UI
proc pbc { args } {
    proc usage {} {
	vmdcon -info {usage: pbc <command> [args...]

	Setting/getting PBC information:
	  set cell [options...]
	  get [options...]
	  readxst $xstfile [options...]
	
	Drawing a box:
	  box [options...]
	  box_draw [options...]
	
	(Un)Wrapping atoms:
	  wrap [options...]
	  unwrap [options...]
	    }
	return
    }

    if { [llength $args] < 1 } then { usage; return }
    set command [ lindex $args 0 ]
    set args [lrange $args 1 end]
    set fullcommand "::PBCTools::pbc$command"

#     vmdcon -info "command=$command"
#     vmdcon -info "fullcommand=$fullcommand"
#     vmdcon -info "args=$args"

    if { [ string length [namespace which -command $fullcommand]] } then {
	eval "$fullcommand $args"
    } else { usage; return }
}

