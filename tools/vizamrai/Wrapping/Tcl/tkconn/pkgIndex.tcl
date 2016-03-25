# tkcon package is from Jeffrey Hobbs, jeff.hobbs@acm.org
# http://www.purl.org/net/hobbs/tcl/script/tkcon/
# http://www.hobbs.wservice.com/tcl/script/tkcon/
# See ../tkcon directory for information on license
package ifneeded tkcon 1.5 [list tclPkgSetup $dir tkcon 1.5 { \
        {tkcon.tcl source {tkcon}} \
    }]

