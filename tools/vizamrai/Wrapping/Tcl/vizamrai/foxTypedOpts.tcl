#!/bin/false -- do not execute!
#
# typedopts -- parse command line options in TCL
#
# USAGE:
#   typedopts <arg-list> <opt-list> <opt-ret> <arg-ret> <rest>
#     [ <options> ]
#
# OPTIONS:
#   -noinit
#     Don't initialize <opt-ret> and <arg-ret>
#
# typedopts is a command line option parser.  It reads <arg-list> and
# finds all the options that are described in <opt-list>.  If no errors
# are found, then typedopts will <rest> to the command line arguments
# that weren't parsed, set <opt-ret> and <arg-ret> as described below,
# and return 1.  <opt-ret> and <arg-ret> will be set to arrays, indexed
# by the names of the options in <opt-list>.  <opt-ret> will contain
# boolean values indicating which options were found in <arg-list>, and
# <arg-ret> will contain the arguments to the options found.
#
# <opt-list> is a TCL list that describes all of the options that the
# program accepts.  It consists of type-name pairs: { <type> <name>
# <type> <name> ... }.  <type> describes what type of argument the
# option may take; the types are described below.  Some of the types
# require additional parameters; for these types, <type> must be a TCL
# list with the type name followed by the parameters.
#
# <name> is the name of the option; if the program accepts more than
# one option of type <type>, then <name> may be a TCL list of all the
# options of that type.  The option may be abbreviated on the command
# line, as long as the abbreviation uniquely identifies one option.
#
# The types for <type> are listed here.  The type name in <type> may be
# a unique abbreviation of the type name listed.
#
#     boolean
#       A simple flag; no argument is accepted for the option.
#
#     string
#       A string -- no type checking is done.
#
#     float
#       A floating point number.  All of the TCL floating point formats
#       should be accepted for this.
#
#     integer
#       An integer.
#
#     number
#       A floating point number or an integer.
#
#     one-of <value>...
#       One of a specific set of values.  The values may be abbreviated
#       when used on the command line.
#
#     list-of <type>
#       A list of values of type <type>.  The list ends when one of the
#       following is found:  A valid command line option, the string
#       "--", or the end of <arg-list>.
#
#     multiple <type>
#       This option takes one argument of type <type>, but may appear
#       more than once on the command line.  If an option is not
#       specified as being of type multiple... then it may appear only
#       once.
#
# If an option is of type list-of... or multiple..., then the value
# found for that option in <arg-ret> will be a TCL list of all the
# values found for the option, in the order that they appeared on the
# command line.
#
# If typedopts finds an option that is not described in <opt-list>, or
# if it finds an option argument of the wrong type, it will set
# <arg-ret>(_ERROR_) to an error message, set <rest> to the rest of the
# options and arguments, and return 0.
#
# If the -noinit option is given to typedopts, then the <opt-ret> and
# <arg-ret> will _not_ be initialized.  This allows the program to call
# typedopts several times with different <arg-list>s without losing the
# information from previous calls.
#
# if typedopts can't parse its options for any reason, it will print an
# error message to stderr and return a -1 without modifying any other
# variables.
#
# EXAMPLE:
#
# The command line parser for a simple text formatter is given below.
# The formatter accepts the options -top, -bottom, -left, -right, and
# -paragraph to set the margins, -header to set the header string,
# -pagenum to set the page number, and -position to position the page
# number in the footer (the position can be left, right, or center).
# It first parses arguments from the environment variable TFORMAT, then
# from the command line.  The command line can have options and
# filenames intermixed; the options affect all files found after the
# option.
#
# proc parseOpts { } {
#
# # global variables needed:  env = environ variables array,
# #                           argv == command line args
#   global env argv
#
# # The options list: (they have to be declared multiple, because we
# # aren't re-initializing the option arrays each time.
#   set optlist {
#     { multiple integer } { left right top bottom paragraph pagenum }
#     { multiple string } header
#     { multiple one-of left right center } position
#   }
#
# # check if we have a $TFORMAT environment variable to parse
#   if { [ lsearch -exact [ array names $env ] "TFORMAT" ] > -1 } then {
#
# # initialize the options arrays found() and values() with the values
# # from TFORMAT
#     set list $env(TFORMAT)
#     while { ! [ typedopts $list $opts found values list ] } {
#
# # error returned from typedopts:  print the error message and
# # continue parsing
#       puts stderr "Parsing \$TFORMAT: $values(_ERROR_)"
#     }
#
# # check if there are any arguments left; if so, give a warning.
#   if { [ llength $list ] } then {
#     puts stderr "Warning:  \$TFORMAT has non-option arguments!"
#   }
#
#   } else {
#
# # initialize options arrays found() and values() from an empty list
#     typedopts { } $opts found values
#   }
#
# # start parsing the command line.  As long as its not empty, we first
# # pass pass it through the option parser, then call the formatter on
# # the files.
#   while { [ llength $argv ] } {
#     while { ! [ typedopts $argv $opts found values argv -noinit ] } {
#       puts stderr "$values(_ERROR_)"
#     }
#     format [ lindex $argv 0 ]
#     set argv [ lrange $argv 1 end ]
#   }
# }
#
# REVISION HISTORY
#
# $LastChangedRevision: 1973 $
#   $LastChangedBy: smithsg $
#     $LastChangedDate: 2008-02-11 16:39:15 -0800 (Mon, 11 Feb 2008) $
#
# Revision 1.1  2001/12/06 23:02:59  smithsg
# Move files for VTK4.0
#
# Revision 1.1  1999/10/01 23:43:03  ssmith
# Added command line parsing
# Load default camera position
# Add low memory option and set polydatamappers to save memory if requested
#
#   Revision 1.0  1994/02/19  22:04:23  darkfox
#   Initial revision
#

proc typedopts { args } {

  proc abbr { s1 s2 } {
    if { [ set len [ string length $s1 ]] } then {
      if { ! [ string compare $s1 [ string range $s2 0 [ expr $len - 1 ]]] } then {
        return 1
      }
    }
    return 0
  }

  proc findabbr { list val } {
    set list [ lsort $list ]
    if { [ set pos [ lsearch -exact $list $val ]] > -1 } then {
      return [ lindex $list $pos ]
    }
    if { [ set pos [ lsearch -glob $list "$val*" ]] > -1 } then {
      if { [ abbr $val [ set realval [ lindex $list $pos ]]] } then {
        if { ! [ abbr $val [ lindex $list [ incr pos ]]] } then {
          return $realval
        }
      }
    }
    return ""
  }

  proc shift { listname } {
    upvar $listname list
    set ret [ lindex $list 0 ]
    set list [ lrange $list 1 end ]
    return $ret
  }

  proc extract { list args } {
    foreach arg $args {
      upvar $arg var
      set var [ shift list ]
    }
    return $list
  }

  proc parseFormats { fmts var } {
    foreach fmt $fmts {
      if { [ regexp $fmt $var ] } then {
        return 1
      }
    }
    return 0
  }

  proc parseOptionType { type listname retname } {
    upvar $listname args
    upvar $retname var

    set ifmt {
      "^\[+-\]?0x\[0-9a-fA-F\]+\$"
      "^\[+-\]?0\[0-7\]+\$"
      "^\[+-\]?\[0-9\]+\$"
    }

    set ffmt {
      "^\[+-\]?\.\[0-9\]+(\[Ee\]\[+-\]?\[0-9\]*)?\$"
      "^\[+-\]?\[0-9\]+\.\[0-9\]*(\[Ee\]\[+-\]?\[0-9\]*)?\$"
      "^\[+-\]?\[0-9\]+\[Ee\]\[+-\]?\[0-9\]*\$"
    }

    set nfmt [ concat $ifmt $ffmt ]

    set otype $type
    switch -exact [ shift type ] {
      b {
        set var ""
        return 1
      }
      i {
        if { [ llength $args ] } then {
          set val [ shift args ]
          if { [ parseFormats $ifmt $val ] } then {
            set var $val
            return 1
          }
        }
        set var "requires an integer argument."
        return 0
      }
      f {
        if { [ llength $args ] } then {
          set val [ shift args ]
          if { [ parseFormats $ffmt $val ] } then {
            set var $val
            return 1
          }
        }
        set var "requires a floating-point argument."
        return 0
      }
      n {
        if { [ llength $args ] } then {
          set val [ shift args ]
          if { [ parseFormats $nfmt $val ] } then {
            set var $val
            return 1
          }
        }
        set var "requires a numeric argument."
        return 0
      }
      s {
        if { [ llength $args ] } then {
          set var [ shift args ]
          return 1
        }
        set var "requires a string argument."
        return 0
      }
      o {
        if { [ llength $args ] } then {
          if { [ string length [ set val [ findabbr $type [ shift args ]]]] } then {
            set var $val
            return 1
          }
        }
        set var "requires a string argument."
        return 0
      }
      m {
        return [ parseOptionType $type args var ]
      }
      l {
        set val ""
        while { [ llength $args ] && ! [ string match "-*" $args ] } {
          if { ! [ parseOptionType $type args ret ] } then {
            set var $ret
            return 0
          }
          lappend val $ret
        }
        set var $val
        return 1
      }
      default {
        puts stderr "Eek!  Option type <$otype> not supported yet!"
        set var "isn't a supported type."
        return 0
      }
    }
  }

  proc parseOption { optlist } {
    set type [ shift optlist ]

    switch -exact [ findabbr { "booleans" "integers" "numbers" "floats" "strings" "one-of" "list-of" "multiple" } $type ] {
      "booleans" -
      "integers" -
      "numbers" -
      "floats" -
      "strings" {
        if { [ llength $optlist ] } then {
          puts stderr "typedopts:  Type $type doesn't take arguments"
          return ""
        }
        return [ string index $type 0 ]
      }
      "one-of" {
        if { ! [ llength $optlist ] } then {
          puts stderr "typedopts:  No arguments given to type $type"
          return ""
        }
        return [ concat [ string index $type 0 ] $optlist ]
      }
      "list-of" -
      "multiple" {
        if { ! [ llength $optlist ] } then {
          puts stderr "typedopts:  No arguments given to type $type"
          return ""
        }
        if { ! [ string length [ set subtype [ parseOption $optlist ]]] } then {
          return ""
        }
        return [ concat [ string index $type 0 ] $subtype ]
      }
      default {
        puts stderr "typedopts:  Unknown option type $type"
        return ""
      }
    }
  }

  set doinit 1

  if { [ llength $args ] < 5 } then {
    puts stderr "typedopts: bad number of arguments."
    return -1
  }

  set args [ extract $args arglist optlist optret argret restret ]

  while { [ llength $args ] } {
    set opt [ shift args ]
    switch -exact [ findabbr { -noinitialize } $opt ] {
      -noinitialize {
        set doinit 0
      }
      default {
        puts stderr "typedopts: bad option \"$opt\": should be -noinitialize or --"
        return -1
      }
    }
  }

  upvar $optret _opts
  upvar $argret _args
  upvar $restret _rest

  set allopts ""

  set type ""

  foreach word $optlist {
    set word [ string trim $word ]
    if { [ string length $type ] } then {
      foreach arg $word {
        if { [ lsearch -exact $arg $allopts ] > -1 } then {
          puts stderr "typedopts: option -$arg multiply declared."
          return -1
        }
        lappend allopts $arg
        set opttype($arg) $type
      }
      set type ""
    } else {
      if { ! [ string length [ set type [ parseOption $word ]]] } then {
        return -1
      }
    }
  }

  if { $doinit } then {
    foreach opt $allopts {
      set _opts($opt) 0
      set _args($opt) ""
    }
  }

set _args(_ERROR_) ""

  set retval 1

  while { [ llength $arglist ] } {
    switch -glob -- $arglist {
      -- {
        shift arglist
        break
      }
      -* {
      }
      * {
        break
      }
    }
    set opt [ string range [ shift arglist ] 1 end ]
    if { [ string length [ set fnd [ findabbr $allopts $opt ]]] } then {
      set type $opttype($fnd)
      if { [ parseOptionType $opttype($fnd) arglist arg ] } then {
        if { $_opts($fnd) && ! [ string match "m*" $type ] } then {
          set _args(_ERROR_) "Found multiple occurrences of option -$fnd"
          set retval 0
          break
        }
        set _opts($fnd) 1
        set _args($fnd) $arg
      } else {
        set _args(_ERROR_) "Option -$fnd $arg"
        set retval 0
        break
      }
    } else {
      set _args(_ERROR_) "Unknown option: -$opt"
      set retval 0
      break
    }
  }

  set _rest $arglist

  return $retval
}
