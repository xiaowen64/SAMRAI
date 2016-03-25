# This class taken from "[incr Tcl/TK] from the Ground Up" 
# by Chad Smith

::itcl::class cvisTrace {
  # Class constructor/destructor
  constructor {} {print EXE}
  destructor {print EXE "" 1}

  # Class Methods
  public method print {type_ {string_ ""} {destructing_ 0}}

  # Static Methods
  public proc ACTIVATE {flags_} {SET_FLAGS $flags_ 1}
  public proc DEACTIVATE {flags_} {SET_FLAGS $flags_ 0}
  private proc SET_FLAGS {flags_ val_}

  # Static Data Members
  private common _ON
    array set _ON "x 0 e 0 w 0 c 0"
  private common _INDENT ""
  private common _PAD "  "
  private common _CHARS
    array set _CHARS "e * w ! c #"
}


#--------------------------------------------------------------------
# PUBLIC PROCEDURE print
#
# Description:
#   Prints each of the four message types: error, warning, comment,
#   and tracing.  The message type is dependant on the type_
#   parameter, which must be one of ERR, WARN, COMM, or EXE.  If it's
#   not EXE, then the string specified by the string_ parameter is
#   printed following the appropriate characters.  Error messages are
#   preceded by * characters, warnings are preceded by ! characters,
#   and comments/informational messages are preceded by # characters.
#   Otherwise, if type_ is EXE, the execution tracing output is
#   printed.  Unlike the Trace class from Chapter 5, the print method
#   is not manually called for tracing output.  Rather, this is
#   handled by Trace2's constructor and destructor, which when invoked
#   call print with the appropriate parameters.  The indentation is
#   adjusted accordingly.
#
# Input Parameters:
#   type_: type of error message to print - ERR, WARN, COMM, or EXE
#   string_: text to print if message type is ERR, WARN, or COMM
#   destructing_: only used by the constructor/destructor to determine
#                 whether to print open/close braces during tracing
#
# Return Values:
#   None
#--------------------------------------------------------------------
::itcl::body cvisTrace::print {type_ {string_ ""} {destructing_ 0}} {
  switch -- $type_ {
    EXE {
      if {$_ON(x)} {
        if {$destructing_} {
          set _INDENT [string range $_INDENT 2 end]
          puts "$_INDENT\}"
        } else {
          # Grab the invoking method name 2 levels down the stack.
          set methodName \
            [lindex [info level [expr [info level] - 2]] 0]
          set fullyQualifiedName \
            [uplevel 2 [list namespace which $methodName]]
          puts "${_INDENT}$fullyQualifiedName"
          puts "$_INDENT\{"
          set _INDENT "$_INDENT$_PAD"
        }
      }
      return
    }
    ERR     {if {$_ON(e)} {set char e} else {return}}
    WARN    {if {$_ON(w)} {set char w} else {return}}
    COMM    {if {$_ON(c)} {set char c} else {return}}
    default {set char c}
  }
  puts "${_INDENT}$_CHARS($char)$_CHARS($char) $type_: $string_"
}


#--------------------------------------------------------------------
# PRIVATE PROCEDURE SET_FLAGS
#
# Description:
#   Invoked by the ACTIVATE and DEACTIVATE procedures to enable and
#   disable different message types.  
#
# Input Parameters:
#   flags_: one or more of x, w, e, or c or can be "all"
#   val_: either 1 or 0
#
# Return Values:
#   None
#--------------------------------------------------------------------
::itcl::body cvisTrace::SET_FLAGS {flags_ val_} {
  if {$flags_ == "all"} {
    foreach type "x e w c" {
      set _ON($type) $val_ 
    }
  } else {
    foreach type "x e w c" {
      if {[string first $type $flags_] != -1} {
        set _ON($type) $val_
      }
    } 
  }
}
