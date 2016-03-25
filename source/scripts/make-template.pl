##
## File:        make-template.pl
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2005 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: perl script to create template files for the SAMRAI repository
##

#
# Usage: perl make-template.pl <filename file> <directory> <template-class> <type>
#        perl make-template.pl standard        ./          tbox::Array       double
#        perl make-template.pl char            ./          tbox::Pointer     tbox::Array char*
#
# <filename file> is used to accumulate the names of the files generated for later
# use in constructing the makefile.

# This script should probably be reworked.  It does not handle templates
# on more than one type very well and the NDIM handling is a 
# bit of a kludge.

#
# Read the input directory, calculate the file name, and create the class
#

$OBJFILE = shift(@ARGV);
$DIR   = shift(@ARGV);
$PACKAGE  = shift(@ARGV);
$FILE  = getFileName(@ARGV);
$CLASS = getClassName(@ARGV);
$_     = $CLASS;





# Append the file generated to list of files
open(OUTFILE, ">>$DIR/$OBJFILE") || die "Cannot open filename file $DIR/$OBJFILE";
print OUTFILE "$FILE.C\n";
close(OUTFILE);

#
# List classes require special static data member instantiation
#

if (/^tbox::List</ | /^List</) {
   ($ITER = $CLASS) =~ s/List/ListIterator/;
   ($NODE = $CLASS) =~ s/List/ListNode/;
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   if ($FILE =~ /HDFDatabase|VisItDataWriter/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   print OUTFILE "#ifdef LACKS_STATIC_DATA_INSTANTIATION\n";
   print OUTFILE "$NODE *$NODE\::s_free_list=0;\n";
   print OUTFILE "bool $NODE\::s_registered_callback=false;\n";
   print OUTFILE "#endif\n";
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "template class $ITER;\n";
   print OUTFILE "template class $NODE;\n";
   if ($FILE =~ /HDFDatabase|VisItDataWriter/) {
      print OUTFILE  "#endif\n";
   }
   printFooter();
   close(OUTFILE);

#
# Special template instantiation for PETSc vector classes
#

} elsif (/^(PETScAbstractVectorReal|PETSc_SAMRAIVectorReal|SNES_SAMRAIContext)/) {
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   print OUTFILE "#ifdef HAVE_PETSC\n"; 
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "#endif\n";
   printFooter();
   close(OUTFILE);


#
# Special template instantiation for Hyper classes
#

} elsif (/^(CellPoissonHypreSolver)/) {
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   print OUTFILE "#ifdef HAVE_HYPRE\n"; 
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "#endif\n";
   printFooter();
   close(OUTFILE);


#
# Special template instantiation for KINSOL classes
#

} elsif (/^(KINSOL_SAMRAIContext)/) {
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   print OUTFILE "#ifdef HAVE_KINSOL\n"; 
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "#endif\n";
   printFooter();
   close(OUTFILE);

#
# Special template instantiation for PVODE classes
#

} elsif (/^(PVodeTrio_SAMRAIVector)/) {
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   print OUTFILE "#ifdef HAVE_PVODE\n"; 
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "#endif\n";
   printFooter();
   close(OUTFILE);


} elsif (/^Box\</) {
   ($ITER = $CLASS) =~ s/Box/BoxIterator/;
   ($BASECLASS = $CLASS) =~ s/^(.*)\<(.*)\>/\1/;
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "template class $ITER;\n";
   print OUTFILE "template istream& operator >> (istream& s, $BASECLASS<NDIM>& box);\n";
   print OUTFILE "template ostream& operator << (ostream& s, const $BASECLASS<NDIM>& box);\n";

   if ($FILE =~ /HDF/) {
      print OUTFILE  "#endif\n";
   }
   if ($FILE =~ /appu_VisIt/) {
      print OUTFILE  "#endif\n";
   }
   printFooter();
   close(OUTFILE);
} elsif (/^IntVector\</) {
   ($BASECLASS = $CLASS) =~ s/^(.*)\<(.*)\>/\1/;
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "template istream& operator >> (istream& s, $BASECLASS<NDIM>& box);\n";
   print OUTFILE "template ostream& operator << (ostream& s, const $BASECLASS<NDIM>& box);\n";

   if ($FILE =~ /HDF/) {
      print OUTFILE  "#endif\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#endif\n";
   }
   printFooter();
   close(OUTFILE);
} elsif (/^PatchLevel\</) {
   ($ITER = $CLASS) =~ s/PatchLevel/PatchLevelIterator/;
   ($BASECLASS = $CLASS) =~ s/^(.*)\<(.*)\>/\1/;
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "template class $ITER;\n";
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#endif\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#endif\n";
   }
   printFooter();
   close(OUTFILE);
} elsif (/^LocallyActiveDataPatchLevelManager\</) {
   ($ITER = $CLASS) =~ s/LocallyActiveDataPatchLevelManager/LocallyActiveDataPatchLevelIterator/;
   ($BASECLASS = $CLASS) =~ s/^(.*)\<(.*)\>/\1/;
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   print OUTFILE "template class $CLASS;\n";
   print OUTFILE "template class $ITER;\n";
   printFooter();
   close(OUTFILE);
#
# Other templates do not require special treatment (except for HDF classes and VisIt classes, which use HDF)
#
} else {
   open(OUTFILE, ">$DIR/$FILE.C") || die "Cannot open output file $DIR/$FILE.C";
   printHeader($FILE, @ARGV);
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#ifdef HAVE_HDF5\n";
   }
   print OUTFILE "template class $CLASS;\n";
   if ($FILE =~ /HDF/) {
      print OUTFILE  "#endif\n";
   }
   if ($FILE =~ /VisIt/) {
      print OUTFILE  "#endif\n";
   }
   printFooter();
   close(OUTFILE);
}

#
# Exit the script
#

exit(0);

#
# Construct a filename from the collection of class names
#

sub getFileName {
   my $FILE = shift(@_);

   $FILE =~ s/::/__/g;
   $FILE =~ s/\,/_/g;

   for (@_) {
      (my $TYPE = $_) =~ s/\*/Pointer/g;
      $TYPE =~ s/::/__/g;
      $TYPE =~ s/\<|\>|\,/_/g;

      $FILE .= "-$TYPE";

   }

   # If the file is NDIM based add X to the filename
   # to pickup the correct build rule.
   if ($FILE =~ /NDIM/) {
       $FILE .= X;
   }

   return $FILE;
}

#
# Generate the template class name (with brackets) from the name list
#

sub getClassName {
   my $NAME = shift;
   if ($#_ >= 0) {
      $NAME .= "< ";
      $NAME .= getClassName(@_);
      $NAME .= " >";
   }
   return($NAME);
}

#
# Print header information to OUTFILE (banner, includes, namespace)
#

sub printHeader {
   my $FILE = shift(@_);
   my $CLASS = shift(@_);

   # tbox classes are in a subdirectory in the include dir
   ($CLASSFILE = $CLASS) =~ s/tbox::/tbox\//g;

   # classname comes after the package prefix (if any)
   $CLASSFILE =~ s/(.*)::(.*)/\2/g;

   print OUTFILE "//\n";
   print OUTFILE "// File:\t$FILE.C\n";
   print OUTFILE "// Package:\tSAMRAI templates\n";
   print OUTFILE "// Copyright:\t(c) 1997-2005 The Regents of the University of California\n";
   print OUTFILE "// Revision:\t\$Revision: 1.32 \n";
   print OUTFILE "// Modified:\t\$Date: 2003/01/22 01:12:43 \n";
   print OUTFILE "// Description:\tAutomatically generated template file\n";
   print OUTFILE "//\n";
   print OUTFILE "\n";

   print OUTFILE "#include \"$CLASSFILE.h\"\n";
   print OUTFILE "#include \"$CLASSFILE.C\"\n";

   # There may be comma seperated class lists so split them 
   # apart so we can include each of the headers
   @classes=split(/[\s,]/, join(" ", @_));
   for (@classes) {
       # strip off pointer for classnames
       $_ =~ s/\*//g;
       # tbox classes are in a subdirectory in the include dir
       if (/tbox::/) {
	   $tbox_dir="tbox/";
       } else {
	   $tbox_dir="";
       }
       if (/^NDIM$/) {
	   # NDIM is an int, no header required
       } elsif (/dcomplex|NDIM\,dcomplex/) {
	   print OUTFILE "#include \"tbox/Complex.h\"\n";
       } elsif (/string|NDIM\,string/) {
	   print OUTFILE "#include <string>\nusing namespace std;\n";
       } elsif (/NDIM\,.*/) {
       } elsif (/^(.*)::(.*)\<NDIM\>$/) {
	   print OUTFILE "#include \"$tbox_dir$2.h\"\n";
       } elsif (/^(.*)::(.*)\<NDIM\>::.*/) {
	   print OUTFILE "#include \"$tbox_dir$2.h\"\n";
       } elsif (/^([^:]+)\<NDIM\>::.*/) {
	   print OUTFILE "#include \"$tbox_dir$1.h\"\n";
       } elsif (/^([^:]+)::(.*)::(.*)/) {
	   print OUTFILE "#include \"$tbox_dir$2.h\"\n";
       } elsif (/^([^:]+)::(.*)/) {
	   print OUTFILE "#include \"$tbox_dir$2.h\"\n";
       } elsif (!/^(bool)|(char)|(double)|(float)|(int)$/) {
	   print OUTFILE "#include \"$tbox_dir$_.h\"\n";
       }
   }
   print OUTFILE "\n";

   print OUTFILE "namespace SAMRAI {\n";
   print OUTFILE "   namespace $PACKAGE {\n";
}

#
# Print footer information to OUTFILE (closes the SAMRAI namespace)
#

sub printFooter {
   print OUTFILE "}\n";
   print OUTFILE "}\n";
}
