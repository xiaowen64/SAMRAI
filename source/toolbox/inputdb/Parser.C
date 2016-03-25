//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/Parser.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Parser that reads the input database grammar
//

#include "tbox/Parser.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/Parser.I"
#endif

#ifndef NULL
#define NULL (0)
#endif

extern int yyparse();
extern void yyrestart(FILE*);

extern void parser_static_table_initialize();


namespace SAMRAI {
   namespace tbox {

Parser *Parser::s_default_parser = NULL;
bool Parser::s_static_tables_initialized = 0;

/*
*************************************************************************
*									*
* The constructor creates an unitialized parser object.  All of the	*
* interesting work is done by member function parse().			*
*									*
*************************************************************************
*/

Parser::Parser()
{
   if (!s_static_tables_initialized) {
      parser_static_table_initialize();
      s_static_tables_initialized = 1;
   }
}

/*
*************************************************************************
*									*
* The destructor automatically deallocates the parser object data.	*
*									*
*************************************************************************
*/

Parser::~Parser()
{
}

/*
*************************************************************************
*									*
* Begin parsing the input database file.  Return the number of errors	*
* encountered in the parse.						*
*									*
*************************************************************************
*/

int Parser::parse(
   const std::string& filename,
   FILE* fstream,
   Pointer<Database> database)
{
   d_errors     = 0;
   d_warnings   = 0;

   // Find the path in the filename, if one exists
   std::string::size_type slash_pos = filename.find_last_of( '/' );
   if(slash_pos == std::string::npos) {
      d_pathname   = "";
   } else {
      d_pathname   = filename.substr(0, slash_pos+1);      
   }
   
   ParseData pd;
   pd.d_filename   = filename;
   pd.d_fstream    = fstream;
   pd.d_linenumber = 1;
   pd.d_cursor     = 1;
   pd.d_nextcursor = 1;
   d_parse_stack.clearItems();
   d_parse_stack.addItem(pd);

   d_scope_stack.clearItems();
   d_scope_stack.addItem(database);

   s_default_parser = this;
   yyrestart(NULL);
   if (yyparse() && (d_errors == 0)) {
      error("Unexpected parse error");
   }
   s_default_parser = NULL;

   d_parse_stack.clearItems();
   d_scope_stack.clearItems();

   return(d_errors);
}

/*
*************************************************************************
*									*
* Advance the cursor to the next line in the current input file.	*
*									*
*************************************************************************
*/

void Parser::advanceLine(const int nline)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();
   pd.d_linenumber += nline;
   pd.d_cursor      = 1;
   pd.d_nextcursor  = 1;
}

/*
*************************************************************************
*									*
* Advance the cursor position by the token in the specified string.	*
* Tabs are expanded assuming tab stops at eight character markers.	*
*									*
*************************************************************************
*/

void Parser::advanceCursor(const std::string& token)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();
   pd.d_cursor = pd.d_nextcursor;
   for (std::string::const_iterator i = token.begin(); i != token.end(); i++) {
      if (*i == '\t') {
         pd.d_nextcursor = ((pd.d_nextcursor + 7) & (~7)) + 1;
      } else {
         pd.d_nextcursor++;
      }
   }
}

/*
*************************************************************************
*									*
* Print out errors to pout and track the number of errors.		*
*									*
*************************************************************************
*/

void Parser::error(const std::string& message)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();

   pout << "Error in " << pd.d_filename << " at line " << pd.d_linenumber 
	<< " column " << pd.d_cursor
        << " : " << message << std::endl << std::flush;

   pout << pd.d_linebuffer << std::endl << std::flush;

   for(int i=0; i < pd.d_cursor; i++)
      pout << " ";
   pout << "^\n";

   d_errors++;
}

/*
*************************************************************************
*									*
* Print out warnings to pout and track the number of warnings.		*
*									*
*************************************************************************
*/

void Parser::warning(const std::string& message)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();

   pout << "Warning in " << pd.d_filename << " at line " << pd.d_linenumber 
	<< " column " << pd.d_cursor
        << " : " << message << std::endl << std::flush;

   pout << pd.d_linebuffer << std::endl << std::flush;

   for(int i=0; i < pd.d_cursor; i++)
      pout << " ";
   pout << "^\n";

   d_warnings++;
}

/*
*************************************************************************
*									*
* Set the input line which is currently being parsed.    		*
*									*
*************************************************************************
*/

void Parser::setLine(const std::string& line)
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();
   pd.d_linebuffer            = line; 
}

/*
*************************************************************************
*									*
* Iterate through the database scopes, looking for the first match on	*
* the key value.							*
*									*
*************************************************************************
*/

Pointer<Database> Parser::getDatabaseWithKey(const std::string& key)
{
   List< Pointer<Database> >::Iterator i(d_scope_stack);
   for ( ; i; i++) {
      if (i()->keyExists(key)) return(i());
   }
   return(NULL);
}

/*
*************************************************************************
*									*
* Create a new parse state on the parse stack and open the specified	*
* new file for reading.							*
*									*
*************************************************************************
*/

bool Parser::pushIncludeFile(const std::string& filename)
{
   FILE *fstream = NULL;

   std::string filename_with_path;

   // If this is not a fully qualified pathname use 
   // current search path
   std::string::size_type slash_pos;
   slash_pos = filename.find_first_of( '/' );
   if ( slash_pos == 0 ) {
      filename_with_path = filename;
   } else {
      filename_with_path = d_pathname;
      filename_with_path += filename;
   }

   if (SAMRAI_MPI::getRank() == 0) {
      fstream = fopen(filename_with_path.c_str(), "r");
   }

   int worked = (fstream ? 1 : 0);

#ifdef HAVE_MPI
   worked = SAMRAI_MPI::bcast(worked, 0);
#endif

   if (!worked) {
      error("Could not open include file ``" + filename_with_path + "''");
   } else {
      ParseData pd;
      pd.d_filename   = filename_with_path;
      pd.d_fstream    = fstream;
      pd.d_linenumber = 1;
      pd.d_cursor     = 1;
      pd.d_nextcursor = 1;
      d_parse_stack.addItem(pd);
   }

   return(worked ? true : false);
}

/*
*************************************************************************
*									*
* Close the current input file and pop the parse stack.			*
*									*
*************************************************************************
*/

void Parser::popIncludeFile()
{
   Parser::ParseData& pd = d_parse_stack.getFirstItem();
   if (pd.d_fstream) fclose(pd.d_fstream);
   d_parse_stack.removeFirstItem();
}

/*
*************************************************************************
*									*
* Manage the input reading for the flex scanner.  If running with MPI,	*
* the node zero reads the data and broadcasts the length and the data	*
* to all processors.							*
*									*
*************************************************************************
*/

int Parser::yyinput(char *buffer, const int max_size)
{
   int byte = 0;
   if (SAMRAI_MPI::getRank() == 0) {
      byte = fread(buffer, 1, max_size, d_parse_stack.getFirstItem().d_fstream);
   }
#ifdef HAVE_MPI
   byte = SAMRAI_MPI::bcast(byte, 0);
   if (byte > 0) {
      SAMRAI_MPI::bcast(buffer, byte, 0);
   }
#endif
   return(byte);
}

}
}
