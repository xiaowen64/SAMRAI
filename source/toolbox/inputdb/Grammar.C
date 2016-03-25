/* A Bison parser, made by GNU Bison 1.875c.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_AND = 258,
     T_ASSIGN = 259,
     T_CHAR = 260,
     T_CLOSE_CURLY = 261,
     T_CLOSE_PAREN = 262,
     T_CLOSE_SQBKT = 263,
     T_COMMA = 264,
     T_DIV = 265,
     T_DOUBLE = 266,
     T_ELSE = 267,
     T_EXP = 268,
     T_EQUALS = 269,
     T_GREATER_EQUALS = 270,
     T_GREATER = 271,
     T_LESS_EQUALS = 272,
     T_LESS = 273,
     T_FALSE = 274,
     T_INTEGER = 275,
     T_KEYWORD = 276,
     T_MINUS = 277,
     T_MULT = 278,
     T_NOT = 279,
     T_NOT_EQUALS = 280,
     T_OR = 281,
     T_OPEN_CURLY = 282,
     T_OPEN_PAREN = 283,
     T_OPEN_SQBKT = 284,
     T_PLUS = 285,
     T_QUESTION = 286,
     T_SEMI = 287,
     T_STRING = 288,
     T_TRUE = 289,
     T_NEGATION = 290
   };
#endif
#define T_AND 258
#define T_ASSIGN 259
#define T_CHAR 260
#define T_CLOSE_CURLY 261
#define T_CLOSE_PAREN 262
#define T_CLOSE_SQBKT 263
#define T_COMMA 264
#define T_DIV 265
#define T_DOUBLE 266
#define T_ELSE 267
#define T_EXP 268
#define T_EQUALS 269
#define T_GREATER_EQUALS 270
#define T_GREATER 271
#define T_LESS_EQUALS 272
#define T_LESS 273
#define T_FALSE 274
#define T_INTEGER 275
#define T_KEYWORD 276
#define T_MINUS 277
#define T_MULT 278
#define T_NOT 279
#define T_NOT_EQUALS 280
#define T_OR 281
#define T_OPEN_CURLY 282
#define T_OPEN_PAREN 283
#define T_OPEN_SQBKT 284
#define T_PLUS 285
#define T_QUESTION 286
#define T_SEMI 287
#define T_STRING 288
#define T_TRUE 289
#define T_NEGATION 290




/* Copy the first part of user declarations.  */


//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/inputdb/Grammar.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC


// Description:	Yacc grammar description for the input database
//

#include "SAMRAI_config.h"
#include <math.h>

#include STL_SSTREAM_HEADER_FILE

using namespace std;

#if !defined(OSTRINGSTREAM_TYPE_IS_BROKEN) && defined(OSTRSTREAM_TYPE_IS_BROKEN)
typedef ostringstream ostrstream;
#endif

#include "tbox/Array.h"
#include "tbox/Complex.h"
#include "tbox/Database.h"
#include "tbox/Parser.h"
#include <string>
using namespace std;

using namespace SAMRAI;
using namespace tbox;


#ifndef NULL
#define NULL (0)
#endif

extern int yylex();
void yyerror(const char *const error)
{
   Parser::getParser()->error(error);

}

// Do not change the numbering of keys without checking promotion logic

#define KEY_COMPLEX (0)
#define KEY_DOUBLE  (1)
#define KEY_INTEGER (2)
#define KEY_BOOL    (3)
#define KEY_BOX     (4)
#define KEY_CHAR    (5)
#define KEY_STRING  (6)

static string type_names[] = {
   "complex", "double", "int", "bool", "box", "char", "string"
};

#define IS_NUMBER(X) (((X) >= 0) && ((X) < KEY_BOOL))
#define PROMOTE(X,Y) ((X) < (Y) ? (X) : (Y))

struct KeyData
{
   int             d_node_type;	// KEYDATA node type (see defines above)
   int             d_array_type;// array type (numbers may be promoted)
   int             d_array_size;// total size of the array if head element
   KeyData*   d_next;	// pointer to next (key,data) pair
   bool            d_bool;	// boolean if node is KEY_BOOL
   DatabaseBox     d_box;	// box if node is KEY_BOX
   char            d_char;	// character if node is KEY_CHAR
   dcomplex        d_complex;	// complex if node is KEY_COMPLEX
   double          d_double;	// double if node is KEY_DOUBLE
   int             d_integer;	// integer if node is KEY_INTEGER
   string          d_string;	// string if node is KEY_STRING
};

static void delete_list(KeyData*);
static void to_boolean(KeyData*);
static void to_integer(KeyData*);
static void to_double(KeyData*);
static void to_complex(KeyData*);
static KeyData* binary_op(KeyData*, KeyData*, const int);
static KeyData* compare_op(KeyData*, KeyData*, const int);
static KeyData* eval_function(KeyData*, const string&);
static KeyData* lookup_variable(const string&, const int, const bool);



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)

typedef union YYSTYPE {
  char          u_char;
  double        u_double;
  int           u_integer;
  KeyData* u_keydata;
  string*       u_keyword;
  string*       u_string;
} YYSTYPE;
/* Line 191 of yacc.c.  */

# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */


#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   define YYSTACK_ALLOC alloca
#  endif
# else
#  if defined (alloca) || defined (_ALLOCA_H)
#   define YYSTACK_ALLOC alloca
#  else
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do {  } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   249

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  36
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  13
/* YYNRULES -- Number of rules. */
#define YYNRULES  44
/* YYNRULES -- Number of states. */
#define YYNSTATES  83

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   290

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     5,     6,     9,    10,    16,    17,    22,
      24,    26,    30,    32,    34,    39,    43,    48,    49,    58,
      61,    65,    69,    73,    77,    81,    85,    89,    93,    97,
     101,   105,   109,   113,   116,   118,   120,   122,   124,   126,
     128,   130,   132,   138,   144
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      37,     0,    -1,    38,    -1,    -1,    38,    39,    -1,    -1,
      21,    27,    40,    38,     6,    -1,    -1,    21,     4,    41,
      42,    -1,    32,    -1,    43,    -1,    42,     9,    43,    -1,
      45,    -1,    21,    -1,    21,    29,    43,     8,    -1,    28,
      43,     7,    -1,    21,    28,    43,     7,    -1,    -1,    28,
      43,    31,    44,    43,    12,    43,     7,    -1,    24,    43,
      -1,    43,    26,    43,    -1,    43,     3,    43,    -1,    43,
      14,    43,    -1,    43,    25,    43,    -1,    43,    15,    43,
      -1,    43,    16,    43,    -1,    43,    17,    43,    -1,    43,
      18,    43,    -1,    43,    30,    43,    -1,    43,    22,    43,
      -1,    43,    23,    43,    -1,    43,    10,    43,    -1,    43,
      13,    43,    -1,    22,    43,    -1,    34,    -1,    19,    -1,
      47,    -1,     5,    -1,    46,    -1,    11,    -1,    20,    -1,
      33,    -1,    28,    43,     9,    43,     7,    -1,    29,    48,
       9,    48,     8,    -1,    28,    42,     7,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   155,   155,   162,   164,   174,   174,   193,   193,   278,
     294,   297,   373,   376,   380,   386,   389,   393,   393,   407,
     412,   419,   426,   429,   433,   437,   440,   444,   447,   458,
     461,   464,   467,   470,   494,   502,   510,   513,   521,   524,
     532,   540,   556,   574,   609
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "T_AND", "T_ASSIGN", "T_CHAR",
  "T_CLOSE_CURLY", "T_CLOSE_PAREN", "T_CLOSE_SQBKT", "T_COMMA", "T_DIV",
  "T_DOUBLE", "T_ELSE", "T_EXP", "T_EQUALS", "T_GREATER_EQUALS",
  "T_GREATER", "T_LESS_EQUALS", "T_LESS", "T_FALSE", "T_INTEGER",
  "T_KEYWORD", "T_MINUS", "T_MULT", "T_NOT", "T_NOT_EQUALS", "T_OR",
  "T_OPEN_CURLY", "T_OPEN_PAREN", "T_OPEN_SQBKT", "T_PLUS", "T_QUESTION",
  "T_SEMI", "T_STRING", "T_TRUE", "T_NEGATION", "$accept",
  "P_SPECIFICATION", "P_DEFINITION_LIST", "P_DEFINITION", "@1", "@2",
  "P_EXPRESSION_LIST", "P_EXPRESSION", "@3", "P_PRIMITIVE_TYPE",
  "P_COMPLEX", "P_BOX", "P_INTEGER_VECTOR", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    36,    37,    38,    38,    40,    39,    41,    39,    39,
      42,    42,    43,    43,    43,    43,    43,    44,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,    43,    43,    43,    45,    45,    45,    45,    45,    45,
      45,    45,    46,    47,    48
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     1,     0,     2,     0,     5,     0,     4,     1,
       1,     3,     1,     1,     4,     3,     4,     0,     8,     2,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     5,     5,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       3,     0,     2,     1,     0,     9,     4,     7,     5,     0,
       3,    37,    39,    35,    40,    13,     0,     0,     0,     0,
      41,    34,     8,    10,    12,    38,    36,     0,     0,     0,
      33,    19,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     6,
       0,     0,    15,     0,    17,     0,     0,    11,    21,    31,
      32,    22,    24,    25,    26,    27,    29,    30,    23,    20,
      28,    16,    14,     0,     0,    44,     0,    42,     0,    43,
       0,     0,    18
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,     1,     2,     6,    10,     9,    22,    23,    74,    24,
      25,    26,    34
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -21
static const short yypact[] =
{
     -21,    10,   -17,   -21,     7,   -21,   -21,   -21,   -21,    28,
     -21,   -21,   -21,   -21,   -21,   -20,    28,    28,    28,   -14,
     -21,   -21,     8,   183,   -21,   -21,   -21,    69,    28,    28,
     -21,   219,    56,    28,     9,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28,   -21,
      81,   102,   -21,    28,   -21,    -2,   -14,   183,   219,    22,
      22,    70,    70,    70,    70,    70,    -7,    22,    70,   201,
      -7,   -21,   -21,   123,    28,   -21,    30,   -21,   144,   -21,
      28,   165,   -21
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -21,   -21,    26,   -21,   -21,   -21,    11,   -16,   -21,   -21,
     -21,   -21,   -15
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      30,    31,    32,    37,     4,    75,    38,    35,    28,    29,
       3,     7,    50,    51,    33,     5,    45,    35,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    11,     8,    38,    27,    73,    79,    12,
       0,    76,     0,     0,    55,     0,     0,    13,    14,    15,
      16,     0,    17,     0,     0,     0,    18,    19,    78,    36,
       0,    20,    21,    52,    81,    53,    37,     0,     0,    38,
      39,    40,    41,    42,    43,    49,     0,     0,    44,    45,
      37,    46,    47,    38,    36,     0,    48,    54,    71,     0,
       4,    37,    44,    45,    38,    39,    40,    41,    42,    43,
      48,     5,     0,    44,    45,    36,    46,    47,     0,     0,
      72,    48,    37,     0,     0,    38,    39,    40,    41,    42,
      43,     0,     0,     0,    44,    45,    36,    46,    47,     0,
      77,     0,    48,    37,     0,     0,    38,    39,    40,    41,
      42,    43,     0,     0,     0,    44,    45,    36,    46,    47,
       0,     0,     0,    48,    37,     0,    80,    38,    39,    40,
      41,    42,    43,     0,     0,     0,    44,    45,    36,    46,
      47,     0,    82,     0,    48,    37,     0,     0,    38,    39,
      40,    41,    42,    43,     0,     0,    36,    44,    45,     0,
      46,    47,     0,    37,     0,    48,    38,    39,    40,    41,
      42,    43,     0,     0,    36,    44,    45,     0,    46,    47,
       0,    37,     0,    48,    38,    39,    40,    41,    42,    43,
       0,     0,     0,    44,    45,     0,    46,     0,     0,    37,
       0,    48,    38,    39,    40,    41,    42,    43,     0,     0,
       0,    44,    45,     0,    46,     0,     0,     0,     0,    48
};

static const yysigned_char yycheck[] =
{
      16,    17,    18,    10,    21,     7,    13,     9,    28,    29,
       0,     4,    28,    29,    28,    32,    23,     9,     9,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,     5,    27,    13,    10,    53,     8,    11,
      -1,    56,    -1,    -1,    33,    -1,    -1,    19,    20,    21,
      22,    -1,    24,    -1,    -1,    -1,    28,    29,    74,     3,
      -1,    33,    34,     7,    80,     9,    10,    -1,    -1,    13,
      14,    15,    16,    17,    18,     6,    -1,    -1,    22,    23,
      10,    25,    26,    13,     3,    -1,    30,    31,     7,    -1,
      21,    10,    22,    23,    13,    14,    15,    16,    17,    18,
      30,    32,    -1,    22,    23,     3,    25,    26,    -1,    -1,
       8,    30,    10,    -1,    -1,    13,    14,    15,    16,    17,
      18,    -1,    -1,    -1,    22,    23,     3,    25,    26,    -1,
       7,    -1,    30,    10,    -1,    -1,    13,    14,    15,    16,
      17,    18,    -1,    -1,    -1,    22,    23,     3,    25,    26,
      -1,    -1,    -1,    30,    10,    -1,    12,    13,    14,    15,
      16,    17,    18,    -1,    -1,    -1,    22,    23,     3,    25,
      26,    -1,     7,    -1,    30,    10,    -1,    -1,    13,    14,
      15,    16,    17,    18,    -1,    -1,     3,    22,    23,    -1,
      25,    26,    -1,    10,    -1,    30,    13,    14,    15,    16,
      17,    18,    -1,    -1,     3,    22,    23,    -1,    25,    26,
      -1,    10,    -1,    30,    13,    14,    15,    16,    17,    18,
      -1,    -1,    -1,    22,    23,    -1,    25,    -1,    -1,    10,
      -1,    30,    13,    14,    15,    16,    17,    18,    -1,    -1,
      -1,    22,    23,    -1,    25,    -1,    -1,    -1,    -1,    30
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    37,    38,     0,    21,    32,    39,     4,    27,    41,
      40,     5,    11,    19,    20,    21,    22,    24,    28,    29,
      33,    34,    42,    43,    45,    46,    47,    38,    28,    29,
      43,    43,    43,    28,    48,     9,     3,    10,    13,    14,
      15,    16,    17,    18,    22,    23,    25,    26,    30,     6,
      43,    43,     7,     9,    31,    42,     9,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,     7,     8,    43,    44,     7,    48,     7,    43,     8,
      12,    43,     7
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(SAMRAI_yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (SAMRAI_yychar == YYEMPTY && yylen == 1)				\
    {								\
      SAMRAI_yychar = (Token);						\
      SAMRAI_yylval = (Value);						\
      yytoken = YYTRANSLATE (SAMRAI_yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)		\
   ((Current).first_line   = (Rhs)[1].first_line,	\
    (Current).first_column = (Rhs)[1].first_column,	\
    (Current).last_line    = (Rhs)[N].last_line,	\
    (Current).last_column  = (Rhs)[N].last_column)
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) do {} while (0)
# define YYDSYMPRINT(Args) do {} while (0)
# define YYDSYMPRINTF(Title, Token, Value, Location) do {} while (0)
# define YY_STACK_PRINT(Bottom, Top) do {} while (0)
# define YY_REDUCE_PRINT(Rule) do {} while (0)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if defined (YYMAXDEPTH) && YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  if(0) {char *temp = (char *)&yyvaluep; temp++;}

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  if(0) {char *temp = (char *)&yyvaluep; temp++;}

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int SAMRAI_yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE SAMRAI_yylval;

/* Number of syntax errors so far.  */
int SAMRAI_yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  SAMRAI_yynerrs = 0;
  SAMRAI_yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (SAMRAI_yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      SAMRAI_yychar = YYLEX;
    }

  if (SAMRAI_yychar <= YYEOF)
    {
      SAMRAI_yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (SAMRAI_yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &SAMRAI_yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (SAMRAI_yychar != YYEOF)
    SAMRAI_yychar = YYEMPTY;

  *++yyvsp = SAMRAI_yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 5:

    {

   /* This is a hack to make a warning message go away from
      a symbol flex defines but does not use */
   if(0) {
      goto yyerrlab1;
   }

      if (Parser::getParser()->getScope()->keyExists(*yyvsp[-1].u_keyword)) {
	 string tmp("Redefinition of key ``");
         tmp += *yyvsp[-1].u_keyword;
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
      Parser::getParser()->enterScope(*yyvsp[-1].u_keyword);
   }
    break;

  case 6:

    {
      Parser::getParser()->leaveScope();
      delete yyvsp[-4].u_keyword;
   }
    break;

  case 7:

    {
      if (Parser::getParser()->getScope()->keyExists(*yyvsp[-1].u_keyword)) {
	 string tmp("Redefinition of key ``");
         tmp += *yyvsp[-1].u_keyword;
         tmp += "''";
         Parser::getParser()->warning(tmp);
      }
   }
    break;

  case 8:

    {
      KeyData* list = yyvsp[0].u_keydata;
      const int n = list->d_array_size;

      switch (list->d_array_type) {
         case KEY_BOOL: {
            Array<bool> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_bool;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putBoolArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_BOX: {
            Array<DatabaseBox> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_box;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDatabaseBoxArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_CHAR: {
            Array<char> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_char;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putCharArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_COMPLEX: {
            Array<dcomplex> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_complex(list);
               data[i] = list->d_complex;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putComplexArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_DOUBLE: {
            Array<double> data(n);
            for (int i = n-1; i >= 0; i--) {
               to_double(list);
               data[i] = list->d_double;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putDoubleArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_INTEGER: {
            Array<int> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_integer;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putIntegerArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         case KEY_STRING: {
            Array<string> data(n);
            for (int i = n-1; i >= 0; i--) {
               data[i] = list->d_string;
               list = list->d_next;
            }
            Parser::getParser()->getScope()->putStringArray(*yyvsp[-3].u_keyword, data);
            break;
         }
         default:
            Parser::getParser()->error("Internal parser error!");
            break;
      }

      delete_list(yyvsp[0].u_keydata);
      delete yyvsp[-3].u_keyword;
   }
    break;

  case 9:

    {
      Parser::getParser()->warning(
         "Semicolon found in keyword phrase (ignored)");
   }
    break;

  case 10:

    {
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 11:

    {
      switch(yyvsp[-2].u_keydata->d_array_type) {
         case KEY_BOOL:
         case KEY_CHAR:
         case KEY_STRING:
            if (yyvsp[0].u_keydata->d_node_type != yyvsp[-2].u_keydata->d_array_type) {
               Parser::getParser()->error("Type mismatch in array");
               delete yyvsp[0].u_keydata;
               yyval.u_keydata = yyvsp[-2].u_keydata;
            } else {
               yyvsp[0].u_keydata->d_array_size = yyvsp[-2].u_keydata->d_array_size + 1;
               yyvsp[0].u_keydata->d_next       = yyvsp[-2].u_keydata;
               yyval.u_keydata               = yyvsp[0].u_keydata;
            }
            break;
         case KEY_BOX:
            if (yyvsp[0].u_keydata->d_node_type != KEY_BOX) {
               Parser::getParser()->error("Type mismatch in box array");
               delete yyvsp[0].u_keydata;
               yyval.u_keydata = yyvsp[-2].u_keydata;
            } else if (yyvsp[0].u_keydata->d_box.getDimension() != yyvsp[-2].u_keydata->d_box.getDimension()) {
               Parser::getParser()->error("Box array dimension mismatch");
               delete yyvsp[0].u_keydata;
               yyval.u_keydata = yyvsp[-2].u_keydata;
            } else {
               yyvsp[0].u_keydata->d_array_size = yyvsp[-2].u_keydata->d_array_size + 1;
               yyvsp[0].u_keydata->d_next       = yyvsp[-2].u_keydata;
               yyval.u_keydata               = yyvsp[0].u_keydata;
            }
            break;
         case KEY_COMPLEX:
         case KEY_DOUBLE:
         case KEY_INTEGER:
            if (!IS_NUMBER(yyvsp[0].u_keydata->d_node_type)) {
               Parser::getParser()->error("Type mismatch in number array");
               delete yyvsp[0].u_keydata;
               yyval.u_keydata = yyvsp[-2].u_keydata;
            } else {
               yyvsp[0].u_keydata->d_array_type = PROMOTE(yyvsp[-2].u_keydata->d_array_type, yyvsp[0].u_keydata->d_node_type);
               yyvsp[0].u_keydata->d_array_size = yyvsp[-2].u_keydata->d_array_size + 1;
               yyvsp[0].u_keydata->d_next       = yyvsp[-2].u_keydata;
               yyval.u_keydata               = yyvsp[0].u_keydata;
            }
            break;
      }
   }
    break;

  case 12:

    {
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 13:

    {
      yyval.u_keydata = lookup_variable(*yyvsp[0].u_keyword, 0, false);
      delete yyvsp[0].u_keyword;
   }
    break;

  case 14:

    {
      to_integer(yyvsp[-1].u_keydata);
      yyval.u_keydata = lookup_variable(*yyvsp[-3].u_keyword, yyvsp[-1].u_keydata->d_integer, true);
      delete yyvsp[-3].u_keyword;
      delete yyvsp[-1].u_keydata;
   }
    break;

  case 15:

    {
      yyval.u_keydata = yyvsp[-1].u_keydata;
   }
    break;

  case 16:

    {
      yyval.u_keydata = eval_function(yyvsp[-1].u_keydata, *yyvsp[-3].u_keyword);
      delete yyvsp[-3].u_keyword;
   }
    break;

  case 17:

    {
      if (yyvsp[-1].u_keydata->d_node_type != KEY_BOOL) {
         Parser::getParser()->error("X in (X ? Y : Z) is not a boolean");
      }
   }
    break;

  case 18:

    {
      if ((yyvsp[-6].u_keydata->d_node_type == KEY_BOOL) && (yyvsp[-6].u_keydata->d_bool)) {
         yyval.u_keydata = yyvsp[-3].u_keydata;
         delete yyvsp[-1].u_keydata;
      } else {
         yyval.u_keydata = yyvsp[-1].u_keydata;
         delete yyvsp[-3].u_keydata;
      }
      delete yyvsp[-6].u_keydata;
   }
    break;

  case 19:

    {
      to_boolean(yyvsp[0].u_keydata);
      yyvsp[0].u_keydata->d_bool = !yyvsp[0].u_keydata->d_bool;
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 20:

    {
      to_boolean(yyvsp[-2].u_keydata);
      to_boolean(yyvsp[0].u_keydata);
      yyvsp[-2].u_keydata->d_bool = yyvsp[-2].u_keydata->d_bool || yyvsp[0].u_keydata->d_bool;
      delete yyvsp[0].u_keydata;
      yyval.u_keydata = yyvsp[-2].u_keydata;
   }
    break;

  case 21:

    {
      to_boolean(yyvsp[-2].u_keydata);
      to_boolean(yyvsp[0].u_keydata);
      yyvsp[-2].u_keydata->d_bool = yyvsp[-2].u_keydata->d_bool && yyvsp[0].u_keydata->d_bool;
      delete yyvsp[0].u_keydata;
      yyval.u_keydata = yyvsp[-2].u_keydata;
   }
    break;

  case 22:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_EQUALS);
   }
    break;

  case 23:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_EQUALS);
      yyval.u_keydata->d_bool = !(yyval.u_keydata->d_bool);
   }
    break;

  case 24:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_LESS);
      yyval.u_keydata->d_bool = !(yyval.u_keydata->d_bool);
   }
    break;

  case 25:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_GREATER);
   }
    break;

  case 26:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_GREATER);
      yyval.u_keydata->d_bool = !(yyval.u_keydata->d_bool);
   }
    break;

  case 27:

    {
      yyval.u_keydata = compare_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_LESS);
   }
    break;

  case 28:

    {
      if ((yyvsp[-2].u_keydata->d_node_type == KEY_STRING) && (yyvsp[0].u_keydata->d_node_type == KEY_STRING)) {
	 string tmp(yyvsp[-2].u_keydata->d_string);
	 tmp += yyvsp[0].u_keydata->d_string;
         yyvsp[-2].u_keydata->d_string = tmp;
         delete yyvsp[0].u_keydata;
         yyval.u_keydata = yyvsp[-2].u_keydata;
      } else {
         yyval.u_keydata = binary_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_PLUS);
      }
   }
    break;

  case 29:

    {
      yyval.u_keydata = binary_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_MINUS);
   }
    break;

  case 30:

    {
      yyval.u_keydata = binary_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_MULT);
   }
    break;

  case 31:

    {
      yyval.u_keydata = binary_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_DIV);
   }
    break;

  case 32:

    {
      yyval.u_keydata = binary_op(yyvsp[-2].u_keydata, yyvsp[0].u_keydata, T_EXP);
   }
    break;

  case 33:

    {
      switch (yyvsp[0].u_keydata->d_node_type) {
         case KEY_INTEGER:
            yyvsp[0].u_keydata->d_integer = -(yyvsp[0].u_keydata->d_integer);
            break;
         case KEY_DOUBLE:
            yyvsp[0].u_keydata->d_double = -(yyvsp[0].u_keydata->d_double);
            break;
         case KEY_COMPLEX:
            yyvsp[0].u_keydata->d_complex = -(yyvsp[0].u_keydata->d_complex);
            break;
         default:
            Parser::getParser()->error("X in -X is not a number");
            break;
      }
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 34:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_BOOL;
      yyval.u_keydata->d_array_type = KEY_BOOL;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_bool       = true;
   }
    break;

  case 35:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_BOOL;
      yyval.u_keydata->d_array_type = KEY_BOOL;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_bool       = false;
   }
    break;

  case 36:

    {
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 37:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_CHAR;
      yyval.u_keydata->d_array_type = KEY_CHAR;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_char       = yyvsp[0].u_char;
   }
    break;

  case 38:

    {
      yyval.u_keydata = yyvsp[0].u_keydata;
   }
    break;

  case 39:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_DOUBLE;
      yyval.u_keydata->d_array_type = KEY_DOUBLE;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_double     = yyvsp[0].u_double;
   }
    break;

  case 40:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_INTEGER;
      yyval.u_keydata->d_array_type = KEY_INTEGER;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_integer    = yyvsp[0].u_integer;
   }
    break;

  case 41:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_STRING;
      yyval.u_keydata->d_array_type = KEY_STRING;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;
      yyval.u_keydata->d_string     = *yyvsp[0].u_string;
      delete yyvsp[0].u_string;
   }
    break;

  case 42:

    {
      to_double(yyvsp[-3].u_keydata);
      to_double(yyvsp[-1].u_keydata);
      yyvsp[-3].u_keydata->d_complex    = dcomplex(yyvsp[-3].u_keydata->d_double, yyvsp[-1].u_keydata->d_double);
      yyvsp[-3].u_keydata->d_node_type  = KEY_COMPLEX;
      yyvsp[-3].u_keydata->d_array_type = KEY_COMPLEX;
      delete yyvsp[-1].u_keydata;
      yyval.u_keydata = yyvsp[-3].u_keydata;
   }
    break;

  case 43:

    {
      yyval.u_keydata = new KeyData;
      yyval.u_keydata->d_node_type  = KEY_BOX;
      yyval.u_keydata->d_array_type = KEY_BOX;
      yyval.u_keydata->d_array_size = 1;
      yyval.u_keydata->d_next       = NULL;

      if (yyvsp[-3].u_keydata->d_array_size != yyvsp[-1].u_keydata->d_array_size) {
         Parser::getParser()->error("Box lower/upper dimension mismatch");
      } else if (yyvsp[-3].u_keydata->d_array_size > 3) {
         Parser::getParser()->error("Box dimension too large (> 3)");
      } else {
         const int n = yyvsp[-3].u_keydata->d_array_size;
         yyval.u_keydata->d_box.setDimension(n);

         KeyData* list_lower = yyvsp[-3].u_keydata;
         KeyData* list_upper = yyvsp[-1].u_keydata;
         for (int i = n-1; i >= 0; i--) {
            yyval.u_keydata->d_box.lower(i) = list_lower->d_integer;
            yyval.u_keydata->d_box.upper(i) = list_upper->d_integer;
            list_lower = list_lower->d_next;
            list_upper = list_upper->d_next;
         }

         delete_list(yyvsp[-3].u_keydata);
         delete_list(yyvsp[-1].u_keydata);
      }
   }
    break;

  case 44:

    {
      KeyData* list = yyvsp[-1].u_keydata;
      while (list) {
         to_integer(list);
         list = list->d_next;
      }
      yyval.u_keydata = yyvsp[-1].u_keydata;
   }
    break;


    }

/* Line 1000 of yacc.c.  */


  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++SAMRAI_yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (SAMRAI_yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (SAMRAI_yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (SAMRAI_yychar == YYEOF)
	     for (;;)
	       {
		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
		 yydestruct (yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  YYDSYMPRINTF ("Error: discarding", yytoken, &SAMRAI_yylval, &yylloc);
	  yydestruct (yytoken, &SAMRAI_yylval);
	  SAMRAI_yychar = YYEMPTY;

	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

  yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = SAMRAI_yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}





/*
 * Delete all elements in a keyword list.
 */

static void delete_list(KeyData* list)
{
   while (list) {
      KeyData* byebye = list;
      list = list->d_next;
      delete byebye;
   }
}

/*
 * Verify that the number is a boolean; otherwise, report an error and
 * convert the argument into a boolean false.
 */

static void to_boolean(KeyData* keydata)
{
   if (keydata->d_node_type != KEY_BOOL) {
      Parser::getParser()->error("Cannot convert type into boolean");
      keydata->d_bool       = false;
      keydata->d_node_type  = KEY_BOOL;
      keydata->d_array_type = KEY_BOOL;
   }
}

/*
 * Convert the number into an integer.  If the conversion cannot be
 * performed, then print an error and return an integer zero.
 */

static void to_integer(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         break;
      case KEY_DOUBLE:
         Parser::getParser()->warning("Double truncated to integer");
         keydata->d_integer = (int)(keydata->d_double);
         break;
      case KEY_COMPLEX:
         Parser::getParser()->warning("Complex truncated to integer");
         keydata->d_integer = (int)(keydata->d_complex.real());
         break;
      default:
         Parser::getParser()->error("Cannot convert type into integer");
         keydata->d_integer = 0;
         break;
   }
   keydata->d_node_type  = KEY_INTEGER;
   keydata->d_array_type = KEY_INTEGER;
}

/*
 * Convert the number in the keydata structure to a double.  If the
 * conversion cannot be performed, then print an error and return zero.
 */

static void to_double(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         keydata->d_double = (double)(keydata->d_integer);
         break;
      case KEY_DOUBLE:
         break;
      case KEY_COMPLEX:
         Parser::getParser()->warning("Complex truncated to double");
         keydata->d_double = keydata->d_complex.real();
         break;
      default:
         Parser::getParser()->error("Cannot convert type into double");
         keydata->d_double = 0.0;
         break;
   }
   keydata->d_node_type  = KEY_DOUBLE;
   keydata->d_array_type = KEY_DOUBLE;
}

/*
 * Convert the number in the keydata structure to a complex.  If the
 * conversion cannot be performed, then print an error and return zero.
 */

static void to_complex(KeyData* keydata)
{
   switch (keydata->d_node_type) {
      case KEY_INTEGER:
         keydata->d_complex = dcomplex((double) keydata->d_integer, 0.0);
         break;
      case KEY_DOUBLE:
         keydata->d_complex = dcomplex(keydata->d_double, 0.0);
         break;
      case KEY_COMPLEX:
         break;
      default:
         Parser::getParser()->error("Cannot convert type into complex");
         keydata->d_complex = dcomplex(0.0, 0.0);
         break;
   }
   keydata->d_node_type  = KEY_COMPLEX;
   keydata->d_array_type = KEY_COMPLEX;
}

/*
 * Perform one of the standard binary operations +, -, *, /, or ^ on numeric
 * types.  Concatenation for strings is implemented above.  Return an integer
 * zero if there is a type mismatch.
 */

static KeyData* binary_op(KeyData* a, KeyData* b, const int op)
{
   if (!IS_NUMBER(a->d_node_type) || !IS_NUMBER(b->d_node_type)) {
      Parser::getParser()->error(
         "Cannot perform numerical operations on non-numeric types");
      a->d_integer    = 0;
      a->d_node_type  = KEY_INTEGER;
      a->d_array_type = KEY_INTEGER;
   } else {
      const int result_type = PROMOTE(a->d_node_type, b->d_node_type);
      switch (result_type) {
         case KEY_INTEGER:
            switch (op) {
               case T_DIV:
                  a->d_integer = a->d_integer / b->d_integer;
                  break;
               case T_EXP:
                  a->d_integer =
                     (int) pow((double) a->d_integer, (double) b->d_integer);
                  break;
               case T_MINUS:
                  a->d_integer = a->d_integer - b->d_integer;
                  break;
               case T_MULT:
                  a->d_integer = a->d_integer * b->d_integer;
                  break;
               case T_PLUS:
                  a->d_integer = a->d_integer + b->d_integer;
                  break;
            }
            break;
         case KEY_DOUBLE:
            to_double(a);
            to_double(b);
            switch (op) {
               case T_DIV:
                  a->d_double = a->d_double / b->d_double;
                  break;
               case T_EXP:
                  a->d_double = pow(a->d_double, b->d_double);
                  break;
               case T_MINUS:
                  a->d_double = a->d_double - b->d_double;
                  break;
               case T_MULT:
                  a->d_double = a->d_double * b->d_double;
                  break;
               case T_PLUS:
                  a->d_double = a->d_double + b->d_double;
                  break;
            }
            break;
         case KEY_COMPLEX:
            to_complex(a);
            to_complex(b);
            switch (op) {
               case T_DIV:
                  a->d_complex = a->d_complex / b->d_complex;
                  break;
               case T_EXP:
                  /*
		   * SGS this is broken for insure++ and gcc 3.3.2
		   * a->d_complex = pow(a->d_complex, b->d_complex);
		   * replaced with the defn from the header file.
		   */
		  a->d_complex = exp(a->d_complex * log(b->d_complex));
                  break;
               case T_MINUS:
                  a->d_complex = a->d_complex - b->d_complex;
                  break;
               case T_MULT:
                  a->d_complex = a->d_complex * b->d_complex;
                  break;
               case T_PLUS:
                  a->d_complex = a->d_complex + b->d_complex;
                  break;
            }
            break;
      }
   }
   delete b;
   return(a);
}

/*
 * Perform one of the standard comparison operations ==, <, or >.  The other
 * operators !=, >=, and <= are computed above by using one of the first three
 * and then negating the result.  Return a boolean false if there is a type
 * mismatch problem.
 */

static KeyData* compare_op(KeyData* a, KeyData* b, const int op)
{
   if (!IS_NUMBER(a->d_node_type) || !IS_NUMBER(b->d_node_type)) {
      if (a->d_node_type != b->d_node_type) {
         Parser::getParser()->error(
            "Cannot compare different non-numeric types");
         a->d_bool = false;
      } else if (op != T_EQUALS) {
         Parser::getParser()->error(
            "Cannot apply <, >, <=, or >= to non-numeric types");
         a->d_bool = false;
      } else {
         switch(a->d_node_type) {
            case KEY_BOOL:
               a->d_bool = (a->d_bool == b->d_bool);
               break;
            case KEY_BOX:
               a->d_bool = (a->d_box == b->d_box);
               break;
            case KEY_CHAR:
               a->d_bool = (a->d_char == b->d_char);
               break;
            case KEY_STRING:
               a->d_bool = (a->d_string == b->d_string);
               break;
         }
      }
   } else {
      const int promoted = PROMOTE(a->d_node_type, b->d_node_type);
      switch (promoted) {
         case KEY_INTEGER:
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_integer == b->d_integer);
                  break;
               case T_LESS:
                  a->d_bool = (a->d_integer < b->d_integer);
                  break;
               case T_GREATER:
                  a->d_bool = (a->d_integer > b->d_integer);
                  break;
            }
            break;
         case KEY_DOUBLE:
// Intel warns about comparison of floating point numbers
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
            to_double(a);
            to_double(b);
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_double == b->d_double);
                  break;
               case T_LESS:
                  a->d_bool = (a->d_double < b->d_double);
                  break;
               case T_GREATER:
                  a->d_bool = (a->d_double > b->d_double);
                  break;
            }
            break;
         case KEY_COMPLEX:
            to_complex(a);
            to_complex(b);
            switch (op) {
               case T_EQUALS:
                  a->d_bool = (a->d_complex == b->d_complex);
                  break;
               case T_LESS:
               case T_GREATER:
                  Parser::getParser()->error(
                     "Operators <, >, <=, and >= are not defined for complex");
                  a->d_bool = false;
                  break;
            }
            break;
      }
   }
   a->d_node_type  = KEY_BOOL;
   a->d_array_type = KEY_BOOL;
   delete b;
   return(a);
}

/*
 * Perform a function evaluation on the specified argument.
 */

struct arith_functions {
   string     d_name;
   double   (*d_r2r_func)(double);
   dcomplex (*d_c2c_func)(const dcomplex&);
   double   (*d_c2r_func)(const dcomplex&);
};



#if 0
// Static initialization was not working with SGI
// compiler; so use an initialization function to 
// create the table
static arith_functions af[] = {
   { "abs"  , fabs , NULL, abs  },
   { "acos" , acos , NULL, NULL },
   { "asin" , asin , NULL, NULL },
   { "atan" , atan , NULL, NULL },
   { "ceil" , ceil , NULL, NULL },
   { "conj" , NULL , conj, NULL },
   { "cos"  , cos  , cos , NULL },
   { "cosh" , cosh , cosh, NULL },
   { "exp"  , exp  , exp , NULL },
   { "fabs" , fabs , NULL, NULL },
   { "floor", floor, NULL, NULL },
   { "imag" , NULL , NULL, imag },
   { "log10", log10, NULL, NULL },
   { "log"  , log  , log , NULL },
   { "real" , NULL , NULL, real },
   { "sin"  , sin  , sin , NULL },
   { "sinh" , sinh , sinh, NULL },
   { "sqrt" , sqrt , sqrt, NULL },
   { "tan"  , tan  , NULL, NULL },
   { ""     , NULL , NULL, NULL }
};
#endif

static arith_functions af[20];

// These are needed to deal with imag/real returning a reference
// under GCC 3.4.x
static double imag_thunk(const dcomplex &a)
{
   return imag(a);
}

static double real_thunk(const dcomplex &a)
{
   return real(a);
}

void parser_static_table_initialize()
{
   af[0].d_name =    "abs";
   af[0].d_r2r_func = fabs;
   af[0].d_c2c_func = NULL;
   af[0].d_c2r_func = std::abs;


   af[1].d_name =    "acos";
   af[1].d_r2r_func = acos;
   af[1].d_c2c_func = NULL;
   af[1].d_c2r_func = NULL;
   
   af[2].d_name =    "asin";
   af[2].d_r2r_func = asin;
   af[2].d_c2c_func = NULL;
   af[2].d_c2r_func = NULL;
   
   af[3].d_name =    "atan";
   af[3].d_r2r_func = atan;
   af[3].d_c2c_func = NULL;
   af[3].d_c2r_func = NULL;
   
   af[4].d_name =    "ceil";
   af[4].d_r2r_func = ceil;
   af[4].d_c2c_func = NULL;
   af[4].d_c2r_func = NULL;

   af[5].d_name =    "conj";
   af[5].d_r2r_func = NULL;
   af[5].d_c2c_func = conj;
   af[5].d_c2r_func = NULL;


   af[6].d_name =    "cos";
   af[6].d_r2r_func = ::cos;
   af[6].d_c2c_func = std::cos;
   af[6].d_c2r_func = NULL;

   af[7].d_name =    "cosh";
   af[7].d_r2r_func = ::cosh;
   af[7].d_c2c_func = std::cosh;
   af[7].d_c2r_func = NULL;

   af[8].d_name =    "exp";
   af[8].d_r2r_func = ::exp;
   af[8].d_c2c_func = std::exp;
   af[8].d_c2r_func = NULL;

   af[9].d_name =    "fabs";
   af[9].d_r2r_func = fabs;
   af[9].d_c2c_func = NULL;
   af[9].d_c2r_func = NULL;

   af[10].d_name =    "floor";
   af[10].d_r2r_func = floor;
   af[10].d_c2c_func = NULL;
   af[10].d_c2r_func = NULL;

   af[11].d_name =    "imag";
   af[11].d_r2r_func = NULL;
   af[11].d_c2c_func = NULL;
   af[11].d_c2r_func = imag_thunk;

   af[12].d_name =    "log10";
   af[12].d_r2r_func = ::log10;
   af[12].d_c2c_func = NULL;
   af[12].d_c2r_func = NULL;

   af[13].d_name =    "log";
   af[13].d_r2r_func = ::log;
   af[13].d_c2c_func = std::log;
   af[13].d_c2r_func = NULL;

   af[14].d_name =    "real";
   af[14].d_r2r_func = NULL;
   af[14].d_c2c_func = NULL;
   af[14].d_c2r_func = real_thunk;

   af[15].d_name =    "sin";
   af[15].d_r2r_func = ::sin;
   af[15].d_c2c_func = std::sin;
   af[15].d_c2r_func = NULL;

   af[16].d_name =    "sinh";
   af[16].d_r2r_func = ::sinh;
   af[16].d_c2c_func = std::sinh;
   af[16].d_c2r_func = NULL;

   af[17].d_name =    "sqrt";
   af[17].d_r2r_func = ::sqrt;
   af[17].d_c2c_func = std::sqrt;
   af[17].d_c2r_func = NULL;

   af[18].d_name =    "tan";
   af[18].d_r2r_func = tan;
   af[18].d_c2c_func = NULL;
   af[18].d_c2r_func = NULL;

   af[19].d_name =    "";
   af[19].d_r2r_func = NULL;
   af[19].d_c2c_func = NULL;
   af[19].d_c2r_func = NULL;
}

static KeyData* eval_function(KeyData* arg, const string& func)
{
   if (!IS_NUMBER(arg->d_node_type)) {
      string tmp("Unknown function ");
      tmp += func;
      tmp += "(";
      tmp += type_names[arg->d_node_type];
      tmp += ")";
      Parser::getParser()->error(tmp);
   } else if (func == "int") {
      to_double(arg);
      arg->d_integer    = (int) arg->d_double;
      arg->d_node_type  = KEY_INTEGER;
      arg->d_array_type = KEY_INTEGER;
   } else {
      for (int f = 0; af[f].d_name.length() > 0; f++) {
         if (af[f].d_name == func) {
            if (arg->d_node_type == KEY_COMPLEX) {
               if (af[f].d_c2c_func) {
                  arg->d_complex = (*af[f].d_c2c_func)(arg->d_complex);
               } else if (af[f].d_c2r_func) {
                  arg->d_double     = (*af[f].d_c2r_func)(arg->d_complex);
                  arg->d_node_type  = KEY_DOUBLE;
                  arg->d_array_type = KEY_DOUBLE;
               } else {
                  to_double(arg);
                  arg->d_double = (*af[f].d_r2r_func)(arg->d_double);
               }
            } else {
               if (af[f].d_r2r_func) {
                  to_double(arg);
                  arg->d_double = (*af[f].d_r2r_func)(arg->d_double);
               } else if (af[f].d_c2r_func) {
                  to_complex(arg);
                  arg->d_double     = (*af[f].d_c2r_func)(arg->d_complex);
                  arg->d_node_type  = KEY_DOUBLE;
                  arg->d_array_type = KEY_DOUBLE;
               } else {
                  to_complex(arg);
                  arg->d_complex = (*af[f].d_c2c_func)(arg->d_complex);
               }
            }
            return(arg);
         }
      }

      string tmp("Unknown function ");
      tmp += func;
      tmp += "(";
      tmp += type_names[arg->d_node_type];
      tmp += ")";
      Parser::getParser()->error(tmp);
   }
   return(arg);
}

/*
 * Fetch a variable in the database.  If there is an error, then print
 * an error message and return an integer zero as result.
 */

static KeyData* lookup_variable(
   const string& key, const int index, const bool is_array)
{
   KeyData* result = new KeyData;
   result->d_node_type  = KEY_INTEGER;
   result->d_array_type = KEY_INTEGER;
   result->d_array_size = 1;
   result->d_next       = NULL;
   result->d_integer    = 0;

   Parser *parser = Parser::getParser();
   Pointer<Database> db = parser->getDatabaseWithKey(key);

   if (db.isNull()) {
      string tmp("Variable ``");
      tmp += key;
      tmp += "'' not found in database";
      parser->error(tmp);
   } else if (!is_array && (db->getArraySize(key) > 1)) {
      string tmp("Variable ``");
      tmp += key;
      tmp += "'' is not a scalar value";
      parser->error(tmp);
   } else if ((index < 0) || (index >= db->getArraySize(key))) {
      ostrstream oss;
      oss << index;
      string tmp("Variable ``");
      tmp += key;
      tmp += "[";
      tmp += oss.str();
      tmp += "]'' out of range";
      parser->error(tmp);
   } else if (db->isInteger(key)) {
      result->d_integer    = db->getIntegerArray(key)[index];
      result->d_node_type  = KEY_INTEGER;
      result->d_array_type = KEY_INTEGER;

   } else if (db->isDouble(key)) {
      result->d_double     = db->getDoubleArray(key)[index];
      result->d_node_type  = KEY_DOUBLE;
      result->d_array_type = KEY_DOUBLE;

   } else if (db->isComplex(key)) {
      result->d_complex    = db->getComplexArray(key)[index];
      result->d_node_type  = KEY_COMPLEX;
      result->d_array_type = KEY_COMPLEX;

   } else if (db->isBool(key)) {
      result->d_bool       = db->getBoolArray(key)[index];
      result->d_node_type  = KEY_BOOL;
      result->d_array_type = KEY_BOOL;

   } else if (db->isDatabaseBox(key)) {
      result->d_box        = db->getDatabaseBoxArray(key)[index];
      result->d_node_type  = KEY_BOX;
      result->d_array_type = KEY_BOX;

   } else if (db->isChar(key)) {
      result->d_char       = db->getCharArray(key)[index];
      result->d_node_type  = KEY_CHAR;
      result->d_array_type = KEY_CHAR;

   } else if (db->isString(key)) {
      result->d_string     = db->getStringArray(key)[index];
      result->d_node_type  = KEY_STRING;
      result->d_array_type = KEY_STRING;

   } else {
      parser->error("Unknown type for variable="+key);
   }

   return(result);
}

