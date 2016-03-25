#!/usr/bin/perl

## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/grammer_fixup.pl $
## Package:     SAMRAI tests
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Code-generating script in inputdb package.

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    
    # substitution to replace [yynerrs,yychar,yylval] with SAMRAI_[yynerrs,yychar,yylval]
    s/yynerrs/SAMRAI_yynerrs/g;
    s/yychar/SAMRAI_yychar/g;
    s/yylval/SAMRAI_yylval/g;
    
    # These fixup some warning messages coming from insure++

    s/^# define YYDPRINTF\(Args\)$/# define YYDPRINTF(Args) do {} while (0)/;
    s/^# define YYDSYMPRINT\(Args\)$/# define YYDSYMPRINT(Args) do {} while (0)/;
    s/^# define YYDSYMPRINTF\(Title, Token, Value, Location\)$/# define YYDSYMPRINTF(Title, Token, Value, Location) do {} while (0)/;
    s/^# define YY_STACK_PRINT\(Bottom, Top\)$/# define YY_STACK_PRINT(Bottom, Top) do {} while (0)/;
    s/^# define YY_REDUCE_PRINT\(Rule\)$/# define YY_REDUCE_PRINT(Rule) do {} while (0)/;

    s/(\s+);}/$1\}/;

    # a more silent null use
    s/\(void\) yyvaluep;/if(0) {char *temp = (char *)&yyvaluep; temp++;}/;

    s/\/\* empty \*\/;//;

    print $_;
}

