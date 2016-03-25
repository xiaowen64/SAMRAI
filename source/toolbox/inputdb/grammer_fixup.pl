#!/usr/bin/perl

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    
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

