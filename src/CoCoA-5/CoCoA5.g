//   Copyright (c) 2009 Giovanni Lagorio (lagorio@disi.unige.it)
//
//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

grammar CoCoA5;

options {
  k=3;
}

// Notes:
// 1) The names "parse<Something>" in the following comments refer to the corresponding methods of 
//    class CoCoA::ParserNS::Parser
// 2) All keywords are case insensitive so, for instance, 'Use' should be specified as
//    (('U'|'u') ('S'|'s') ('E'|'e')).
//    For the time being, there is an exception: 'RECORD' (all uppercase)
//    is recognized as an Identifier instead of (the uppercase version of) the keyword 'Record',
//    while all the other casing (e.g., ReCoRD, recoRD, RECORd ...) are treated as the keyword 
//    Needless to say, this is a backward compatibility hack and this issue should be rectified
// 3) Since CoCoA 4 is inconsistent in the treatment of semicolons, the current parser for
//    CoCoA 5 is quite forgiving when semicolons are omitted (the extra-pedantic warning level
//    has to be turned on to get a feedback about this issue). Future versions might become strictier.

/* parseTopLevelStatement */
topLevelStatement
    : /* parseUseStatement         */ useStatement 
    | /* parsePackageStatement     */ packageDefinition
    | /* parseSourceStatement      */ sourceStatement     
    | /* parseAliasStatement       */ aliasStatement      
    | /* parseDefineStatement      */ defineStatement     
    |                                 funBodyStatement
    ;

/* parsePackageDefinition */
packageDefinition
    : 'Package' PackageName exportList packageMember* ('EndPackage'|'End')
    ;

/* parsePackageExportList */
exportList
    : ('export' Identifier (',' Identifier)* ';')*
    ;

packageMember
    : /* parseDefineStatement      */ defineStatement     
    | /* parseAliasStatement       */ aliasStatement      
    | Identifier ':=' expression
    ;

funBodyStatement
    : /* parseIfStatment               */ ifStatment          
    | /* parseReturnStatement          */ returnStatement     
    | /* parseUsing                    */ usingStatement            /* ??? not sure about this one ??? */
    | /* parsePrint                    */ printStatement      
    | /* parseBreakOrContinueStatement */ breakStatement      
    | /* parseBreakOrContinueStatement */ continueStatement      
    | /* parseDescribeStatement        */ describeStatement   
    | /* parseCiaoOrQuitStatement      */ ciaoOrQuitStatement 
    | /* parseSkipStatement            */ skipStatement       
    | /* parseBlockStatement           */ blockStatement      
    |                                     helpStatement       
    | /* parseTryStatement             */ tryStatement        
    | /* parseProtect                  */ protectStatement    
    | /* parseUnprotect                */ unprotectStatement  
    |                                     timeStatement       
    |                                     (Identifier ':')? loopStatements
    |                                     evalOrAssignmentStmt
    /* OBSOLESCENT STATEMENTS: */
    | /* parseSetStatement             */ setStatement        
    | /* parseUnsetStatement           */ unsetStatement      
    | /* parseClearStatement           */ clearStatement      
    | /* parseDeleteStatement          */ deleteStatement     
    | /* parseDestroyStatement         */ destroyStatement    
    | /* parseCatchStatement           */ catchStatement      
    ;

loopStatements
    : /* parseForStatement         */ forStatement        
    | /* parseForeachStatement     */ foreachStatement    
    | /* parseWhileStatement       */ whileStatement      
    | /* parseRepeatUntilStatement */ repeatUntilStatement
    ;

useStatement         : 'Use' (ringDefinition | Identifier '::=' ringDefinition  )';';
setStatement         : 'Set' Identifier ( ':=' expression )? ';';
unsetStatement       : 'Unset' Identifier ';';
forStatement         : 'For' Identifier ':=' expression 'To' expression ('Step' expression)? 'Do' statements ('EndFor'|'End');
foreachStatement     : 'Foreach' Identifier 'In' expression 'Do' statements ('EndForeach'|'End');
ifStatment           : 'If' expression 'Then' statements ('Elif' expression 'Then' statements)* ('Else' statements) ('Endif'|'End');
returnStatement      : 'Return' expression? ';';
usingStatement       : 'Using' Identifier 'Do' statements ('EndUsing'|'End');
printStatement       : ('Print'|'PrintLn') (expressionOrNewline (',' expressionOrNewline)*)? ('On' expression)? ';';
breakStatement       : 'Break' Identifier? ';';
continueStatement    : 'Continue' Identifier? ';';
repeatUntilStatement : 'Repeat' statements ( ('End'|'EndRepeat')|'Until' expression ';');
whileStatement       : 'While' expression 'Do' statements ('EndWhile'|'End');
describeStatement    : 'Describe' (PackageName|expression) ';';
clearStatement       : 'Clear' Identifier (',' Identifier)* ';';
deleteStatement      : 'Delete' Identifier (',' Identifier)* ';';
destroyStatement     : 'Destroy' Identifier (',' Identifier)* ';';
aliasStatement       : 'Alias' bindingList (('In' statements ('EndAlias'|'End')) |';'); /* the form Alias...In...EndAlias is obsolescent */
ciaoOrQuitStatement  : ('Ciao'|'Quit') ';';
skipStatement        : 'Skip' ';';
catchStatement       : 'Catch' statements ('In' Identifier)? ('EndCatch'|'End');
blockStatement       : 'Block' statements ('EndBlock'|'End');
helpStatement        :  HelpStatement;
sourceStatement      : ('<<'|'Source'|'Load') expression ';';
defineStatement      : 'Define' Identifier funBody ('EndDefine'|'End');
tryStatement         : 'Try' statements 'UponError' Identifier 'Do' statements ('EndTry'|'End');
protectStatement     : 'Protect' Identifier (':' expression)? ';';
unprotectStatement   : 'Unprotect' Identifier ';';
evalOrAssignmentStmt :  expression (':=' expression | '::=' ringDefinition)? ';';
timeStatement        : 'Time' funBodyStatement /* but not an expression-statement */;

/* parseStatements */
statements
    : funBodyStatement*
    ;

/* parseFunctionDeclaration */    
funBody
    : '(' (paramList|'...') ')' importStatement* statements
    ;

/* parseImports */
importStatement
    : (('ImportByRef')|('TopLevel')|('ImportByValue')) Identifier (',' Identifier)* ';'
    ;

/* parseBindingList */
bindingList
    : Identifier ':=' PackageName (',' Identifier ':=' PackageName)*
    ;
        
/* parseExpOrNewLine */
expressionOrNewline
    : expression
    // | 'NewLine' /* not implemented (yet?) */
    ;
    
/* parseIdOrVarList */
idList
    : (Identifier (',' Identifier)*)?
    ;
    
/* parseParamList */
paramList // note: all Opt parameter must be at the end
    : ( (('Ref')|('Opt'))? Identifier (',' (('Ref')|('Opt'))? Identifier)*)?
    ;
        
/* parseRingDefinition */
ringDefinition
    :  Identifier ringDefQuotient? ringDefPoly?
    ;
    
ringDefQuotient
    : '/' '(' expression ')'
    ;
    
ringDefPoly
    : '[' indeterminateDecl ( ',' indeterminateDecl )* ']'
      ( ',' ('ToPos' | 'PosTo' | ('Weights' '(' expression (',' expression)* ')')|'Lex'|'Xel'|'DegRevLex'|'DegLex'|('Ord' '(' expression ')')|
      (Identifier /* Compatibility hack: Identifier must be 'Elim' or 'Ord' */ '(' elimIndeterminateDecl ('..' elimIndeterminateDecl)? ')')))*
    ;

/* parseIndeterminateDecl(isInsideElim=false) */
indeterminateDecl    
    : Identifier ('[' additiveExpression ('..'additiveExpression)? (',' additiveExpression ('..'additiveExpression)?)* ']')?
        // here we use additiveExpression because '..' creates ambiguities
    ;

/* parseIndeterminateDecl(isInsideElim=true) */
elimIndeterminateDecl    
    : Identifier ('[' additiveExpression  (',' additiveExpression)* ']')?
        // here we use additiveExpression because '..' creates ambiguities
    ;
    
    
/* parseExpressionList */
expressionList
    : expression (',' expression)*
    ;
  
/* parseArgumentList */
argumentList
    : 'ref'? expression (',' 'ref'? expression)*
    ;
    
/* parseExpression */
expression
    : conditionalOrExpression
    ;

/* parseConditionOrExpression */
conditionalOrExpression
    : conditionalAndExpression ( 'Or' conditionalAndExpression )*
    ;

/* parseConditionalAndExpression */
conditionalAndExpression
    : equalityExpression ( 'And' equalityExpression )*
    ;

/* parseEqualityExpression */
equalityExpression
    : relationalExpression ( ('=' | '<>' | 'IsIn') relationalExpression )*
    ;

/* parseRelationalExpression */
relationalExpression
    : cartesianProductExpression ( ('<='|'<'|'>'|'>=') cartesianProductExpression )*
    ;
    
/* parseCartesianProductExpression */
cartesianProductExpression
    : listExpression ( '><' listExpression )* // note: this is an n-ary operator
        ;
  
/* parseListExpression */
listExpression
    : additiveExpression ('..' additiveExpression)?
    ;

/* parseAdditiveExpression */
additiveExpression
    : multiplicativeExpression ( ('+' | '-') multiplicativeExpression )*
    ;

/* parseMultiplicativeExpression */
multiplicativeExpression
    : unaryExpression ( ( '*' | '/' | '%' ) unaryExpression )*
    ;
    
/* parseUnaryExpression */
unaryExpression
	: '+' unaryExpression
	| '-' unaryExpression
	| 'Not' unaryExpression
	| powerExpression 
	;
	
/* parsePowerExpression */
powerExpression
    : primary selector* ('^' powerExpression)?
    ;

/* parsePrimary */
primary
    : '(' expression ')'
    | 'IsDefined' '(' Identifier ')'
    | 'True'
    | 'False'
    | (PackageName '.')? Identifier ('::' primary)?
    | FLOAT_LITERAL
    | INT_LITERAL
    | STRING_LITERAL STRING_LITERAL*
    | /* parseListPrimary */ '[' ']'
    | /* parseListPrimary */ '[' Identifier 'In' expression '|' expression ']' 
    | /* parseListPrimary */ '[' expression ('|' Identifier 'In' equalityExpression ('And' expression)? | (',' expression)*) ']'
                                // here we use equalityExpression because 'And' creates ambiguities
    | /* parseRecord       */ 'Record' '[' ((Identifier ('='|':=') expression) (',' Identifier ('='|':=') expression)*)? ']'
    | 'Ord' '(' expressionList? ')'
    | '***' expression '***'
    | 'Lambda' funBody (('End')|('EndLambda'))
    ;
     
/* parseSelectors */
selector
    :   '.' Identifier
    |   '[' expressionList ']'
    |   '(' argumentList? ')'
    |	'(' '...' ')'
    ;

// Lexer:
 
Identifier 
    : ('a'..'z'|'A'..'Z'|'_')('a'..'z'|'A'..'Z'|'_'|'0'..'9')*
    ;

PackageName
    : '$' ('/' | 'a'..'z'|'A'..'Z'|'_'|'0'..'9')*
    ;

WS
    : (' '|'\r'|'\t'|'\u000C'|'\n') {$channel=HIDDEN;}
    ;

COMMENT
    : '/*' ( options {greedy=false;} : . )* '*/' {$channel=HIDDEN;}
    ;

HelpStatement
    : '?' ~('\n'|'\r')* '\r'? '\n'
    ;
    
LINE_COMMENT
    : ('//'|'--') ~('\n'|'\r')* '\r'? '\n' {$channel=HIDDEN;}
    ;
    
INT_LITERAL
    : ('0'..'9')+
    ;
    
FLOAT_LITERAL
    : ('0'..'9')+ '.' ('0'..'9'+)
    ;
    
STRING_LITERAL // note: newlines are not allowed inside string literals (users must use the escape sequences)
    : '"'  ( ESCAPE_SEQUENCE | ~('\\'|'"' ) )* '"'
    | '\'' ( ESCAPE_SEQUENCE | ~('\\'|'\'') )* '\''
    ;

fragment
ESCAPE_SEQUENCE
    : '\\' ('a'|'t'|'n'|'r'|'\"'|'\''|'\\')
    ;

