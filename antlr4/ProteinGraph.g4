grammar ProteinGraph;

// Parser rules
// Full protein header
protein: (ACCESSION position COMMA misscleavages COMMA? RPAR) | (ACCESSION position COMMA misscleavages COMMA features RPAR);

//Parse Features
features: term_feature ( ( COMMA | AND ) term_feature)*;
term_feature: variant | mutagen | conflict | signal | init_met | propep | peptide | chain | fixmod | varmod | or;

// Single features
or: OR LBRAC ( features PIPE | PIPE | NONE PIPE ) (features PIPE | NONE PIPE | features | NONE )* RBRAC;
variant: (VARIANT LBRAC position COMMA qualifier_info RBRAC) | (VARIANT LBRAC position COMMA qualifier_info COMMA FEATURE_ID RBRAC);
mutagen: (MUTAGEN LBRAC position COMMA qualifier_info RBRAC) | (MUTAGEN LBRAC position COMMA qualifier_info COMMA FEATURE_ID RBRAC);
conflict: (CONFLICT LBRAC position COMMA qualifier_info RBRAC) | (CONFLICT LBRAC position COMMA qualifier_info COMMA FEATURE_ID RBRAC);
signal: SIGNAL LBRAC position (COMMA FEATURE_ID)? RBRAC;
init_met: INIT_MET LBRAC position (COMMA FEATURE_ID)? RBRAC;
propep: PROPEP LBRAC position (COMMA FEATURE_ID)? RBRAC;
peptide: PEPTIDE LBRAC position (COMMA FEATURE_ID)? RBRAC;
chain: CHAIN LBRAC position (COMMA FEATURE_ID)? RBRAC;
fixmod: FIXMOD LBRAC position COMMA qualifier_info RBRAC;
varmod: VARMOD LBRAC position COMMA qualifier_info RBRAC;

// "Terminal" rules
qualifier_info: replacement | shift | MISSING;
replacement: AMINOS REPL AMINOS;
shift: AMINOS DD NUMBER;
position: NUMBER DD NUMBER;
misscleavages: MSSCLVG DD NUMBER;


// Lexer rules
// Fixed Strings
VARIANT: 'VARIANT';
MUTAGEN: 'MUTAGEN';
CONFLICT: 'CONFLICT';
INIT_MET: 'INIT_MET';
SIGNAL: 'SIGNAL';
PROPEP: 'PROPEP';
PEPTIDE: 'PEPTIDE';
CHAIN: 'CHAIN';
OR: 'OR';
VARMOD: 'VARMOD';
FIXMOD: 'FIXMOD';
MISSING: 'Missing';
MSSCLVG: 'mssclvg';
COMMA: ',';
DD: ':';
REPL: '->';
LPAR: '(';  // Not needed due to no positive lookahead
RPAR: ')';
LBRAC: '[';
RBRAC: ']';
AND: '&';
PIPE: '|';
NONE: 'None';

// Regex to fit parsable info
FEATURE_ID: [A-Z]+ '_' [0-9]+;
ACCESSION: [a-zA-Z0-9-]+ '('; // There is no positive lookahead in antlr
AMINOS: [A-Z]+;
NUMBER: [-<>]?[0-9]+ | '?' | [-]?[0-9]+ '.' [0-9]+; // Special case: Numbers can be '?' --> no knowledge about position 
WS: [ \t\r\n]+ -> skip;
