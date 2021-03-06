/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  April 2015 (updated July 2016, Apr. 2017)

  Header file for stitch.c.
*/
#define VERSION     "0.2"

#define MAX_SIZE    1024   // maximum length of input lines
#define NOTMATCH    1.5f   // stitch failure
#define GZEXT       ".gz"  // file extension for gzip compression
#define COM         ", "   // separator for input file names

// command-line parameters
#define HELP        "-h"
#define HELP2       "--help"
#define FIRST       "-1"
#define SECOND      "-2"
#define OUTFILE     "-o"
#define UNFILE      "-u"
#define LOGFILE     "-l"
#define OVERLAP     "-m"
#define MISMATCH    "-p"
#define DOVEOPT     "-d"
#define DOVEOVER    "-e"
#define DOVEFILE    "-f"
#define MAXOPT      "-n"
#define ADAPTOPT    "-a"
#define VERBOSE     "-v"
#define VERSOPT     "--version"

// extensions for output files
#define ONEEXT      "_1.fastq"
#define TWOEXT      "_2.fastq"

// default parameter values
#define DEFOVER     20    // min. overlap
#define DEFDOVE     50    // min. overlap for dovetailed reads
#define DEFMISM     0.1f  // mismatch fraction
#define NA          "NA"  // n/a

// third parameter to copyStr()
#define FWD         0
#define RC          1
#define REV         2

// error messages
#define ERROPEN     0
#define MERROPEN    ": cannot open file for reading"
#define ERRCLOSE    1
#define MERRCLOSE   "Cannot close file"
#define ERROPENW    2
#define MERROPENW   ": cannot open file for writing"
#define ERRUNK      3
#define MERRUNK     ": unknown nucleotide"
#define ERRMEM      4
#define MERRMEM     "Cannot allocate memory"
#define ERRSEQ      5
#define MERRSEQ     "Cannot load sequence"
#define ERRQUAL     6
#define MERRQUAL    "Sequence/quality scores do not match"
#define ERRHEAD     7
#define MERRHEAD    ": not matched in input files"
#define ERRINT      8
#define MERRINT     ": cannot convert to int"
#define ERRFLOAT    9
#define MERRFLOAT   ": cannot convert to float"
#define ERRPARAM    10
#define MERRPARAM   ": unknown command-line parameter"
#define ERROVER     11
#define MERROVER    "Overlap must be greater than 0"
#define ERRMISM     12
#define MERRMISM    "Mismatch must be in [0,1)"
#define DEFERR      "Unknown error"

typedef union file {
  FILE* f;
  gzFile gzf;
} File;
