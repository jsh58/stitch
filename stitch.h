/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  April 2015 (updated July 2016, Apr. 2017, June 2017)

  Header file for stitch.c.
*/
#define VERSION     "0.6"

// constants
#define MAX_SIZE    1024    // maximum length of input lines
#define NOTMATCH    1.5f    // stitch failure
#define COM         ", "    // separator for input file names
#define MINQUAL     2       // minimum quality score (Illumina)
#define MAXQUAL     41      // maximum quality score
#define NA          "NA"    // n/a

// default parameter values
#define DEFOVER     20      // min. overlap
#define DEFDOVE     50      // min. overlap for dovetailed reads
#define DEFMISM     0.1f    // mismatch fraction
#define OFFSET      33      // fastq quality offset (Sanger = 33)
#define DEFTHR      1       // number of threads

// fastq parts
enum fastq { HEAD, SEQ, PLUS, QUAL, FASTQ };  // lines of a fastq read
#define BEGIN       '@'     // beginning of header line
#define PLUSCHAR    '+'     // beginning of 3rd line
#define EXTRA       2       // save 2 extra strings for 2nd read:
                            //   revComp(seq) and rev(qual)

// command-line parameters
#define HELP        "-h"
#define HELP2       "--help"
#define FIRST       "-1"
#define SECOND      "-2"
#define OUTFILE     "-o"
#define UNFILE      "-f"
#define LOGFILE     "-l"
#define OVERLAP     "-m"
#define MISMATCH    "-p"
#define DOVEOPT     "-d"
#define DOVEOVER    "-e"
#define DOVEFILE    "-c"
#define MAXOPT      "-s"
#define ADAPTOPT    "-a"
#define ALNFILE     "-j"
#define GZOPT       "-g"
#define UNGZOPT     "-u"
#define DIFFOPT     "-b"
#define INTEROPT    "-t"
#define QUALITY     "-q"
#define THREADS     "-n"
#define VERBOSE     "-v"
#define VERSOPT     "--version"

// extensions for output files
#define ONEEXT      "_1.fastq"
#define TWOEXT      "_2.fastq"
#define GZEXT       ".gz"   // for gzip compression

// OMP locks
enum omp_locks { OUT, UN, LOG, DOVE, ALN, OMP_LOCKS };

// error messages
enum errCode { ERROPEN, ERRCLOSE, ERROPENW, ERRUNK, ERRMEM, ERRSEQ,
  ERRQUAL, ERRHEAD, ERRINT, ERRFLOAT, ERRPARAM, ERRPARAM2, ERROVER,
  ERRMISM, ERRFASTQ, ERROFFSET, ERRUNGET, ERRGZIP, ERRTHREAD, DEFERR };
const char* errMsg[] = { ": cannot open file for reading",
  "Cannot close file",
  ": cannot open file for writing",
  ": unknown nucleotide",
  "Cannot allocate memory",
  "Cannot load sequence",
  "Sequence/quality scores do not match",
  ": not matched in input files",
  ": cannot convert to int",
  ": cannot convert to float",
  ": unknown command-line parameter",
  ": unknown command-line parameter with no arg",
  "Overlap must be greater than 0",
  "Mismatch must be in [0,1)",
  "Input file does not follow fastq format",
  "Quality scores outside of set range",
  "Failure in ungetc() call",
  "Cannot pipe in gzip compressed file (use zcat instead)",
  "Number of threads must be >= 1",
  "Unknown error" };

typedef union file {
  FILE* f;
  gzFile gzf;
} File;
