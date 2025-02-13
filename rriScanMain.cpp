#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<getopt.h>
extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
#include "duplex.h"
#include "cofold.h"
}
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
#include <zlib.h>

using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "samFile.h"
#include "rriScan.h"

char version[] = "rriScan version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *outfile       = NULL;
  char *bamFile       = NULL;
  char *junFile       = NULL;
  char *readFile      = NULL;
  char *fastaFile     = NULL;
  char *faiFile       = NULL;
  FILE *outfp         = NULL;
  FILE *genomefp      = NULL;
  FILE *faifp         = NULL;
  FILE *junfp         = NULL;
  int showVersion     = 0;
  int showHelp        = 0;
  int i               = 0;
  int c               = 0;
  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVbtPSo:l:m:d:p:M:s:f:F:j:r:g:" ;

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "peak" , no_argument , NULL, 'P' },
    { "small" , no_argument , NULL, 'S' },
    { "output" , required_argument , NULL, 'o' },
    { "min-seg-len" , required_argument, NULL, 'l' },
    { "min-read-num" , required_argument, NULL, 'm' },
    { "max-dist" , required_argument, NULL, 'd' },
    { "min-pair" , required_argument , NULL, 'p' },
    { "max-mfe" , required_argument , NULL, 'M' },
    { "min-score" , required_argument , NULL, 's' },
    { "min-gap" , required_argument , NULL, 'g' },
    { "fa" , required_argument , NULL, 'f' },
    { "fai" , required_argument , NULL, 'F' },
    { "bam" , required_argument , NULL, 'b' },
    { "jun" , required_argument , NULL, 'j' },
    { "read" , required_argument , NULL, 'r' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.minSegmentLen  = 15;
  paraInfo.verbose        = 0;
  paraInfo.circRNA        = 0;
  paraInfo.minReadNum     = 1;
  paraInfo.peak           = 0;
  paraInfo.minPair        = 0;
  paraInfo.maxMFE         = -5.0;
  paraInfo.minScore       = 5.0;
  paraInfo.minGapDist     = 1;
  paraInfo.smallGenome    = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'P':
      paraInfo.peak = 1;
      break;
    case 'S':
      paraInfo.smallGenome = 1;
      break;
    case 'o':
      outfile   = optarg;
      break;
    case 'f':
      fastaFile = optarg;
      break;
    case 'F':
      faiFile   = optarg;
      break;
    case 'b':
      bamFile   = optarg;
      break;
    case 'j':
      junFile   = optarg;
      break;
    case 'r':
      readFile  = optarg;
      break;
    case 'l':
      paraInfo.minSegmentLen = atoi(optarg);
      break;
    case 'm':
      paraInfo.minReadNum    = atoi(optarg);
      break;
    case 'p':
      paraInfo.minPair       = atoi(optarg);
      break;
    case 'M':
      paraInfo.maxMFE        = atof(optarg);
      break;
    case 's':
      paraInfo.minScore      = atof(optarg);
      break;
    case 'g':
      paraInfo.minGapDist    = atoi(optarg);
      break;      
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (argc != optind) usage();

  if (fastaFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --fa <genome fasta file>\n");
    usage();
  }
  genomefp = (FILE *) fopen(fastaFile, "r");
  if (genomefp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open genome file: %s\n", fastaFile);
    usage();
  }
  if (faiFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --fai <fai file>\n");
    usage();
  }
  faifp = (FILE *) fopen(faiFile, "r");
  if (faifp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open fai file: %s\n", faiFile);
    usage();
  }

  if (junFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --jun <junction file>\n");
    usage();
  }
  junfp = (FILE *) fopen(junFile, "r");
  if (junfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open junction file: %s\n", junFile);
    usage();
  }

  if (bamFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --bam <mapped alignments, BAM format>\n");
    usage();
  }

  if (outfile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outfile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outfile);
      usage();
    }
  }

  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  fprintf(stderr, "#program start\n");
  scanRRI(&paraInfo, genomefp, faifp, bamFile, junfp, readFile, outfp);
  fprintf(stderr, "#program end\n");
  
  fclose(faifp);
  fclose(genomefp);
  fclose(outfp);

  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  rriScan [options] --fa <fasta file> --fai <fai file> --bam <mapped alignments> --jun <junctions>\n\
File format for mapped alignments is BAM\n\
[options]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : rriScan version\n\
-h/--help                      : help informations\n\
-S/--small                     : small genome\n\
--fa                           : genome FASTA file.[required]\n\
--fai                          : genome fai file, an index file for fasta file.[required]\n\
--bam                          : alignment file, BAM format.[required]\n\
--jun                          : junction file from STAR software, junction format.[required]\n\
--read                         : read file[fastq or fasta].[optional]\n\
-o/--output <string>           : output file\n\
-l/--min-seg-len <int>         : minimum length of segments in a chimera [default>=15]\n\
-m/--min-read-number <int>     : minimum read number for chimera [default>=1]\n\
-M/--max-mfe <double>          : maximum MFE in duplex[default<=-5.0]\n\
-p/--min-pair <int>            : minimum pair number in duplex [default>=0]\n\
-s/--min-score <int>           : minimum alignment score in duplex [default>=5]\n\
-g/--min-gap <int>             : minimum gaps between two segments [default>=1]\n\
");

  exit(1);
}
