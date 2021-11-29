/* rriScan head file */

#ifndef rriScan_HEAD_H
#define rriScan_HEAD_H

#ifndef MIN
#define MIN(x, y) (x)<(y)?(x):(y)
#endif

#ifndef MAX
#define MAX(x, y) (x)>(y)?(x):(y)
#endif

#define MATCH 2.0
#define MISMATCH -3.0
#define GU_MATCH -0.5

#define GAP_OPEN -3.0
#define GAP_CONT -1.0

#define MAX_CLIP_LEN 0

// is not used
#define SEEKMAXLOOP 20

typedef struct alignInfo
{
  char *mirSeq;
  char *tarSeq;
  char *pairStr;
  int   mirStart;
  int   mirEnd;
  int   tarStart;
  int   tarEnd;
  int   maxPairNum;
  double alignScore;
} alignInfo;

struct chimeraJunctionInfo
{
  char *lChrom;
  char *rChrom;
  int lStart;
  int lEnd;
  int rStart;
  int rEnd;
  int lLen;
  int rLen;
  char *name;
  int score;
  char lStrand;
  char rStrand;
  int alnNum;
  int gapDist;
  int isReverse;
};

typedef struct chimeraJunctionInfo junctionInfo;

typedef vector<junctionInfo *> junctionVector;

typedef struct parameterInfo
{
  int minSegmentLen;
  int verbose;
  int circRNA;
  int minReadNum;
  int peak;
  int minGapDist;
  int minPair;
  int smallGenome;
  double maxMFE;
  double minScore;
} parameterInfo;

struct readInfo
{
  char *readName;
  char *readSeq;
  char *qualityName;
  char *quality;
};

typedef struct readInfo readInfo;

typedef map<string, char*> seqMap;

typedef map<string, int> repeatMap;

extern seqMap seqHash;

void scanRRI(parameterInfo *paraInfo, FILE *gfp, FILE *faifp, char *bamFile, FILE *junfp, char *readFile, FILE *outfp);

double readJunctionToChimera(parameterInfo *paraInfo, BamReader &reader, FILE *junfp,
                             FILE *gfp, faidxMap &faiHash, FILE *outfp);

double readBamToChimera(parameterInfo *paraInfo, BamReader &globalReader, BamReader &reader,
                        FILE *gfp, faidxMap &faiHash, FILE *outfp);

int getSegmentLen(char *cigar);

void bedBlockToJunction(bed6Vector &bedBlocks, junctionInfo *jInfo);

void outputJunction(parameterInfo * paraInfo, BamReader &reader,
                    FILE * gfp, faidxMap &faiHash,
                    junctionInfo *bestJunc, FILE * outfp);

void freeJunction(junctionInfo *jInfo);

void freeJunctionVector(junctionVector &jList);

int getClipLen(char *cigar, int *leftClipLen, int *rightClipLen);

int getGenomePos(junctionInfo *jptr);

const char *getChimeraTypes(parameterInfo *paraInfo, junctionInfo *jInfo);

int padStructure(char *chimeraStruct, int lLen, int seqLen, const duplexT *dup);

int drawPairs(parameterInfo *paraInfo, const char *seq, char *structure, alignInfo * align, char *rna1, char *rna2);

double scoreAlignments(alignInfo * align);

double scorePair(char a, char b);

double scoreGap(int gapFlag);

void freeFloatMatrix(double **matrix);

void freeIntMatrix(int **matrix);

int argmax(double m[], int len);

double max(double m[], int len);

void freeAlignInfo(alignInfo *align);

int filterPairs(char *pairStr);

int encodeIntChar (char ch);

int RNApair(char bp1, char bp2);

char getPairChar(char bp1, char bp2);

double smithWatermanScorePair(parameterInfo *paraInfo, char *rna1, char *rna2, alignInfo *swAlign);

int maxPairNum(char *pairStr);

int cmpJunctionSite(const junctionInfo *x, const junctionInfo *y);

void printTagType(BamAlignment &bam, string tagName);

readInfo *getOneRead(FILE *fp, char *line, char format);

void freeRead(readInfo *read);;

int readGzipSeqToMap(char *readFile, seqMap &seqHash);

char *getGzipLine(gzFile gzf);

readInfo *getGzipOneRead(gzFile gzf, char *line, char format);

int readSeqToMap(char *readFile, seqMap &seqHash);

void freeRead(readInfo *read);

readInfo *getOneRead(FILE *fp, char *line, char format);

void freeSeqMap(seqMap & seqHash);

void outputHeader(FILE * outfp);

int cmpChromLocus(const junctionInfo *x, const junctionInfo *y);

void exchangeJunctionSite(junctionInfo *junc);

void exchangeSite(junctionInfo *junc);

#endif /* End rriScan_HEAD_H */
