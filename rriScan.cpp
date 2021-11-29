/* API for bed format */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<limits.h>
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
#include "homer_statistics.h"
#include "statistic.h"
#include "rriScan.h"

// global parameters for all functions
seqMap seqHash;

void scanRRI(parameterInfo *paraInfo, FILE *gfp, FILE *faifp, char *bamFile, FILE *junfp, char *readFile, FILE *outfp)
{
	junctionVector juncList;
	faidxMap faiHash;
	BamReader reader;
	BamReader globalReader;
	fprintf(stderr, "read fai file\n");
	readFai(faifp, faiHash);
	fprintf(stderr, "read bam file\n");
	openBamFile(bamFile, reader);
	openBamFile(bamFile, globalReader);

	if (readFile != NULL)
	{
		if (strstr(readFile, ".gz") != NULL)
		{
			fprintf(stderr, "read sequencing reads from zip file\n");
			int readNum = readGzipSeqToMap(readFile, seqHash);
			fprintf(stderr, "get %d reads\n", readNum);
		}
		else
		{
			fprintf(stderr, "read sequencing reads from file\n");
			int readNum = readSeqToMap(readFile, seqHash);
			fprintf(stderr, "get %d reads\n", readNum);
		}
	}

	outputHeader(outfp);
	fprintf(stderr, "read bam to chimera\n");
	double totalRead1 = readBamToChimera(paraInfo, globalReader, reader, gfp, faiHash, outfp);
	fprintf(stderr, "identify %.0f chimeric reads from bam file\n", totalRead1);
	fprintf(stderr, "read junction to chimera\n");
	double totalRead2 = readJunctionToChimera(paraInfo, globalReader, junfp, gfp, faiHash, outfp);
	fprintf(stderr, "identify %.0f chimeric reads from junction file\n", totalRead2);
	fprintf(stderr, "identify total: %.0f chimeric reads\n", totalRead1 + totalRead2);

	if (readFile != NULL) freeSeqMap(seqHash);
	reader.Close();
	globalReader.Close();
}

void outputHeader(FILE * outfp)
{
	fprintf(outfp, "lChrom\tlChromStart\tlChromEnd\tlName\tlScore\tlStrand\trChrom\trChromStart\trChromEnd\trName\trScore\trStrand\t");
	fprintf(outfp, "lociNum\tgapDist\t");
	if (seqHash.size() > 0) fprintf(outfp, "readSeq\t");
	fprintf(outfp, "chimericSeq\tchimericStruct\tMFE\trriType\t");
	fprintf(outfp, "lAlignSeq\tpairs\trAlignSeq\tpairNum\talignScore\tloReadNum\troReadNum\n");
}

void outputJunction(parameterInfo * paraInfo, BamReader &globalReader,
                    FILE * gfp, faidxMap &faiHash,
                    junctionInfo *bestJunc, FILE * outfp)
{
	if (paraInfo->verbose) fprintf(stderr, "output Junctions\n");

	if (paraInfo->smallGenome)
	{
		getGenomePos(bestJunc);
	}

	if (bestJunc->isReverse) exchangeSite(bestJunc);

	string lchrom(bestJunc->lChrom);
	faidx *lfai = faiHash[lchrom];
	string rchrom(bestJunc->rChrom);
	faidx *rfai = faiHash[rchrom];

	if (paraInfo->verbose) fprintf(stderr, "output best junction %s %d %d %s %d %d %d %d\n",
		                               bestJunc->lChrom, bestJunc->lStart, bestJunc->lEnd,
		                               bestJunc->rChrom, bestJunc->rStart, bestJunc->rEnd, bestJunc->lLen, bestJunc->rLen);
	char *lseq = faidxFetchSeq(gfp, lfai, bestJunc->lStart, bestJunc->lEnd, bestJunc->lStrand);
	char *rseq = faidxFetchSeq(gfp, rfai, bestJunc->rStart, bestJunc->rEnd, bestJunc->rStrand);
	int duplex = 1;
	//int lds = 0, lde = 0, rds = 0, rde = 0;
	double mfe = 0;
	int lLen = strlen(lseq);
	int rLen = strlen(rseq);
	int seqLen = lLen + rLen + 1;
	char *chimeraSeq    = (char *)safeMalloc(seqLen + 1);
	char *chimeraStruct = (char *)safeMalloc(seqLen + 1);
	strcpy(chimeraSeq, lseq);
	strcat(chimeraSeq, "&");
	strcat(chimeraSeq, rseq);
	chimeraSeq[seqLen] = '\0';
	if (duplex == 0)
	{
		mfe = fold(chimeraSeq, chimeraStruct);
	}
	else
	{
		duplexT dup = duplexfold(lseq, rseq);
		mfe = dup.energy;
		padStructure(chimeraStruct, lLen, seqLen, &dup);
		safeFree(dup.structure);
	}
	char *chimeraType  = (char *)getChimeraTypes(paraInfo, bestJunc);
	alignInfo *swAlign = (alignInfo *)safeMalloc(sizeof(alignInfo));
	int pairScore = drawPairs(paraInfo, chimeraSeq, chimeraStruct, swAlign, lseq, rseq);
	if (swAlign->alignScore >= paraInfo->minScore
	        && swAlign->maxPairNum >= paraInfo->minPair)
	{
		int lReadNum = 0;
		int rReadNum = 0;
		//int lReadNum = getRegionCounts(globalReader, bestJunc->lChrom, bestJunc->lStart, bestJunc->lEnd, bestJunc->lStrand);
		//int rReadNum = getRegionCounts(globalReader, bestJunc->rChrom, bestJunc->rStart, bestJunc->rEnd, bestJunc->rStrand);
		fprintf(outfp, "%s\t%d\t%d\t%s_lr\t%d\t%c\t%s\t%d\t%d\t%s_rr\t%d\t%c\t%d\t%d\t",
		        bestJunc->lChrom, bestJunc->lStart, bestJunc->lEnd, bestJunc->name,
		        bestJunc->lLen, bestJunc->lStrand,
		        bestJunc->rChrom, bestJunc->rStart, bestJunc->rEnd, bestJunc->name,
		        bestJunc->rLen, bestJunc->rStrand, bestJunc->alnNum, bestJunc->gapDist);
		if (seqHash.size() > 0)
		{
			string readName(bestJunc->name);
			if (seqHash.find(readName) != seqHash.end())
			{
				char *seq = seqHash[readName];
				fprintf(outfp, "%s\t", seq);
			}
		}
		fprintf(outfp, "%s&%s\t%s\t%.5f\t%s\t", lseq, rseq, chimeraStruct, mfe, chimeraType);
		fprintf(outfp, "RNA1  5'>%s>3'\tPair     %s\tRNA2  3'<%s<5'\t%d\t%.2f\t", swAlign->mirSeq, swAlign->pairStr, swAlign->tarSeq, swAlign->maxPairNum, swAlign->alignScore);
		fprintf(outfp, "%d\t%d\n", lReadNum, rReadNum);
		fflush(outfp);
	}
	free_arrays;
	//if (duplex) safeFree(dup.structure);
	safeFree(chimeraSeq);
	safeFree(chimeraStruct);
	safeFree(lseq);
	safeFree(rseq);
	freeAlignInfo(swAlign);
}

// for small genome
int getGenomePos(junctionInfo *jptr)
{
	int fieldNum  = 0;
	char **fields = NULL;
	char delims[] = ":";

	// for left chrom
	int start   = 0;
	int end     = 0;
	char strand = '+';
	int tag = 0;
	fields = splitString(jptr->lChrom, delims, &fieldNum);
	if (fieldNum >= 4) {
		tag = 1;
		safeFree(jptr->lChrom);
		jptr->lChrom = strClone(fields[1]);
		char *posStr = fields[2];
		char *cp = strchr(posStr, '-');
		if (cp != NULL) {
			char cs = *cp;
			*cp   = '\0';
			start  = atoi(posStr) - 1;
			end    = atoi(cp + 1);
			*cp   = cs;
		}
		strand = fields[3][0];
		int lStart = 0;
		int lEnd = 0;
		if (strand == '+') {
			lStart = start + jptr->lStart;
			lEnd = start + jptr->lEnd;
		}
		else {
			lStart = end - jptr->lEnd;
			lEnd = end - jptr->lStart;
		}
		if (lStart < 0) lStart = 0;
		if (lEnd < 0) lEnd = 0;
		jptr->lStart = lStart;
		jptr->lEnd = lEnd;
		if (strand == '-')
		{
			if (jptr->lStrand == '+')
			{
				jptr->lStrand = '-';
			}
			else
			{
				jptr->lStrand = '+';
			}
		}
	}
	freeWords(fields, fieldNum);

	start   = 0;
	end     = 0;
	strand  = '+';
	tag     = 0;

	// for right chrom
	fields = splitString(jptr->rChrom, delims, &fieldNum);
	if (fieldNum >= 4) {
		tag = 1;
		safeFree(jptr->rChrom);
		jptr->rChrom = strClone(fields[1]);
		char *posStr = fields[2];
		char *cp = strchr(posStr, '-');
		if (cp != NULL) {
			char cs = *cp;
			*cp   = '\0';
			start  = atoi(posStr) - 1;
			end    = atoi(cp + 1);
			*cp   = cs;
		}
		strand = fields[3][0];
		int rStart = 0;
		int rEnd = 0;
		if (strand == '+') {
			rStart = start + jptr->rStart;
			rEnd = start + jptr->rEnd;
		}
		else {
			rStart = end - jptr->rEnd;
			rEnd = end - jptr->rStart;
		}
		if (rStart < 0) rStart = 0;
		if (rEnd < 0) rEnd = 0;
		jptr->rStart = rStart;
		jptr->rEnd = rEnd;
		if (strand == '-')
		{
			if (jptr->rStrand == '+')
			{
				jptr->rStrand = '-';
			}
			else
			{
				jptr->rStrand = '+';
			}
		}
	}
	freeWords(fields, fieldNum);

	return tag;
}

int readGzipSeqToMap(char *readFile, seqMap &seqHash)
{
	int readNum = 0;
	gzFile gzf;
	gzf = gzopen(readFile, "rb");
	if (gzf == NULL) {
		fprintf(stderr, "Error: gzopen error\n");
		exit(1);
	}

	char *line = NULL;
	while (line = getGzipLine(gzf))
	{
		if (gzeof(gzf) || line == NULL)
		{
			safeFree(line);
			break;
		}
		if (line[0] != '@' && line[0] != '>')
		{
			fprintf(stderr, "error read format: %c\n", line[0]);
			safeFree(line);
			break;
		}
		readInfo *read = getGzipOneRead(gzf, line, line[0]);
		string name(read->readName);
		seqHash[name] = strClone(read->readSeq);
		freeRead(read);
		readNum++;
	}
	gzclose(gzf);
	return readNum;
}

char *getGzipLine(gzFile gzf)
/* get whole line from gzip
*/
{
	const int lineLen = 512;
	char s[lineLen];
	char *line = NULL;
	char *cp = NULL;
	int done = FALSE;
	line = NULL;

	do {
		if (gzgets(gzf, s, lineLen) == Z_NULL)
			break; /* EOF */
		/* for unix OS */
		cp = strchr(s, '\n');
		if (cp != NULL) {
			*cp   = '\0';
			done  = TRUE;
		}
		/* for window OS */
		cp = strchr(s, '\r');
		if (cp != NULL) {
			*cp  = '\0';
			done = TRUE;
		}
		if (line == NULL)
			line = (char *)safeMalloc(strlen(s) + 1); /* don't use malloc, for we will using strcat function, so we must initilized the line with '0' */
		else
			line = (char *)safeRealloc(line, strlen(s) + strlen(line) + 1);
		strcat(line, s);
	} while (!done);

	return line;
}

readInfo *getGzipOneRead(gzFile gzf, char *line, char format)
{
	char *seq = NULL;
	char *qualityName = NULL;
	char *quality = NULL;
	seq = getGzipLine(gzf);
	if (format == '@')
	{
		qualityName = getGzipLine(gzf);
		quality = getGzipLine(gzf);
	}
	readInfo *read = (readInfo *)safeMalloc(sizeof(readInfo));
	if (format == '>')
	{
		read->readName = strClone(line + 1);
	}
	else
	{
		read->readName = strClone(line + 1);
	}
	read->readSeq     = seq;
	read->qualityName = qualityName;
	read->quality     = quality;
	return read;
}


int readSeqToMap(char *readFile, seqMap &seqHash)
{
	int readNum = 0;
	FILE *readfp = NULL;
	readfp = (FILE *) fopen(readFile, "r");
	if (readfp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open read file: %s\n", readFile);
		exit(1);
	}

	char *line = NULL;
	while (line = getLine(readfp))
	{
		if (feof(readfp) || line == NULL)
		{
			safeFree(line);
			break;
		}
		if (line[0] != '@' && line[0] != '>')
		{
			fprintf(stderr, "error read format: %c\n", line[0]);
			safeFree(line);
			break;
		}
		readInfo *read = getOneRead(readfp, line, line[0]);
		string name(read->readName);
		seqHash[name] = strClone(read->readSeq);
		freeRead(read);
		readNum++;
	}
	fclose(readfp);

	return readNum;
}

void freeRead(readInfo *read)
{
	safeFree(read->readName);
	safeFree(read->readSeq);
	if (read->qualityName != NULL && read->quality != NULL)
	{
		safeFree(read->qualityName);
		safeFree(read->quality);
	}
	safeFree(read);
}

readInfo *getOneRead(FILE *fp, char *line, char format)
{
	char *seq = NULL;
	char *qualityName = NULL;
	char *quality = NULL;
	seq = getLine(fp);
	if (format == '@')
	{
		qualityName = getLine(fp);
		quality = getLine(fp);
	}
	readInfo *read = (readInfo *)safeMalloc(sizeof(readInfo));
	if (format == '>')
	{
		read->readName = strClone(line + 1);
	}
	else
	{
		read->readName = strClone(line + 1);
	}
	read->readSeq     = seq;
	read->qualityName = qualityName;
	read->quality     = quality;
	return read;
}

double readJunctionToChimera(parameterInfo *paraInfo, BamReader &globalReader, FILE *junfp,
                             FILE *gfp, faidxMap &faiHash, FILE *outfp)
// must at same chromosome
{
	int i = 0;
	int fieldNum = 0;
	char **fields = NULL;
	int minTermMis  = MAX_CLIP_LEN;
	bool skipDeletion = 0, skipSplice = 0, collapser = 0;
	double totalNum  = 0;
	char *line = NULL;

	char readSeqName[1204];
	char repSeqName[1204];
	repSeqName[0] = 'N';
	repSeqName[1] = '\0';
	junctionVector  jVector;

	while (line = getLine(junfp))
	{
		int leftClipLen = 0, rightClipLen = 0;
		if (feof(junfp) || line == NULL)
		{
			safeFree(line);
			break;
		}
		fields = splitWhitespace(line, &fieldNum);
		if (fieldNum < 14)
		{
			freeWords(fields, fieldNum);
			fprintf(stderr, "the format of annotation file must be junction generated from STAR\n");
			exit(1);
		}
		if (fieldNum > 18) // kept best junctions
		{
			int thisScore = atoi(fields[17]);
			int bestScore = atoi(fields[18]);
			if (thisScore != bestScore)
			{
				safeFree(line);
				freeWords(fields, fieldNum);
				continue;
			}
		}

		junctionInfo *jInfo = (junctionInfo *)safeMalloc(sizeof(junctionInfo));
		jInfo->lChrom  = strClone(fields[0]);
		jInfo->lStrand = fields[2][0];
		jInfo->rChrom  = strClone(fields[3]);
		jInfo->rStrand = fields[5][0];
		jInfo->name    = strClone(fields[9]);
		jInfo->lStart  = atoi(fields[10]) - 1;
		jInfo->rStart  = atoi(fields[12]) - 1;

		char *lCigar   = fields[11];
		char *rCigar   = fields[13];

		jInfo->alnNum = 1;
		if (fieldNum == 15)
			jInfo->alnNum = atoi(fields[14]);

		int lLen = getSegmentLen(lCigar);
		int rLen = getSegmentLen(rCigar);

		jInfo->lEnd    = jInfo->lStart + lLen;
		jInfo->rEnd    = jInfo->rStart + rLen;

		jInfo->lLen    = lLen;
		jInfo->rLen    = rLen;

		if (strchr(lCigar, 'N') != NULL || strchr(rCigar, 'N') != NULL)
		{
			safeFree(line);
			freeWords(fields, fieldNum);
			freeJunction(jInfo);
			continue;
		}

		getClipLen(lCigar, &leftClipLen, &rightClipLen);
		if (leftClipLen > minTermMis && rightClipLen > minTermMis)
		{
			safeFree(line);
			freeWords(fields, fieldNum);
			freeJunction(jInfo);
			continue;
		}

		leftClipLen = 0, rightClipLen = 0;
		getClipLen(rCigar, &leftClipLen, &rightClipLen);
		if (leftClipLen > minTermMis && rightClipLen > minTermMis)
		{
			safeFree(line);
			freeWords(fields, fieldNum);
			freeJunction(jInfo);
			continue;
		}

		if (strcmp(jInfo->lChrom, jInfo->rChrom) == 0)
			jInfo->gapDist  = overlapLength(jInfo->lStart, jInfo->lEnd, jInfo->rStart, jInfo->rEnd);
		else
			jInfo->gapDist  = -10000000;

		if (jInfo->lLen >= paraInfo->minSegmentLen
		        && jInfo->rLen >= paraInfo->minSegmentLen
		        && abs(jInfo->gapDist) >= paraInfo->minGapDist)
		{
			exchangeJunctionSite(jInfo);
			if (strcmp(repSeqName, jInfo->name) == 0) // almost full length
			{
				jVector.push_back(jInfo);
				strcpy(repSeqName, jInfo->name);
			}
			else
			{
				if (jVector.size() > 0)
				{
					sort(jVector.begin(), jVector.end(), cmpChromLocus);
					jVector[0]->alnNum = jVector.size();
					outputJunction(paraInfo, globalReader, gfp, faiHash, jVector[0], outfp); //  same read, only kept first
					freeJunctionVector(jVector);
					totalNum += 1;
				}
				jVector.push_back(jInfo);
				strcpy(repSeqName, jInfo->name);
			}
		}
		else
		{
			freeJunction(jInfo);
		}
		// free line and fields
		safeFree(line);
		freeWords(fields, fieldNum);
	} // while end
	if (jVector.size() > 0)
	{
		sort(jVector.begin(), jVector.end(), cmpChromLocus);
		jVector[0]->alnNum = jVector.size();
		outputJunction(paraInfo, globalReader, gfp, faiHash, jVector[0], outfp); //  same read, only kept  first
		freeJunctionVector(jVector);
		totalNum += 1;
	}
	return totalNum;
}

int cmpChromLocus(const junctionInfo *x, const junctionInfo *y)
{
	if (strcmp(x->lChrom, y->lChrom) != 0)
		return strcmp(x->lChrom, y->lChrom) < 0;

	if (strcmp(x->rChrom, y->rChrom) != 0)
		return strcmp(x->rChrom, y->rChrom) < 0;

	if (x->lStart != y->lStart)
		return x->lStart < y->lStart;

	return x->rStart < y->rStart;
}

void exchangeJunctionSite(junctionInfo *junc)
{

	if (strcmp(junc->lChrom, junc->rChrom) > 0)
	{
		junc->isReverse = 1;
		exchangeSite(junc);
	}
	if (strcmp(junc->lChrom, junc->rChrom) == 0
	        && junc->lStart > junc->rStart)
	{
		junc->isReverse = 1;
		exchangeSite(junc);
	}
}

void exchangeSite(junctionInfo *junc)
{
	char *tmpChrom;
	int tmpStart;
	int tmpEnd;
	char strand;
	tmpChrom = junc->lChrom;
	tmpStart = junc->lStart;
	tmpEnd = junc->lEnd;
	strand = junc->lStrand;
	junc->lChrom = junc->rChrom;
	junc->lStart = junc->rStart;
	junc->lEnd = junc->rEnd;
	junc->lStrand = junc->rStrand;
	junc->rChrom = tmpChrom;
	junc->rStart = tmpStart;
	junc->rEnd = tmpEnd;
	junc->rStrand = strand;
}

int getSegmentLen(char *cigar)
{
	int i = 0;
	int start = 0;
	int cigarLen = strlen(cigar);
	int segLen = 0;
	for (i = 0; i < cigarLen; i++)
	{
		char c = cigar[i];
		if (!isdigit(c)) {
			cigar[i] = '\0';
			int len = atoi(cigar + start);
			switch ( c ) {
			//increase end position on CIGAR chars [DMXN=]
			case ('M') :
				segLen  += len;
				break;
			case ('=') :
				segLen  += len;
				break;
			case ('D') :
				segLen  += len;
				break;
			}
			cigar[i] = c;
			start = i + 1;
		}// if end
	}
	return segLen;
}

double readBamToChimera(parameterInfo *paraInfo, BamReader &globalReader, BamReader &reader,
                        FILE *gfp, faidxMap &faiHash, FILE *outfp)
// must at same chromosome
{
	int minTermMis  = MAX_CLIP_LEN;
	bool skipDeletion = 1, skipSplice = 0, collapser = 0; // skip the deletion
	double totalNum  = 0;
	repeatMap repeatHash;
	// get header & reference information
	reader.Rewind();
	string header  = reader.GetHeaderText();
	RefVector refs = reader.GetReferenceData();
	// rip through the BAM file and convert each mapped entry to BED
	BamAlignment bam;
	while (reader.GetNextAlignment(bam))
	{
		if (bam.IsMapped() == true)
		{
			string readName    = bam.Name;
			if (repeatHash.find(readName) != repeatHash.end())
			{
				continue;
			}
			else
			{
				repeatHash[readName] = 1;
			}
			int leftClipLen = 0, rightClipLen = 0;
			string JmbTag = "jM";
			vector<int8_t> jmbVals;
			string cigar = BuildCigarString(bam.CigarData);
			char *cigarPtr = const_cast<char *>(cigar.c_str());
			// skip the reads without gaps
			if (strchr(cigarPtr, 'N') == NULL) continue;

			if (bam.GetTag(JmbTag, jmbVals))
			{
				if (jmbVals[0] >= 20) continue; // annotated splicing sites
			}
			else
			{
				continue;
			}

			getClipLen(cigarPtr, &leftClipLen, &rightClipLen);

			if (leftClipLen > minTermMis || rightClipLen > minTermMis) continue;
			string chrom = refs.at(bam.RefID).RefName;

			bed6Vector bedBlocks;
			bamToBlocks(bam, chrom, bedBlocks, skipDeletion, skipSplice);
			int blockCount = bedBlocks.size();

			if (blockCount != 2)
			{
				freeBed6Vector(bedBlocks);
				continue;
			}

			junctionInfo *jInfo = (junctionInfo *)safeMalloc(sizeof(junctionInfo));
			bedBlockToJunction(bedBlocks, jInfo);

			jInfo->alnNum = 1;
			string nhStr = "NH";
			uint8_t nhVal = 1;
			if (!bam.GetTag(nhStr, nhVal)) {
				fprintf(stderr, "The requested tag NH is not exist in bam file\n");
			}
			else {
				jInfo->alnNum = nhVal;
			}
			jInfo->isReverse = 0;
			if (jInfo->lLen >= paraInfo->minSegmentLen
			        && jInfo->rLen >= paraInfo->minSegmentLen
			        && abs(jInfo->gapDist) >= paraInfo->minGapDist)
				outputJunction(paraInfo, globalReader, gfp, faiHash, jInfo, outfp);

			freeBed6Vector(bedBlocks);
			freeJunction(jInfo);

			totalNum++;
		}
	}
	return totalNum;
}

void printTagType(BamAlignment &bam, string tagName)
{
	char type;
	if ( bam.GetTagType(tagName, type) )
	{
		switch ( type ) {
		case (Constants::BAM_TAG_TYPE_ASCII):
			fprintf(stderr, "BAM_TAG_TYPE_ASCII\n");
			break;
		case (Constants::BAM_TAG_TYPE_INT8):
			fprintf(stderr, "BAM_TAG_TYPE_INT8\n");
			break;
		case (Constants::BAM_TAG_TYPE_UINT8):
			fprintf(stderr, "BAM_TAG_TYPE_UINT8\n");
			break;
		case (Constants::BAM_TAG_TYPE_INT16):
			fprintf(stderr, "BAM_TAG_TYPE_INT16\n");
			break;
		case (Constants::BAM_TAG_TYPE_UINT16):
			fprintf(stderr, "BAM_TAG_TYPE_UINT16\n");
			break;
		case (Constants::BAM_TAG_TYPE_INT32):
			fprintf(stderr, "BAM_TAG_TYPE_UINT16\n");
			break;
		case (Constants::BAM_TAG_TYPE_UINT32):
			fprintf(stderr, "BAM_TAG_TYPE_UINT32\n");
			break;
		case (Constants::BAM_TAG_TYPE_FLOAT):
			fprintf(stderr, "BAM_TAG_TYPE_FLOAT\n");
			break;
		case (Constants::BAM_TAG_TYPE_STRING):
			fprintf(stderr, "BAM_TAG_TYPE_STRING\n");
			break;
		case (Constants::BAM_TAG_TYPE_ARRAY):
			fprintf(stderr, "BAM_TAG_TYPE_ARRAY\n");
			break;
		case (Constants::BAM_TAG_TYPE_HEX):
			fprintf(stderr, "BAM_TAG_TYPE_STRING\n");
			break;
		} // swith option
	}
}

void bedBlockToJunction(bed6Vector &bedBlocks, junctionInfo *jInfo)
{
	CBed6* lsam = bedBlocks[0];
	CBed6* rsam = bedBlocks[1];
	if (bedBlocks[0]->strand == '-')
	{
		lsam = bedBlocks[1];
		rsam = bedBlocks[0];
	}
	jInfo->lChrom  = strClone(lsam->chrom);
	jInfo->rChrom  = strClone(rsam->chrom);
	jInfo->lStart  = lsam->chromStart;
	jInfo->rStart  = rsam->chromStart;
	jInfo->lEnd    = lsam->chromEnd;
	jInfo->rEnd    = rsam->chromEnd;
	jInfo->lStrand = lsam->strand;
	jInfo->rStrand = rsam->strand;
	jInfo->name    = strClone(lsam->name);
	jInfo->lLen    = jInfo->lEnd - jInfo->lStart;
	jInfo->rLen    = jInfo->rEnd - jInfo->rStart;
	if (strcmp(jInfo->lChrom, jInfo->rChrom) == 0)
		jInfo->gapDist  = overlapLength(jInfo->lStart, jInfo->lEnd, jInfo->rStart, jInfo->rEnd);
	else
		jInfo->gapDist  = -10000000;
}

int cmpJunctionSite(const junctionInfo *x, const junctionInfo *y)
{
	return x->score > y->score;
}

int padStructure(char *chimeraStruct, int lLen, int seqLen, const duplexT *dup)
{
	int i = 0;
	int j = 0;
	char *dupStr = dup->structure;
	int dupLen = strlen(dupStr);
	int ss = strchr(dupStr, '&') - dupStr;
	int lds = dup->i + 1 - ss; //left start
	int lde = dup->i; // left end
	int rds = dup->j; // right start
	int rde = dup->j + dupLen - ss - 2; // right end

	for (i = 0; i < seqLen; i++)
		chimeraStruct[i] = '.';

	for (i = lds - 1, j = 0; i < seqLen && j < ss; i++, j++)
		chimeraStruct[i] = dupStr[j];
	chimeraStruct[lLen] = '&';

	for (i = lLen + rds, j = ss + 1; i < seqLen && j < dupLen; i++, j++)
		chimeraStruct[i] = dupStr[j];

	chimeraStruct[seqLen] = '\0';

	return dupLen;
}

const char *getChimeraTypes(parameterInfo * paraInfo, junctionInfo * jInfo)
{
	int maxGeneLen = 50000;
	int sameChrom  = 0;

	const char *chimeraType[] = {"sameSeq", "antiRRI", "intraRRI", "intraRRI", "interRRI"};
	// sameSeq
	if (strcmp(jInfo->lChrom, jInfo->rChrom) == 0) sameChrom = 1;
	if (sameChrom == 1)
	{
		int distance = overlapLength(jInfo->lStart, jInfo->lEnd, jInfo->rStart, jInfo->rEnd);
		int absDist = abs(distance);

		if (distance >= 0 && jInfo->lStrand == jInfo->rStrand)
		{
			return chimeraType[0];
		}
		else if (distance >= 0 && jInfo->lStrand != jInfo->rStrand)
		{
			return chimeraType[1];
		}
		else if  (distance < 0)
		{
			if (jInfo->lStrand == jInfo->rStrand && absDist <= maxGeneLen)
			{
				if (jInfo->lStrand == '+')
				{
					if (jInfo->lStart < jInfo->rStart) return chimeraType[2]; // cis splicing in same gene
					if (jInfo->lStart > jInfo->rStart) return chimeraType[3]; // back splicing, such as circRNA
				}
				else
				{
					if (jInfo->lStart < jInfo->rStart) return chimeraType[3]; // back splicing, such as circRNA
					if (jInfo->lStart > jInfo->rStart) return chimeraType[2]; // cis splicing in same gene
				}
			}
			else
			{
				return chimeraType[4];
			}
		}
	}
	else
	{
		return chimeraType[4];
	}
	return chimeraType[4];
}

void freeJunction(junctionInfo * jInfo)
{
	safeFree(jInfo->lChrom);
	safeFree(jInfo->rChrom);
	safeFree(jInfo->name);
	safeFree(jInfo);
}

void freeJunctionVector(junctionVector & jList)
{
	for (junctionVector::iterator vecItr = jList.begin(); vecItr != jList.end(); vecItr++)
	{
		junctionInfo *junc = *vecItr;
		freeJunction(junc);
	}
	jList.clear();
}

void freeSeqMap(seqMap & seqHash)
{
	for (seqMap::iterator mapItr = seqHash.begin(); mapItr != seqHash.end(); mapItr++) {
		char* seq = mapItr->second;
		safeFree(seq);
	}
	seqHash.clear();
}

int getClipLen(char *cigar, int *leftClipLen, int *rightClipLen)
// get left or right soft length
{
	int i = 0;
	int start = 0;
	int tag = 0;
	for (i = 0; i < strlen(cigar); i++)
	{
		char c = cigar[i];
		if (!isdigit(c))
		{
			tag += 1;
			cigar[i] = '\0';
			int len = atoi(cigar + start);
			switch ( c )
			{
			case 'S' :
				if (tag == 1)
				{
					*leftClipLen = len;
				}
				else
				{
					*rightClipLen = len;
				}
				break;
			}
			cigar[i] = c;
			start = i + 1;
		}
	}
	if (*leftClipLen > 0 || *rightClipLen > 0) return 1;
	return 0;
}

int drawPairs(parameterInfo *paraInfo, const char *seq, char *structure, alignInfo * align, char *rna1, char *rna2)
{
	int i, j, k, flag, m;
	short int *pairTable = make_pair_table(structure);
	int seqLen = strlen(seq);
	char *qSeq = NULL;
	char *mSeq = NULL;
	char *tSeq = NULL;
	int allocSpace = seqLen * 2;
	qSeq = (char *)safeMalloc(sizeof(char) * (allocSpace));
	mSeq = (char *)safeMalloc(sizeof(char) * (allocSpace));
	tSeq = (char *)safeMalloc(sizeof(char) * (allocSpace));
	int preLeftPos  = 0;
	int preRightPos = 0;
	int linkLen = 1;
	int pairNum = 0;
	i = 0;
	j = 0;
	k = 0;
	m = 0;
	flag = 0; // start flag
	int extendLen = strlen(rna1);
	int matchLen = extendLen + linkLen;
	for (i = 0; i < extendLen && i < strlen(structure); i++)
	{
		if (k >= allocSpace)
		{
			qSeq = (char *)safeRealloc(qSeq, k + allocSpace);
			mSeq = (char *)safeRealloc(mSeq, k + allocSpace);
			tSeq = (char *)safeRealloc(tSeq, k + allocSpace);
			allocSpace = k + allocSpace;
		}
		j = i + 1;
		if (structure[i] == '(' && pairTable[j] > matchLen)
		{
			if (flag == 1) // for unpaired
			{
				int lLen = i - preLeftPos - 1;
				int rLen = preRightPos - pairTable[j] - 1;
				int maxLen = MAX(lLen, rLen);
				for (m = 0; m < maxLen; m++)
				{
					if (m < lLen)
					{
						qSeq[k] = seq[preLeftPos + m + 1];
					}
					else
					{
						qSeq[k] = '-';
					}
					if (m < rLen)
					{
						tSeq[k] = seq[pairTable[j] + rLen - 1 - m];
					}
					else
					{
						tSeq[k] = '-';
					}
					mSeq[k] = '.';
					k++;
				}
			}
			qSeq[k] = seq[i];
			tSeq[k] = seq[pairTable[j] - 1];
			mSeq[k] = getPairChar(qSeq[k], tSeq[k]);
			pairNum += RNApair(qSeq[k], tSeq[k]);
			k++;
		}
		if (structure[i] == '(' && pairTable[j] < matchLen)
		{
			if (flag == 1)
			{
				int lLen = i - preLeftPos - 1;
				for (m = 0; m < lLen; m++)
				{
					qSeq[k] = seq[preLeftPos + m + 1];
					mSeq[k] = '.';
					tSeq[k] = '-';
					k++;
				}
			}
			int inPairNum = 0;
			for (m = i; m < pairTable[j]; m++)
			{
				qSeq[k] = seq[m];
				mSeq[k] = '.';
				tSeq[k] = '-';
				inPairNum++;
				k++;
			}
			i = i + inPairNum;
			preLeftPos  = i - 1;
		}
		if (structure[i] == '(')
		{
			preLeftPos  = i;
			preRightPos = pairTable[j];
			flag = 1;
		}
	}
	safeFree(pairTable);
	qSeq[k] = '\0';
	mSeq[k] = '\0';
	tSeq[k] = '\0';

	// follows: padding the rna pairs in two ends
	int rna1Len = strlen(rna1);
	int rna2Len = strlen(rna2);

	int rna1LeftLen  = 0;
	int rna2LeftLen  = 0;
	int rna1RightLen = 0;
	int rna2RightLen = 0;

	for (i = 0; i < rna1Len; i++)
	{
		if (structure[i] == '.')
		{
			rna1LeftLen++;
		}
		else
		{
			break;
		}
	}
	for (i = rna1Len - 1; i >= 0; i--)
	{
		if (structure[i] == '.')
		{
			rna1RightLen++;
		}
		else
		{
			break;
		}
	}

	for (i = rna1Len + 1; i < seqLen; i++)
	{
		if (structure[i] == '.')
		{
			rna2RightLen++;
		}
		else
		{
			break;
		}
	}
	for (i = seqLen - 1; i >= rna1Len + 1; i--)
	{
		if (structure[i] == '.')
		{
			rna2LeftLen++;
		}
		else
		{
			break;
		}
	}

	int leftLen  = MAX(rna1LeftLen, rna2LeftLen);
	int rightLen = MAX(rna1RightLen, rna2RightLen);

	char *qlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));
	char *mlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));
	char *tlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));

	char *qrSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));
	char *mrSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));
	char *trSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));

	for (i = 0; i < leftLen; i++)
	{
		qlSeq[i] = '-';
		mlSeq[i] = '.';
		tlSeq[i] = '-';
	}
	for (i = 0; i < rightLen; i++)
	{
		qrSeq[i] = '-';
		mrSeq[i] = '.';
		trSeq[i] = '-';
	}
	for (i = rna1LeftLen, j = leftLen; i > 0 && j > 0; i--, j--)
	{
		qlSeq[j - 1] = rna1[i - 1];
	}
	for (i = rna2Len - rna2LeftLen, j = leftLen; i < rna2Len && j > 0; i++, j--)
	{
		tlSeq[j - 1] = rna2[i];
	}
	for (i = 0, j = rna1Len - rna1RightLen; i < rna1RightLen && j < rna1Len; i++, j++)
	{
		qrSeq[i] = rna1[j];
	}
	for (i = 0, j = rna2RightLen; i < rna2RightLen && j > 0; i++, j--)
	{
		trSeq[i] = rna2[j - 1];
	}

	qlSeq[leftLen] = '\0';
	mlSeq[leftLen] = '\0';
	tlSeq[leftLen] = '\0';
	qrSeq[rightLen] = '\0';
	mrSeq[rightLen] = '\0';
	trSeq[rightLen] = '\0';

	int pairLen = strlen(qSeq);
	int extLen = pairLen + leftLen + rightLen;
	char *nqSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));
	char *nmSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));
	char *ntSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));

	strcpy(nqSeq, qlSeq);
	strcpy(nmSeq, mlSeq);
	strcpy(ntSeq, tlSeq);
	strcat(nqSeq, qSeq);
	strcat(nmSeq, mSeq);
	strcat(ntSeq, tSeq);
	strcat(nqSeq, qrSeq);
	strcat(nmSeq, mrSeq);
	strcat(ntSeq, trSeq);

	nqSeq[extLen] = '\0';
	nmSeq[extLen] = '\0';
	ntSeq[extLen] = '\0';

	align->mirSeq  = nqSeq;
	align->tarSeq  = ntSeq;
	align->pairStr = nmSeq;

	align->alignScore = scoreAlignments(align);
	align->maxPairNum = maxPairNum(align->pairStr);

	safeFree(qSeq);
	safeFree(tSeq);
	safeFree(mSeq);
	safeFree(qlSeq);
	safeFree(tlSeq);
	safeFree(mlSeq);
	safeFree(qrSeq);
	safeFree(trSeq);
	safeFree(mrSeq);

	return pairNum;
}

double scoreAlignments(alignInfo * align)
{
	int i = 0;
	int gapFlag = 0;
	double score = 0;
	double bestScore = 0;
	int aliLen = strlen(align->mirSeq);
	for (i = 0; i < aliLen; i++)
	{
		char a = align->mirSeq[i];
		char b = align->tarSeq[i];

		if (a == '-' || b == '-')
		{
			score += scoreGap(gapFlag);
			gapFlag = 1;
		}
		else
		{
			score += scorePair(a, b);
			gapFlag = 0;
		}
		align->pairStr[i] = getPairChar(a, b); // correct the pairs missed by RNAduplex
		if (score < 0) score = 0;
		if (score >= bestScore)
		{
			bestScore = score;
		}
	}
	return bestScore;
}

void freeAlignInfo(alignInfo * align)
{
	safeFree(align->mirSeq);
	safeFree(align->tarSeq);
	safeFree(align->pairStr);
	safeFree(align);
}

double scoreGap(int gapFlag)
{
	if (gapFlag == 0)
	{
		return (GAP_OPEN + GAP_CONT);
	}
	else return GAP_CONT;
}

double scorePair(char a, char b)
/* score match and mismatch */
{
	int matchNum = RNApair(a, b);
	if (matchNum == 2)
	{
		return MATCH;
	}
	else if (matchNum == 1)
	{
		return GU_MATCH;
	}
	else
	{
		return MISMATCH;
	}
}

void freeFloatMatrix(double **matrix)
/* free float matrix */
{
	free(matrix[0]);
	free(matrix);
}

void freeIntMatrix(int **matrix)
/*free int matrix */
{
	free(matrix[0]);
	free(matrix);
}

double max(double m[], int len)
/* return maxinum value */
{
	int i = 0;
	double maxScore = m[0];

	for (i = 1; i < len; i++)
	{
		if (m[i] > maxScore)
		{
			maxScore = m[i];
		}
	}
	return maxScore;
}

int encodeIntChar (char ch)
{
	ch = toupper(ch);
	if (ch == 'A')
		return 0;
	else if (ch == 'C')
		return 1;
	else if (ch == 'G')
		return 2;
	else if (ch == 'T' || ch == 'U')
		return 3;
	else
		return 4;
}

char getPairChar(char bp1, char bp2)
{
	char s = '.';
	if (RNApair(bp1, bp2) == 2)
	{
		s = '|';
	}
	else if (RNApair(bp1, bp2) == 1)
	{
		s = ':';
	}
	else
	{
		s = '.';
	}
	return s;
}

int RNApair(char bp1, char bp2)
{
	int pairMatrix[5][5] =
	{
		/* A C G T N*/
		{0, 0, 0, 2, 0},
		{0, 0, 2, 0, 0},
		{0, 2, 0, 1, 0},
		{2, 0, 1, 0, 0},
		{0, 0, 0, 0, 0}
	};
	return pairMatrix[(int)(encodeIntChar(bp1))][(int)(encodeIntChar(bp2))];
}

double smithWatermanScorePair(parameterInfo * paraInfo, char *rna1, char *rna2, alignInfo * swAlign)
{
	double **scoreMatrix;
	double tmp[4];
	int i, j, k;
	int rna1Len = strlen(rna1);
	int rna2Len = strlen(rna2);
	int M = rna1Len + 1;
	int N = rna2Len + 1;
	char *newRna1 = strClone(rna1);
	reverseBytes(newRna1, rna1Len);
	char *seq1 = newRna1;
	char *seq2 = rna2;
	int allocSpace = 0;
	int extSpace = 10;
	double bestScore = 0;
	int bestI = 0;
	int bestJ = 0;
	int traceI = 0;
	int traceJ = 0;
	char *qSeq = NULL;
	char *mSeq = NULL;
	char *tSeq = NULL;
	// allocate memory
	scoreMatrix = (double **)safeMalloc(sizeof(double *)*M);
	scoreMatrix[0] = (double *)safeMalloc(sizeof(double) * (M * N));
	for (i = 1; i < M; i++)
		scoreMatrix[i] = scoreMatrix[0] + N * i;

	// initialize scoreMatrix
	for (i = 0; i < N; i++)
		scoreMatrix[0][i] = 0;
	for (i = 0; i < M; i++)
		scoreMatrix[i][0] = 0;

	// matrix score
	for (i = 1; i < M; i++)
	{
		for (j = 1; j < N; j++)
		{
			tmp[0] = scoreMatrix[i - 1][j - 1] + scorePair(seq1[i - 1], seq2[j - 1]); /* sequence index is i-1 and j-1*/
			tmp[1] = scoreMatrix[i - 1][j] + MISMATCH;
			tmp[2] = scoreMatrix[i][j - 1] + MISMATCH;
			tmp[3] = 0;
			scoreMatrix[i][j] = max(tmp, 4);
			if (scoreMatrix[i][j] >= bestScore)
			{
				bestI = i;
				bestJ = j;
				bestScore = scoreMatrix[i][j];
			}
		}
	}
	// trace back
	allocSpace = M > N ? (M + extSpace) : (N + extSpace);
	qSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
	mSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
	tSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));

	i = bestI;
	j = bestJ;
	k = 0;
	for (;;)
	{
		if (i == 0 && j == 0) break;
		if (scoreMatrix[i][j] <= 0)
		{
			traceI = i;
			traceJ = j;
			break;
		}
		if (k >= allocSpace)
		{
			qSeq = (char *)safeRealloc(qSeq, k + extSpace);
			mSeq = (char *)safeRealloc(mSeq, k + extSpace);
			tSeq = (char *)safeRealloc(tSeq, k + extSpace);
			allocSpace = k;
		}
		if (i >= 1 && j >= 1 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + scorePair(seq1[i - 1], seq2[j - 1]))
		{
			qSeq[k] = seq1[i - 1];
			tSeq[k] = seq2[j - 1];
			mSeq[k] = getPairChar(qSeq[k], tSeq[k]);
			i--;
			j--;
		}
		else if (i >= 1 && j >= 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + MISMATCH)
		{
			qSeq[k] = seq1[i - 1];
			mSeq[k] = '.';
			tSeq[k] = '-';
			i--;
		}
		else if (j >= 1 && i >= 0 && scoreMatrix[i][j] == scoreMatrix[i][j - 1] + MISMATCH)
		{
			qSeq[k] = '-';
			mSeq[k] = '.';
			tSeq[k] = seq2[j - 1];
			j--;
		}
		k++;
	}

	qSeq[k] = '\0';
	mSeq[k] = '\0';
	tSeq[k] = '\0';

	swAlign->mirStart = i;
	swAlign->mirEnd = bestI;
	swAlign->tarStart = j;
	swAlign->tarEnd = bestJ;

	int rna1LeftLen = rna1Len - bestI;
	int rna2LeftLen = rna2Len - bestJ;
	int rna1RightLen = i;
	int rna2RightLen = j;
	int leftLen  = MAX(rna1LeftLen, rna2LeftLen);
	int rightLen = MAX(rna1RightLen, rna2RightLen);

	char *qlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));
	char *mlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));
	char *tlSeq = (char *)safeMalloc(sizeof(char) * (leftLen + 1));

	char *qrSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));
	char *mrSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));
	char *trSeq = (char *)safeMalloc(sizeof(char) * (rightLen + 1));

	for (i = 0; i < leftLen; i++)
	{
		qlSeq[i] = '-';
		mlSeq[i] = '.';
		tlSeq[i] = '-';
	}
	for (i = 0; i < rightLen; i++)
	{
		qrSeq[i] = '-';
		mrSeq[i] = '.';
		trSeq[i] = '-';
	}
	for (i = rna1LeftLen, j = leftLen; i > 0 && j > 0; i--, j--)
	{
		qlSeq[j - 1] = rna1[i - 1];
	}
	for (i = bestJ, j = leftLen; i < rna2Len && j > 0; i++, j--)
	{
		tlSeq[j - 1] = rna2[i];
	}
	for (i = 0, j = rna1Len - rna1RightLen; i < rna1RightLen && j < rna1Len; i++, j++)
	{
		qrSeq[i] = rna1[j];
	}
	for (i = 0, j = rna2RightLen; i < rna2RightLen && j > 0; i++, j--)
	{
		trSeq[i] = rna2[j - 1];
	}

	qlSeq[leftLen] = '\0';
	mlSeq[leftLen] = '\0';
	tlSeq[leftLen] = '\0';
	qrSeq[rightLen] = '\0';
	mrSeq[rightLen] = '\0';
	trSeq[rightLen] = '\0';

	int pairLen = strlen(qSeq);
	int extLen = pairLen + leftLen + rightLen;
	char *nqSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));
	char *nmSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));
	char *ntSeq = (char *)safeMalloc(sizeof(char) * (extLen + 1));

	strcpy(nqSeq, qlSeq);
	strcpy(nmSeq, mlSeq);
	strcpy(ntSeq, tlSeq);
	strcat(nqSeq, qSeq);
	strcat(nmSeq, mSeq);
	strcat(ntSeq, tSeq);
	strcat(nqSeq, qrSeq);
	strcat(nmSeq, mrSeq);
	strcat(ntSeq, trSeq);

	nqSeq[extLen] = '\0';
	nmSeq[extLen] = '\0';
	ntSeq[extLen] = '\0';

	swAlign->mirSeq = nqSeq;
	swAlign->tarSeq = ntSeq;
	swAlign->pairStr = nmSeq;
	swAlign->maxPairNum = maxPairNum(swAlign->pairStr);
	swAlign->alignScore = bestScore;
	freeFloatMatrix(scoreMatrix);
	safeFree(newRna1);
	safeFree(qSeq);
	safeFree(tSeq);
	safeFree(mSeq);
	safeFree(qlSeq);
	safeFree(tlSeq);
	safeFree(mlSeq);
	safeFree(qrSeq);
	safeFree(trSeq);
	safeFree(mrSeq);
	return bestScore;
}

int maxPairNum(char *pairStr)
{
	int maxNum = 0;
	int i = 0;
	int pairNum = 0;
	int pairLen = strlen(pairStr);
	for (i = 0; i < pairLen; i++)
	{
		if (pairStr[i] == '|')
		{
			pairNum++;
		}
		else {
			if (pairNum > maxNum) maxNum = pairNum;
			pairNum = 0;
		}
	}
	return maxNum;
}
