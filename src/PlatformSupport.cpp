//////////////////////////////////////////////////////////////////////////////////
//
// STAMP version 1 
//
// Written By: Shaun Mahony
//
// PlatformSupport.cpp
//
// Started: 1st Nov 2005
//
// Copyright 2007 Shaun Mahony
//
// This file is part of STAMP.
//
// STAMP is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// STAMP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with STAMP; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
////////////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include "PlatformSupport.h"
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R.h>

//Constructor
PlatformSupport::PlatformSupport()
{
	matCount=0; 
	matchDBSize=0;
	markov=NULL; 
	charMap=NULL; 
	scoreDistStdDev=NULL;
	scoreDistMean = NULL;
	pairwiseAlign=NULL;
	backgroundOrder=0;
	usingWeighting=false;
	
	charMap=(char***)malloc((MAX_MARKOV+1)*sizeof(char**));
	markov=(double**)malloc((MAX_MARKOV+1)*sizeof(double*));
	
	for(int j=1;j<=MAX_MARKOV;j++)
	{
		markov[j]=(double*)malloc(((int)pow(B,j))*sizeof(double));
		charMap[j]=(char**)malloc(((int)pow(B,j))*sizeof(char*));
		for(int i=0;i<pow(B,j);i++) 
			charMap[j][i]=(char*)malloc((j+1)*sizeof(char));
	}
	backgroundSet=false;
}

//Read in a Markov background
void PlatformSupport::ReadBackground(char* fn)
{
	char *xmer;
	int i,j=0;
	double px;
	
	xmer=(char*)malloc((MAX_MARKOV+1)*sizeof(char));	
	if(fn!=NULL){
		FILE* modelp = fopen(fn, "r");
		if(modelp==NULL){perror("Cannot open background file");}
		
		while(fscanf(modelp,"%d %s %lf\n",&i,xmer,&px) != EOF)
		{
			j=strlen(xmer);
			strcpy(charMap[j][i], xmer);
			markov[j][i]=px;
		}
		backgroundOrder=j;
		
		fclose(modelp);
	}else{
		//Set background to neutral bias
		backgroundOrder=1;
		markov[1][0]=.25; markov[1][1]=.25; markov[1][2]=.25; markov[1][3]=.25;
		strcpy(charMap[1][0], "A"); strcpy(charMap[1][1], "C"); strcpy(charMap[1][2], "G"); strcpy(charMap[1][3], "T");
	}
	
	backgroundSet=true;
	free(xmer);
}

//Read in a TransFac file
int PlatformSupport::ReadTransfacFile( SEXP inputPWM, SEXP inputDB)
{
	int i, j, p;
	double curr_ttl;
	Motif** currMotifs; int currCnt;	
	SEXP inputName, input;
	PROTECT(inputName=NEW_CHARACTER(50));
	
	if(inputPWM!=0)
	{
		currMotifs=inputMotifs;
		currCnt=matCount;
		input=inputPWM;
	}
	else if (inputDB!=0)
	{
		currMotifs=matchMotifs;
		currCnt=matchDBSize;
		input=inputDB;
	}
	else
	{Rprintf("\tERROR.\n");}
	
	if(!backgroundSet)
	{	Rprintf("\tReadBackground not called; exiting");
	}
	
	currCnt=0;
	inputName=GET_NAMES(input);
	for (currCnt=0; currCnt< length(input) ;currCnt++)
	{			
		currMotifs[currCnt] = new Motif(length(VECTOR_ELT(input,currCnt))/4);
		strcpy(currMotifs[currCnt]->name,CHAR(STRING_ELT(AS_CHARACTER(inputName),currCnt)) );
		currMotifs[currCnt]->weighting=1;
		p=0;
		for (i=0; i<length(VECTOR_ELT(input,currCnt))/4; i++)
		{				
			if(NULL!=input)
			{
				curr_ttl=0;
				for(j=0; j<B; j++)
				{							
					currMotifs[currCnt]->n[i][j] = NUMERIC_DATA(VECTOR_ELT(input,currCnt))[p]; 
					curr_ttl+=currMotifs[currCnt]->n[i][j];
					p++;
				}
				
				for (j=0; j<B; j++)
				{
					currMotifs[currCnt]->f[i][j] = (currMotifs[currCnt]->n[i][j] + (SCALE_FACTOR*markov[1][j]))/(curr_ttl+SCALE_FACTOR);
					currMotifs[currCnt]->pwm[i][j] = log_2(currMotifs[currCnt]->f[i][j]/markov[1][j]);						
				}
			}
		}
	}
	
	if(inputPWM!=0)
	{
		currCnt=currCnt;
		matCount=currCnt;
	}
	else
	{matchDBSize=currCnt;}
	
	UNPROTECT(1);
	return (currCnt);
}

//Find the scores distribution in the input Motifs
SEXP PlatformSupport::GetRandDistrib( Alignment* A_man)
{
	int i, j;
	int x, y;
	int align1, align2, aLen;
	bool forward1, forward2;
	double bestScore;
	int len=400;
	int compt=0;
	SEXP score;	
	PROTECT (score =allocMatrix(REALSXP, len,7));
		
	//Set up the mean and std_dev arrays
	double** sum = (double**)malloc(sizeof(double*)*maxLen);
	double** max = (double**)malloc(sizeof(double*)*maxLen);
	double** min = (double**)malloc(sizeof(double*)*maxLen);
	double** std_dev = (double**)malloc(sizeof(double*)*maxLen);
	double** count = (double**)malloc(sizeof(double*)*maxLen);
	double** sampSq = (double**)malloc(sizeof(double*)*maxLen);
	for(i=0; i<maxLen; i++) {
		sum[i] = (double*)malloc(sizeof(double)*maxLen);
		max[i] = (double*)malloc(sizeof(double)*maxLen);
		min[i] = (double*)malloc(sizeof(double)*maxLen);
		std_dev[i] = (double*)malloc(sizeof(double)*maxLen);
		count[i] = (double*)malloc(sizeof(double)*maxLen);
		sampSq[i] = (double*)malloc(sizeof(double)*maxLen);
		for(j=0; j<maxLen; j++) {
			max[i][j]=0;
			min[i][j]=100000;
			sum[i][j]=0;
			std_dev[i][j]=0;
			count[i][j]=0;
			sampSq[i][j]=0;
		}
	}
	
	//Compare each matrix to every other matrix
	Rprintf( "\tGenerate scores :\n");
	for(i=0; i<matCount; i++){
		for(j=0; j<i; j++) 
		{
			if(i!=j)
			{
				bestScore = A_man->AlignMotifs2D(inputMotifs[i], inputMotifs[j], align1, align2, aLen, forward1, forward2);
				
				//Add bestScore to the proper mean and std_dev array
				x=inputMotifs[i]->len;
				if(x<minLen){x=minLen;}
				else if(x>=maxLen){x=maxLen-1;}
				y=inputMotifs[j]->len;
				if(y<minLen){y=minLen;}
				else if(y>=maxLen){y=maxLen-1;}
				sum[x][y]+=bestScore;
				sum[y][x]+=bestScore;
				sampSq[x][y] += (bestScore*bestScore);
				sampSq[y][x] += (bestScore*bestScore);
				count[x][y]++;
				count[y][x]++;
				if(bestScore>max[x][y])
					max[x][y]=bestScore;
				else if(bestScore<min[x][y])
					min[x][y]=bestScore;
				if(bestScore>max[y][x])
					max[y][x]=bestScore;
				else if(bestScore<min[y][x])
					min[y][x]=bestScore;
			}
		}
		if (((i+1)%250)==0)
		{	Rprintf("\t\t%d scores generated\n",i+1);}
	}
	
	for(x=minLen; x<maxLen; x++)
		for(y=minLen; y<maxLen; y++)
		{	
			std_dev[x][y] = sampSq[x][y] -((sum[x][y]*sum[x][y])/count[x][y]);
			std_dev[x][y] = std_dev[x][y]/count[x][y];
			if(std_dev[x][y]!=0)
				std_dev[x][y] = sqrt(std_dev[x][y]);
		}
	
	
	for(x=minLen; x<maxLen; x++)
		for(y=minLen; y<maxLen; y++)
		{	if(count[x][y]>0)
		{
			DOUBLE_DATA(score)[compt]=x;
			DOUBLE_DATA(score)[compt+len]=y;
			DOUBLE_DATA(score)[compt+2*len]=sum[x][y]/count[x][y];
			DOUBLE_DATA(score)[compt+3*len]=std_dev[x][y];
			DOUBLE_DATA(score)[compt+4*len]= count[x][y];
			DOUBLE_DATA(score)[compt+5*len]= min[x][y];
			DOUBLE_DATA(score)[compt+6*len]=max[x][y];
			compt++;			
		}
		else		
		{
			DOUBLE_DATA(score)[compt]=x;
			DOUBLE_DATA(score)[compt+len]=y;
			DOUBLE_DATA(score)[compt+2*len]=0;
			DOUBLE_DATA(score)[compt+3*len]=0;
			DOUBLE_DATA(score)[compt+4*len]= 0;
			DOUBLE_DATA(score)[compt+5*len]=0;
			DOUBLE_DATA(score)[compt+6*len]=0;
			compt++;
		}
		}	
	for(i=0; i<maxLen; i++) {
		free(sum[i]);
		free(std_dev[i]);
		free(count[i]);
		free(sampSq[i]);
	}
	free(sum);
	free(std_dev);
	free(count);
	free(sampSq);
	UNPROTECT(1);
	return score;
}

//Read in the score distance distributions
void PlatformSupport::ReadScoreDists( SEXP inputScores)
{
	int i,j,k;
	int x, y;
	
	//set up the matrices
	scoreDistMean = (double**)malloc(sizeof(double*)*maxLen);
	scoreDistMax = (double**)malloc(sizeof(double*)*maxLen);
	scoreDistMin = (double**)malloc(sizeof(double*)*maxLen);
	scoreDistStdDev = (double**)malloc(sizeof(double*)*maxLen);
	for(i=0; i<maxLen; i++) {
		scoreDistMean[i] = (double*)malloc(sizeof(double)*maxLen);
		scoreDistMax[i] = (double*)malloc(sizeof(double)*maxLen);
		scoreDistMin[i] = (double*)malloc(sizeof(double)*maxLen);
		scoreDistStdDev[i] = (double*)malloc(sizeof(double)*maxLen);
		for(j=0; j<maxLen; j++) {
			scoreDistMax[i][j]=0;
			scoreDistMin[i][j]=0;
			scoreDistMean[i][j]=0;
			scoreDistStdDev[i][j]=0;
		}
	}
	
	for (k=0; k<400 ; k++ )
	{
		x=REAL(inputScores)[k] ;
		y=REAL(inputScores)[k+400] ;
		scoreDistMean[x][y] =REAL(inputScores)[k+2*400] ;
		scoreDistStdDev[x][y]=REAL(inputScores)[k+3*400] ;
		scoreDistMax[x][y] =REAL(inputScores)[k+6*400] ;
		scoreDistMin[x][y] =REAL(inputScores)[k+5*400] ;
	}
}

//Convert a Score to a Z-score given two lengths.
double PlatformSupport::Score2ZScore(int len1, int len2, double score)
{
	int l1=len1, l2=len2;
	double mean, std_dev;
	
	if(len1<minLen)
		l1=minLen;
	else if(len1>maxLen-1)
		l1=maxLen-1;
	if(len2<minLen)
		l2=minLen;
	else if(len2>maxLen-1)
		l2=maxLen-1;
	
	mean = scoreDistMean[l1][l2];
	std_dev=scoreDistStdDev[l1][l2];
	if(std_dev<=0)
		std_dev=1;
	
	return((score-mean)/std_dev);
}


//Private Method: Convert a Score to a P-value given two lengths.
double PlatformSupport::Score2PVal(int len1, int len2, double score)
{
	int l1=len1, l2=len2;
	double mean, std_dev;
	double start;
	double p_val=0;
	
	if(len1<minLen)
		l1=minLen;
	else if(len1>maxLen-1)
		l1=maxLen-1;
	if(len2<minLen)
		l2=minLen;
	else if(len2>maxLen-1)
		l2=maxLen-1;
	
	mean = scoreDistMean[l1][l2];
	std_dev=scoreDistStdDev[l1][l2];
	if(std_dev<=0)
		std_dev=1;
	
	start = score - mean;
	p_val=pnorm(start,0,std_dev,TRUE,FALSE);
	return(p_val);
}

//Private Method: Convert a Score to an approximate distance measure (from Feng & Doolittle)
double PlatformSupport::Score2Dist(int len1, int len2, double score, double maxScore)
{
	int l1=len1, l2=len2;
	double S_eff=0, S_rand;
	
	if(len1<minLen)
		l1=minLen;
	else if(len1>maxLen-1)
		l1=maxLen-1;
	if(len2<minLen)
		l2=minLen;
	else if(len2>maxLen-1)
		l2=maxLen-1;
	
	double std_dev=scoreDistStdDev[l1][l2];
	if(std_dev<=0)
		std_dev=1;
	S_rand = scoreDistMean[l1][l2] - 4*std_dev;
	S_rand = scoreDistMin[l1][l2];
	S_eff = (score - S_rand)/(maxScore - S_rand);
	
	if(S_eff <= 0)
		S_eff=-1*(log(0.001));
	else
		S_eff = -1*(log(S_eff));
	
	return((S_eff));
}



//Align all matrices against all others

void PlatformSupport::PreAlign(Alignment* A_man)
{
	int i, j;
	int i1, i2, aL;
	double curr_score, curr_z_score, curr_p_val, max_score;
	bool forward1, forward2;
	pairwiseAlign = new AlignRec*[matCount];
	for(i=0; i<matCount; i++)
	{	
		pairwiseAlign[i] = new AlignRec[matCount];
	}
	//firstly align each matrix to itself (to get max scores later)
	for(i=0; i<matCount; i++)
	{
		curr_score = A_man->AlignMotifs(inputMotifs[i], inputMotifs[i], i1, i2, aL, forward1);
		pairwiseAlign[i][i].forward1=forward1;
		pairwiseAlign[i][i].forward2=false;
		pairwiseAlign[i][i].i1=i1;
		pairwiseAlign[i][i].i2=i2;
		pairwiseAlign[i][i].score=curr_score;
		curr_z_score = Score2ZScore(inputMotifs[i]->len, inputMotifs[i]->len, curr_score);
		pairwiseAlign[i][i].z_score=curr_z_score;
		curr_p_val = Score2PVal(inputMotifs[i]->len, inputMotifs[i]->len, curr_score);
		pairwiseAlign[i][i].p_value = curr_p_val;
		pairwiseAlign[i][i].CopyAlignSec(A_man->alignSection, aL);
		strcpy(pairwiseAlign[i][i].alignedNames[0], inputMotifs[i]->name);
		strcpy(pairwiseAlign[i][i].alignedNames[1], inputMotifs[i]->name);
		pairwiseAlign[i][i].alignedIDs[0] = i; pairwiseAlign[i][i].alignedIDs[1] = i;
	}

	//Go through each possible set of alignments
	for(i=0; i<matCount; i++)
		for(j=0; j<matCount; j++){//Now possible to set j<i, as long as [i][j] and [j][i] are copied
			if(i!=j){
				curr_score = A_man->AlignMotifs2D(inputMotifs[i], inputMotifs[j], i1, i2, aL, forward1, forward2);
				pairwiseAlign[i][j].forward1=forward1;
				pairwiseAlign[i][j].forward2=forward2;
				pairwiseAlign[i][j].i1=i1;
				pairwiseAlign[i][j].i2=i2;
				pairwiseAlign[i][j].score=curr_score;
				curr_z_score = Score2ZScore(inputMotifs[i]->len, inputMotifs[j]->len, curr_score);
				pairwiseAlign[i][j].z_score=curr_z_score;
				curr_p_val = Score2PVal(inputMotifs[i]->len, inputMotifs[j]->len, curr_score);
				pairwiseAlign[i][j].p_value = curr_p_val;
				pairwiseAlign[i][j].CopyAlignSec(A_man->alignSection, aL);
				strcpy(pairwiseAlign[i][j].alignedNames[0], inputMotifs[i]->name);
				strcpy(pairwiseAlign[i][j].alignedNames[1], inputMotifs[j]->name);
				pairwiseAlign[i][i].alignedIDs[0] = i; pairwiseAlign[i][i].alignedIDs[1] = j;
				max_score = (pairwiseAlign[i][i].score + pairwiseAlign[j][j].score)/2;
				pairwiseAlign[i][j].dist = -1 * log(pairwiseAlign[i][j].p_value);
			}
		}
}

//Print out the pairwise alignments
SEXP PlatformSupport::PrintPairwise()
{
	int i, j;
	SEXP Eval2;
	PROTECT (Eval2 =allocMatrix(REALSXP, matCount,matCount));
	int compt=0;
		
	for(i=0; i<matCount; i++){
		for(j=0; j<matCount; j++){
			if(i!=j){
				double Eval = 1-pairwiseAlign[i][j].p_value;
				DOUBLE_DATA(Eval2)[compt]=Eval;
			}else
				{
				DOUBLE_DATA(Eval2)[compt]=0;
				}
		compt++;		
		}
	}
	
	UNPROTECT(1);
	return Eval2;
}

//Find the best matching motifs in the match set and print the pairs to a file
SEXP PlatformSupport::SimilarityMatching(Alignment* A_man, const int matchTopX)
{
	const char *tmp_strandSeq[1000];
	const char *tmp_strandMatch[1000];
	int compt=0;
	double* topScores;
	int* topIndices;
	bool printAll=false;
	int i, j, x, y;
	double currScore, currPVal;
	Motif* one; Motif* two;
	int i1, i2, aL;
	bool forward1, forward2;
	char currName[STR_LEN];
	char*** topAligns;
	int topX = matchTopX;	
	SEXP tf, eval, seq, match, strandSeq, strandMatch, pwm, name, vec;
	int n=GetMatCount()*topX;
		
	PROTECT (tf=NEW_CHARACTER(n));
	PROTECT (eval=NEW_NUMERIC(n));
	PROTECT (seq=NEW_CHARACTER(n));
	PROTECT (match=NEW_CHARACTER(n));
	PROTECT (strandSeq=NEW_CHARACTER(n));
	PROTECT (strandMatch=NEW_CHARACTER(n));
	PROTECT (pwm=NEW_LIST(n));
	PROTECT(name=NEW_CHARACTER(GetMatCount()));
	
	
	if(topX>GetMatchDBSize()){
		topX=GetMatchDBSize();
	}
	Rprintf("\tMotif matches : %d\n",topX);
	bool inserted=false;
	topScores = new double[topX];
	topIndices = new int[topX];
	topAligns = new char**[topX];
	
	//init
	for(x=0; x<topX; x++){
		topScores[x]=0; topIndices[x]=0; 
		topAligns[x]=new char*[2];
		topAligns[x][0]=new char[STR_LEN];
		topAligns[x][1]=new char[STR_LEN];
		strcpy(topAligns[x][0], "");strcpy(topAligns[x][1], "");
	}
		
	if(printAll){
		Rprintf("\t\t");
		for(j=0; j<GetMatchDBSize(); j++){
			Rprintf("\t%s\t", matchMotifs[j]->GetName());
		}
		Rprintf("\t\n");
	}
		
	for(i=0; i<GetMatCount(); i++){
		
		if(printAll)
			Rprintf("\t%s\t", inputMotifs[i]->GetName());
		
		for(x=0; x<topX; x++){
			topScores[x]=0; topIndices[x]=0; 
			strcpy(topAligns[x][0], "");strcpy(topAligns[x][1], "");
		}
		
		
		for(j=0; j<GetMatchDBSize(); j++)
		{			
			currScore = A_man->AlignMotifs2D(inputMotifs[i], matchMotifs[j], i1, i2, aL, forward1, forward2);
			currPVal = Score2PVal(inputMotifs[i]->len, matchMotifs[j]->len, currScore);
			if(printAll){Rprintf("\t%lf\t", currPVal);/*Rprintf("\t%s\t%lf\n", matchMotifs[j]->GetName(),currPVal);*/}
			
			
			//Check the current score against the topScores
			inserted=false;
			for(x=0; x<topX && !inserted; x++){
			
				if(currPVal>topScores[x]){
									
					//Shift and insert
					for(y=topX-1; y>x; y--){
						topScores[y]=topScores[y-1];
						topIndices[y]=topIndices[y-1];
						strcpy(topAligns[y][0],topAligns[y-1][0]);
						strcpy(topAligns[y][1],topAligns[y-1][1]);
						tmp_strandSeq[y]=tmp_strandSeq[y-1];
						tmp_strandMatch[y]=tmp_strandMatch[y-1];
					}
					topScores[x] = currPVal;
					topIndices[x]=j;
				
					if(forward1){
						one = inputMotifs[i];
						tmp_strandMatch[x]="+";
					}else{
						one = new Motif(inputMotifs[i]->GetLen());
						inputMotifs[i]->RevCompMotif(one);
						 tmp_strandMatch[x]="-";
					}
					if(forward2){
						two = matchMotifs[j];
						tmp_strandSeq[x]="+";
					}else{
						two = new Motif(matchMotifs[j]->GetLen());
						matchMotifs[j]->RevCompMotif(two);
						tmp_strandSeq[x]="-";
					}
										
					A_man->CopyAlignmentConsensus(one, two,topAligns[x][0], topAligns[x][1]);
					if(!forward1){
						delete one;
						
					}if(!forward2){
						delete two;
					}
					inserted=true;
				}
			}
		}
		if(printAll){Rprintf("\t\n");}
						
		SET_STRING_ELT(name, i, mkChar( inputMotifs[i]->GetName()));
		
		for(x=0; x<topX; x++)
		{
			sprintf(currName, "%s", matchMotifs[topIndices[x]]->GetName());
			double Eval =1-topScores[x];
			SET_STRING_ELT(tf, compt, mkChar(currName));
			DOUBLE_DATA(eval)[compt]=Eval;
			SET_STRING_ELT(seq, compt, mkChar( topAligns[x][0]));
			SET_STRING_ELT(match, compt, mkChar( topAligns[x][1]));
			SET_STRING_ELT(strandSeq, compt, mkChar(tmp_strandSeq[x]));
			SET_STRING_ELT(strandMatch, compt, mkChar(tmp_strandMatch[x]));
			SET_VECTOR_ELT(pwm,compt,matchMotifs[topIndices[x]]->PrintMotif(NULL));
			compt++;
		}
	}
	
	delete [] topScores;
	delete [] topIndices;
	for(x=0; x<topX; x++)
	{
		delete [] topAligns[x][0];delete [] topAligns[x][1];
		delete [] topAligns[x];
	}
	delete [] topAligns;
	
	PROTECT(vec=allocVector(VECSXP,8));
	SET_VECTOR_ELT(vec,0,name);
	SET_VECTOR_ELT(vec,1,tf);
	SET_VECTOR_ELT(vec,2,pwm);
	SET_VECTOR_ELT(vec,3,eval);
	SET_VECTOR_ELT(vec,4,seq);
	SET_VECTOR_ELT(vec,5,match);
	SET_VECTOR_ELT(vec,6,strandSeq);
	SET_VECTOR_ELT(vec,7,strandMatch);
	
	UNPROTECT(9);
	return vec;	
}

////////////////////////////////////////////////////////////////////////////////////
///////////  Motif Helpers  ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//Convert f to n (multiply by 100)
void PlatformSupport::f_to_n(Motif* m)
{
	for(int i=0; i<m->GetLen(); i++){
		for(int j=0; j<B; j++){
			m->n[i][j] = floor(m->f[i][j]*DFLT_NUM_INSTANCES);
		}
	}
}
//Convert n to pwm
void PlatformSupport::n_to_pwm(Motif* m)
{
	int i,j;
	double ttl;
	for(i=0; i<m->GetLen(); i++){
		ttl=0;
		for(j=0; j<B; j++){ttl += m->n[i][j];}
		for(j=0; j<B; j++)
			m->pwm[i][j] = log_2(((m->n[i][j] + (SCALE_FACTOR*markov[1][j]))/(ttl+SCALE_FACTOR))/markov[1][j]);
	}
}

//Information content
double PlatformSupport::InfoContent(Motif* m)
{
	double sum=0.0;
	
	for(int j=0;j<m->GetLen();j++) {
		for(int b=0;b<B;b++) {
			if(m->f[j][b]) {
				sum+=m->f[j][b]*log_2(m->f[j][b]);
			}
		}
	}
	return 2+sum;
}


//Log base 2
double PlatformSupport::log_2(double x)
{ return(log(x) / LOG_2);}


//Destructor
PlatformSupport::~PlatformSupport()
{
	int i, j;
	
	if(markov!=NULL && charMap!=NULL){
		for(i=1; i<=MAX_MARKOV; i++){
			for(j=0;j<pow(B,i);j++) 
			{	free(charMap[i][j]); }
			free(charMap[i]);
			free(markov[i]);
		}
		free(charMap);
		free(markov);
	}
	if(scoreDistMean!=NULL){
		for(i=0; i<maxLen; i++) 
			free(scoreDistMean[i]);
		free(scoreDistMean);
	}
	if(scoreDistStdDev!=NULL){
		for(i=0; i<maxLen; i++) 
			free(scoreDistStdDev[i]);
		free(scoreDistStdDev);
	}
	if(pairwiseAlign!=NULL){
		for(i=0; i<matCount; i++)
			delete [] pairwiseAlign[i];
		delete [] pairwiseAlign;
	}
	for(i=0; i<matCount; i++)
		delete inputMotifs[i];
}
