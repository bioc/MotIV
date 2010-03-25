 
//////////////////////////////////////////////////////////////////////////////////
//
// STAMP version 1
//
// Written By: Shaun Mahony
//
// main.cpp
//
// Started: 31st Oct 2005
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

#include <stdio.h>
#include <stdlib.h>
#include "Alignment.h"
#include "ColumnComp.h"
#include "PlatformSupport.h"
#include "RandPSSMGen.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

SEXP generateScoresDB(SEXP cc, SEXP align,SEXP go, SEXP ge, SEXP nRand,   SEXP inputPWM)
{
	SEXP listRandPWM, returnData;
	PROTECT (returnData=allocVector(VECSXP,1));
	PROTECT(listRandPWM=allocVector(VECSXP,1));
	PlatformSupport* Plat = new PlatformSupport();
    ColumnComp* CC;
	Alignment* ALIGN;
	int i;
	// int argc=length(argv);
	// int nRand;

    bool colChosen=false, alignChosen=false;
	//int matchTopX = TOP_MATCH;
	
	//Misc option settings
	char* randMatOut = new char[STR_LEN];
	//Default alignment settings
	double gapOpen = DFLT_GAP_OPEN;
	double gapExtend = DFLT_GAP_EXTEND;
	bool overlapAlign = DFLT_OVLP_ALIGN;
	bool extendOverlap=false;
	bool ungapped=false;
	
	
	
				if((strcmp(CHAR(STRING_ELT(cc,0)), "PCC"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "pcc"))==0){
					CC = new PearsonCorrelation(); //Pearson's correllation coefficient
				}else if((strcmp(CHAR(STRING_ELT(cc,0)), "ALLR"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "allr"))==0){
					CC = new ALLR(); //ALLR
				}else if((strcmp(CHAR(STRING_ELT(cc,0)), "ALLR_LL"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "allr_ll"))==0){
					CC = new ALLR_LL(); //ALLR with lower limit
				}else if((strcmp(CHAR(STRING_ELT(cc,0)), "CS"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "cs"))==0){
					CC = new ChiSq(); //Pearson's Chi Square
				}else if((strcmp(CHAR(STRING_ELT(cc,0)), "KL"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "kl"))==0){
					CC = new KullbackLieber(); //Kullback-Lieber
				}else if((strcmp(CHAR(STRING_ELT(cc,0)), "SSD"))==0 || (strcmp(CHAR(STRING_ELT(cc,0)), "ssd"))==0){
					CC = new SumSqDiff(); //sum of squared difference
				}else{
					CC = new PearsonCorrelation(); //Default = PCC
				}
				colChosen=true;
	
				gapOpen=NUMERIC_VALUE(go); //strtod(CHAR(STRING_ELT(go,0)), NULL);
				gapExtend=NUMERIC_VALUE(ge); //strtod(CHAR(STRING_ELT(ge,0)), NULL);
	
	
	
				if((strcmp(CHAR(STRING_ELT(align,0)), "NW"))==0 || (strcmp(CHAR(STRING_ELT(align,0)), "nw"))==0)
				{
					ALIGN = new NeedlemanWunsch(CC, gapOpen, gapExtend, overlapAlign, extendOverlap);
				}
				if((strcmp(CHAR(STRING_ELT(align,0)), "SWU"))==0 || (strcmp(CHAR(STRING_ELT(align,0)), "swu"))==0)
				{
					ALIGN = new SmithWatermanUngappedExtended(CC); ungapped=true;
				}
				if((strcmp(CHAR(STRING_ELT(align,0)), "SWA"))==0 || (strcmp(CHAR(STRING_ELT(align,0)), "swa"))==0)
				{
					ALIGN = new SmithWatermanAffine(CC, gapOpen, gapExtend, overlapAlign, extendOverlap);
				}
				if((strcmp(CHAR(STRING_ELT(align,0)), "SW"))==0 || (strcmp(CHAR(STRING_ELT(align,0)), "sw"))==0)
				{
					ALIGN = new SmithWaterman(CC, gapOpen, gapExtend, overlapAlign, extendOverlap);
				}
				alignChosen = true;
	
	
	
			//	nRand=INTEGER_VALUE(nRand);
	
	
	
	
	
	
	
////////////////////////////////////////////////////////////////////////////////////
//////// Main Program /////////////////////////////////////////////////////////////

		//Initialise the background
		Plat->ReadBackground();
		//Read in the matrices
		Plat->ReadTransfacFile(  inputPWM, 0);
		if(ungapped)
		Rprintf("\n\tUngapped Alignment\n");
		else
		Rprintf("\tGap open = %.3lf, gap extend = %.3lf\n", gapOpen, gapExtend);

		//Generate some random matrices
		RandPSSMGen* RPG = new RandPSSMGen(Plat->inputMotifs, Plat->GetMatCount(), INTEGER_VALUE(nRand), randMatOut);
		SET_VECTOR_ELT(listRandPWM,0,RPG->RunGenerator());		
		Plat->ReadTransfacFile( VECTOR_ELT(listRandPWM,0),0);
		SET_VECTOR_ELT(returnData,0,Plat->GetRandDistrib( ALIGN));
		delete(CC);
		delete(ALIGN);

	delete(Plat);	
	UNPROTECT(2);
	return  VECTOR_ELT(returnData,0);
}

//R code
extern "C" {
	 SEXP RgenerateScoresDB (SEXP cc, SEXP align,SEXP go, SEXP ge, SEXP nRand, SEXP inputPWM)
		{	generateScoresDB( cc, align, go, ge, nRand, inputPWM );}
} //extern "C"
