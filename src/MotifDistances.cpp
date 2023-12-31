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

#include <iostream>
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

SEXP motifDistances (SEXP  cc, SEXP align, SEXP top, SEXP go, SEXP ge, SEXP inputPWM, SEXP inputScores)
{
	SEXP  returnData;
	PlatformSupport* Plat = new PlatformSupport();
    ColumnComp* CC;
	Alignment* ALIGN=NULL;
    bool colChosen=false, alignChosen=false;
	int matchTopX = TOP_MATCH;
	
	//Default alignment settings
	double gapOpen = DFLT_GAP_OPEN;
	double gapExtend = DFLT_GAP_EXTEND;
	bool overlapAlign = DFLT_OVLP_ALIGN;
	bool extendOverlap=false;
	bool ungapped=false;	
	
	PROTECT (returnData=allocVector(VECSXP,1));
		
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
	
	gapOpen=NUMERIC_VALUE(go);
	gapExtend=NUMERIC_VALUE(ge); 
	matchTopX=INTEGER_VALUE(top); 	
	
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
	
	
	////////////////////////////////////////////////////////////////////////////////////
	//////// Main Program /////////////////////////////////////////////////////////////
	//Initialise the background
	Plat->ReadBackground();
	//Read in the matrices
	Plat->ReadTransfacFile(  inputPWM, 0);
		Plat->ReadScoreDists(inputScores);
	if(Plat->GetMatCount()>1)
	{}	

	Plat->PreAlign(ALIGN);		
	SET_VECTOR_ELT(returnData,0,Plat-> PrintPairwise());
	
	delete(CC);
	delete(ALIGN);
	
	delete(Plat);
	UNPROTECT(1);
	return   VECTOR_ELT(returnData,0);
}

//R code
extern "C" {
	void RmotifDistances (SEXP  cc, SEXP align, SEXP top, SEXP go, SEXP ge,  SEXP inputPWM, SEXP inputScores)
	{	motifDistances(cc, align, top, go, ge, inputPWM, inputScores );	}
}