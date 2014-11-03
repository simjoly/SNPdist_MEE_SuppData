//
//  anal.h
//  
//
//  Created by Simon Joly on 13-11-29.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef _anal_h
#define _anal_h


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "nexusdata.h"
#include "organismsdata.h"
#include "GenFunctions.h"
#include <algorithm>

#include <string.h>


#define	MAX_BUFFER_SIZE	20000000
#define MAX 1000
#define punct "()[]{}/\\,;:=*\'\"`+"			/* these are NEXUS definitions */
#define NL_punct "()[]{}/\\,;:=*\'\"`+\r\n"             /*                             */

#define isNL_NEXUSwhiteSpace(c)  ( strchr(" \t\v\f", (c)) || (((c) <= 6) && ((c) >= 0)))
#define isNL_NEXUSpunct(c) ( strchr(NL_punct,(c)) )
#define isNEXUSpunct(c) ( strchr(punct,(c)) )
#define isNEXUSwhiteSpace(c)	( isspace((c)) )
/* current NEXUS format does not exclude ASCII 0-6 */
#define NULL_RETURN {TempPtr=""; return(TempPtr);}

#define NB_OF_A_INDIVIDUALS 2
#define NB_OF_B_INDIVIDUALS 1
#define NB_OF_C_INDIVIDUALS 1
#define NB_OF_D_INDIVIDUALS 1
#define NB_OF_E_INDIVIDUALS 1
#define NB_OF_F_INDIVIDUALS 1
#define NB_OF_G_INDIVIDUALS 1


using namespace std;


/*** Function declaration ***/


int get_position_nexus(string the_string);
string char_consensus(char char1, char char2);
void read_organisms(string name);
void doOrgDimensions(void);
void doOrgLabels(void);
//string	To_Uppercase(string input_string);
void AttributeSequencesToOrganisms(void);
double get_PBC_distance(int org1, int org2);
double get_MIN_distance(int org1, int org2);


/* Private functions */

static string nextToken(void);
static int parse_assignment2(string target);
static void doUnrecognizedBlock(void);

information info_data;
nexusdata nexus_data;
organismsdata org_data;

#endif
