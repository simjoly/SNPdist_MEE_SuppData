/**********************************************\
|                                              |
|                  Analyses                    | 
|                                              |
|                version 1.00                  |
|                                              |
|----------------------------------------------|
|                                              |
|           released June 3rd, 2006            |
|                                              |
|      (c) Copyright 2006 by Simon Joly        |
|                                              |
\**********************************************/

/*
    Pofad is a free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2.

    Pofad is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <regex>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "anal.h"


unsigned int i;
int j,k,z;                     //integrers used for "for" loops
int samplesize;                //Number of generation sampled
int nb_p, nb_c, nb_i, nb_e, nb_s, nb_l, nb_h;  //number of sequences for each species
ofstream logfile, DistributionFile;
vector<string> FileList;
double within_dist, between_sister_dist, between_clade_dist;
double dist_a_a, dist_a_b, dist_a_c, dist_a_d, dist_a_e, dist_a_f, dist_a_g;
vector<double> genpofad_mean_a_h, genpofad_mean_h_b, genpofad_mean_a_a, genpofad_mean_a_b, genpofad_mean_a_c, genpofad_mean_a_d, genpofad_mean_a_e, genpofad_mean_a_f, genpofad_mean_a_g;
vector<double> mrca_mean_a_h, mrca_mean_h_b, mrca_mean_a_a, mrca_mean_a_b, mrca_mean_a_c, mrca_mean_a_d, mrca_mean_a_e, mrca_mean_a_f, mrca_mean_a_g;
vector<double> machstates_mean_a_h, machstates_mean_h_b, machstates_mean_a_a, machstates_mean_a_b, machstates_mean_a_c, machstates_mean_a_d, machstates_mean_a_e, machstates_mean_a_f, machstates_mean_a_g;
vector<double> FRQ_mean_a_h, FRQ_mean_h_b, FRQ_mean_a_a, FRQ_mean_a_b, FRQ_mean_a_c, FRQ_mean_a_d, FRQ_mean_a_e, FRQ_mean_a_f, FRQ_mean_a_g;
vector<double> nei_mean_a_h, nei_mean_h_b, nei_mean_a_a, nei_mean_a_b, nei_mean_a_c, nei_mean_a_d, nei_mean_a_e, nei_mean_a_f, nei_mean_a_g;
vector<double> PBC_mean_a_h, PBC_mean_h_b, PBC_mean_a_a, PBC_mean_a_b, PBC_mean_a_c, PBC_mean_a_d, PBC_mean_a_e, PBC_mean_a_f, PBC_mean_a_g;
vector<double> MIN_mean_a_h, MIN_mean_h_b, MIN_mean_a_a, MIN_mean_a_b, MIN_mean_a_c, MIN_mean_a_d, MIN_mean_a_e, MIN_mean_a_f, MIN_mean_a_g;
vector<double> ISP_mean_a_h, ISP_mean_h_b, ISP_mean_a_a, ISP_mean_a_b, ISP_mean_a_c, ISP_mean_a_d, ISP_mean_a_e, ISP_mean_a_f, ISP_mean_a_g;
int number_replicates;
double Divergence[7] = {0,0.002,0.004,0.004,0.008,0.012,0.02};


/***** Local variables *****/

static char *bufPtr;
static string aTokenPtr;
static string  LocalToken;

double sumdist;
int thecount;


/***** Main program *****/

int main(int argc, char **argv)
{
    string input_file="";                   //A string use to store the name of the file containning the dataset names
    string output_file="";
    string temp;
    int sim_chars=0;
    bool equal_hybrid=true;
    bool hybrid_analysis=false;

    char *cur;
    char ch;
    char character,selection;
	int w;
	
    for(w=1;w<argc;w++)
    {
        cur=argv[w];
        ch=cur[0];
        
        if (ch == '-')
        {
            ch=toupper(cur[1]);
            switch(ch)
            {
                case 'H':
                    hybrid_analysis=true;
                    w++;
                    if ( w<argc && (argv[w][0] != '-')) {
                        if (atoi(argv[w]) == 1) {
                            equal_hybrid=false;
                        }
                    }
                    else w--;
                case 'C':
                    w++;
                    if ( w<argc && (argv[w][0] != '-')) {
                        sim_chars=atoi(argv[w]);
                        org_data.GetSimChar(atoi(argv[w]));
                        nexus_data.GetSimChar(atoi(argv[w]));
                    }
                    else w--;
                    continue;
                case 'I':
                    w++;
		            if(w<argc && argv[w][0] != '-') {
                        input_file = argv[w];
                    }
		            else {w--;}
		            continue;
                case 'O':
                    w++;
		            if(w<argc && argv[w][0] != '-')
                    {
                        output_file = argv[w];
                    }
		            else {w--;}
		            continue;
                default:
                    continue;
            }
        }
    }

    cout << endl;
    cout << " -------------------------------"<< endl;
    cout << " Welcome in Analysis version 1.0"<< endl;
    cout << " (c) Simon Joly, 2006-2014"<< endl;
    cout << " -------------------------------"<< endl << endl;


    if(input_file=="") {
        cerr << "No input file entered (-i)... exiting!" << endl;
        exit(0);
    }
    if(sim_chars==0) {
        cerr << "You need to enter the umber of characters that were simulated (-c)... exiting!" << endl;
        exit(0);
    }
/*
    cout << endl << " Enter the name of the file containning the names of the simulated datasets: ";
    cin >> input_file;
*/
    
    //Get output file fame
    if (output_file=="") {
        int lastdot = input_file.find_last_of("_");
        if (lastdot == -1) lastdot = input_file.size();
        int lastbar = input_file.find_last_of("/")+1;
        if (lastbar == -1) lastbar = 0;
        output_file = input_file.substr(lastbar, lastdot-lastbar) + "_results.txt";
        if (output_file == "") output_file = "results.txt";
    }
/*
    cout << endl << " Enter the total number of characters simulated: ";
    cin >> i;
    org_data.GetSimChar(i);
    nexus_data.GetSimChar(i);
*/
	read_organisms("organisms.nex");
    info_data.SetIsIgnoreMissingData(true);

    cout<< " Reading dataset names in file " << input_file << endl;
    ifstream infile;
    infile.open(input_file.c_str());
    while (!infile) {
        cout << "Error in opening file " << input_file << endl;
        exit(0);
        }    

    infile >> temp;
    while (!infile.eof()) {
        FileList.push_back(temp);
        infile >> temp;
    }
    
    infile.close();
    number_replicates = FileList.size();
    cout << " Found " << number_replicates << " simulated datasets" << endl;
    cout << " Assuming " << sim_chars << " characters were simulated" << endl;

    if (!hybrid_analysis) {
        cout << " Analysing datasets:" << endl;
        
        for (z=0; z<number_replicates; z++) {
            cout << "  -> Performing replicate: " << (z+1) << "\r";
            cout.flush();
            
            /*** Read File ***/
            
            vector<string> tokens;
            infile.open(FileList[z].c_str());
            stringstream ss;
            
            nexus_data.InitializeMember();
            nexus_data.GetNTaxa(org_data.ReturnNOrg()*4);
            nexus_data.InitDataMatrix();
            org_data.InitDataMatrix();
            int orgnumber=0;
            
            while (!infile.eof()) {
                getline(infile, temp, '\n');
                if (temp=="") continue;
                string buf; // buffer string
                ss.str (temp);
                while (ss >> buf)
                    tokens.push_back(buf);
                regex rr("(.+_.+)");
                if (regex_match(tokens[0], rr) ) {
                    nexus_data.ReadTaxa(tokens[0]);
                    nexus_data.AddCharacter(orgnumber, tokens[2]);
                    orgnumber++;
                }
                tokens.clear();
                ss.clear();
            }
            
            temp = nexus_data.ReturnSequence(0);
            nexus_data.GetNChars(temp.size());
            org_data.GetNChar(temp.size());
            
            string char1="";
            string char2="";
            string char3="";
            string char4="";
            string astring, astring_1, astring_2, astring_3, astring_4;
            for (i=0; i<org_data.ReturnNOrg(); i++)
	        {
                astring = To_Uppercase(org_data.ReturnOrganism(i));
                astring_1 = To_Uppercase(org_data.ReturnOrganism(i)) + "1";
                astring_2 = To_Uppercase(org_data.ReturnOrganism(i)) + "2";
                astring_3 = To_Uppercase(org_data.ReturnOrganism(i)) + "3";
                astring_4 = To_Uppercase(org_data.ReturnOrganism(i)) + "4";
                if (get_position_nexus(astring_3) != -1) /* simulation of polyploids */
	        	{
                    for (j=0; j<nexus_data.ReturnNChars(); j++)
		            {
                        char1 = nexus_data.ReturnChar( get_position_nexus(astring_1), j);
                        char2 = nexus_data.ReturnChar( get_position_nexus(astring_2), j);
                        char3 = nexus_data.ReturnChar( get_position_nexus(astring_3), j);
                        char4 = nexus_data.ReturnChar( get_position_nexus(astring_4), j);
                        temp = char_consensus(char1[0],char2[0]);
                        temp = char_consensus(temp[0],char3[0]);
                        temp = char_consensus(temp[0],char4[0]);
                        org_data.AddCharacter(i,temp);
		            }
	        	}
                else
	        	{
                    for (j=0; j<nexus_data.ReturnNChars(); j++)
		            {
                        char1 = nexus_data.ReturnChar( get_position_nexus(astring_1), j);
                        char2 = nexus_data.ReturnChar( get_position_nexus(astring_2), j);
                        temp = char_consensus(char1[0],char2[0]);
                        org_data.AddCharacter(i,temp);
		            }
	        	}
	        }
            
            //For methods with FRQ matrix
            AttributeSequencesToOrganisms();
            org_data.InitFRQMatrix();
            org_data.BuildFRQMatrix();
            //For distance based methods
            nexus_data.InitMatrix();
            nexus_data.ComputeDistanceMatrix();
            
            
            /***  Calculate Genpofad distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("genpofad");
            
            // A vs A
            
            genpofad_mean_a_a.push_back(org_data.ReturnDist(0,1));
            
            // A vs B
            
            genpofad_mean_a_b.push_back(org_data.ReturnDist(0,2));
            
            
            // A vs C
            
            genpofad_mean_a_c.push_back(org_data.ReturnDist(0,3));
            
            // A vs D
            
            genpofad_mean_a_d.push_back(org_data.ReturnDist(0,4));
            
            // A vs E
            
            genpofad_mean_a_e.push_back(org_data.ReturnDist(0,5));
            
            // A vs F
            
            genpofad_mean_a_f.push_back(org_data.ReturnDist(0,6));
            
            // A vs G
            
            genpofad_mean_a_g.push_back(org_data.ReturnDist(0,7));
            
            // Calculate within species, between sister, and between clade distances
            /*
             within_dist = (dist_a_a + dist_b_b + dist_c_c + dist_d_d) / 4;
             between_sister_dist = (dist_a_b + dist_c_d) / 2;
             between_clade_dist = (dist_a_c + dist_a_d + dist_b_c + dist_b_d) / 4;
             
             genpofad_within_dist_mean.push_back(within_dist);
             genpofad_betweensister_dist_mean.push_back(between_sister_dist);
             genpofad_between_clade_dist_mean.push_back(between_clade_dist);
             */
            
            
            /***  Calculate mrca distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("mrca");
            
            // A vs A
            mrca_mean_a_a.push_back(org_data.ReturnDist(0,1));
            
            // A vs B
            mrca_mean_a_b.push_back(org_data.ReturnDist(0,2));
            
            // A vs C
            mrca_mean_a_c.push_back(org_data.ReturnDist(0,3));
            
            // A vs D
            mrca_mean_a_d.push_back(org_data.ReturnDist(0,4));
            
            // A vs E
            mrca_mean_a_e.push_back(org_data.ReturnDist(0,5));
            
            // A vs F
            mrca_mean_a_f.push_back(org_data.ReturnDist(0,6));
            
            // A vs G
            mrca_mean_a_g.push_back(org_data.ReturnDist(0,7));
            
            
            /***  Calculate matchstates distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("matchstates");
            
            // A vs A
            machstates_mean_a_a.push_back(org_data.ReturnDist(0,1));
            
            // A vs B
            machstates_mean_a_b.push_back(org_data.ReturnDist(0,2));
            
            // A vs C
            machstates_mean_a_c.push_back(org_data.ReturnDist(0,3));
            
            // A vs D
            machstates_mean_a_d.push_back(org_data.ReturnDist(0,4));
            
            // A vs E
            machstates_mean_a_e.push_back(org_data.ReturnDist(0,5));
            
            // A vs F
            machstates_mean_a_f.push_back(org_data.ReturnDist(0,6));
            
            // A vs G
            machstates_mean_a_g.push_back(org_data.ReturnDist(0,7));
            
            
            
            /***  Calculate 2ISP distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("2ISP");
            
            // A vs A
            ISP_mean_a_a.push_back(org_data.ReturnDist(0,1));
            
            // A vs B
            ISP_mean_a_b.push_back(org_data.ReturnDist(0,2));
            
            // A vs C
            ISP_mean_a_c.push_back(org_data.ReturnDist(0,3));
            
            // A vs D
            ISP_mean_a_d.push_back(org_data.ReturnDist(0,4));
            
            // A vs E
            ISP_mean_a_e.push_back(org_data.ReturnDist(0,5));
            
            // A vs F
            ISP_mean_a_f.push_back(org_data.ReturnDist(0,6));
            
            // A vs G
            ISP_mean_a_g.push_back(org_data.ReturnDist(0,7));
            
            
            
            /***  Calculate FRQ distances  ***/
            /*
             org_data.InitDistMatrix();
             org_data.CalculateDistances("FRQ");
             
             // A vs A
             FRQ_mean_a_a.push_back(org_data.ReturnDist(0,1));
             
             // A vs B
             FRQ_mean_a_b.push_back(org_data.ReturnDist(0,2));
             
             // A vs C
             FRQ_mean_a_c.push_back(org_data.ReturnDist(0,3));
             
             // A vs D
             FRQ_mean_a_d.push_back(org_data.ReturnDist(0,4));
             
             // A vs E
             FRQ_mean_a_e.push_back(org_data.ReturnDist(0,5));
             
             // A vs F
             FRQ_mean_a_f.push_back(org_data.ReturnDist(0,6));
             
             // A vs G
             FRQ_mean_a_g.push_back(org_data.ReturnDist(0,7));
             */
            
            
            /***  Calculate Nei's genetic distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("nei");
            
            // A vs A
            nei_mean_a_a.push_back(org_data.ReturnDist(0,1));
            
            // A vs B
            nei_mean_a_b.push_back(org_data.ReturnDist(0,2));
            
            // A vs C
            nei_mean_a_c.push_back(org_data.ReturnDist(0,3));
            
            // A vs D
            nei_mean_a_d.push_back(org_data.ReturnDist(0,4));
            
            // A vs E
            nei_mean_a_e.push_back(org_data.ReturnDist(0,5));
            
            // A vs F
            nei_mean_a_f.push_back(org_data.ReturnDist(0,6));
            
            // A vs G
            nei_mean_a_g.push_back(org_data.ReturnDist(0,7));
            
            
            
            /***  Calculate PBC distances  ***/
            
            // A vs A
            PBC_mean_a_a.push_back(get_PBC_distance(0,1));
            
            // A vs B
            PBC_mean_a_b.push_back(get_PBC_distance(0,2));
            
            // A vs C
            PBC_mean_a_c.push_back(get_PBC_distance(0,3));
            
            // A vs D
            PBC_mean_a_d.push_back(get_PBC_distance(0,4));
            
            // A vs E
            PBC_mean_a_e.push_back(get_PBC_distance(0,5));
            
            // A vs F
            PBC_mean_a_f.push_back(get_PBC_distance(0,6));
            
            // A vs G
            PBC_mean_a_g.push_back(get_PBC_distance(0,7));
            
            
            
            /***  Calculate MIN distances  ***/
            
            // A vs A
            MIN_mean_a_a.push_back(get_MIN_distance(0,1));
            
            // A vs B
            MIN_mean_a_b.push_back(get_MIN_distance(0,2));
            
            // A vs C
            MIN_mean_a_c.push_back(get_MIN_distance(0,3));
            
            // A vs D
            MIN_mean_a_d.push_back(get_MIN_distance(0,4));
            
            // A vs E
            MIN_mean_a_e.push_back(get_MIN_distance(0,5));
            
            // A vs F
            MIN_mean_a_f.push_back(get_MIN_distance(0,6));
            
            // A vs G
            MIN_mean_a_g.push_back(get_MIN_distance(0,7));
            
            //		org_data.DeleteCharMatrix();
            //		nexus_data.DeleteCharMatrix();
            infile.close();
        }
        
        
        
        
        //Write distance values to a file in column format for drawing histograms
        
        cout << endl;
        cout << " Writing results in file \"" << output_file << "\"" << endl;
        
        DistributionFile.open(output_file.c_str());
        while (!DistributionFile) {
            cout << "Error openning distribution file!" << endl;
            exit(0);
        }
        DistributionFile << "distance" << '\t' << "method" << '\t' << "divergence" << endl;
        
        for (i=0;i<number_replicates;i++) {
            for (j=0;j<7;j++) {
                DistributionFile << genpofad_mean_a_a[i] << '\t' << "GENPOFAD" << '\t' << Divergence[0] << endl;
                DistributionFile << genpofad_mean_a_b[i] << '\t' << "GENPOFAD" << '\t' << Divergence[1] << endl;
                DistributionFile << genpofad_mean_a_c[i] << '\t' << "GENPOFAD" << '\t' << Divergence[2] << endl;
                //            DistributionFile << genpofad_mean_a_d[i] << '\t' << "GENPOFAD" << '\t' << Divergence[3] << endl;
                DistributionFile << genpofad_mean_a_e[i] << '\t' << "GENPOFAD" << '\t' << Divergence[4] << endl;
                DistributionFile << genpofad_mean_a_f[i] << '\t' << "GENPOFAD" << '\t' << Divergence[5] << endl;
                DistributionFile << genpofad_mean_a_g[i] << '\t' << "GENPOFAD" << '\t' << Divergence[6] << endl;
                DistributionFile << machstates_mean_a_a[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[0] << endl;
                DistributionFile << machstates_mean_a_b[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[1] << endl;
                DistributionFile << machstates_mean_a_c[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[2] << endl;
                //            DistributionFile << machstates_mean_a_d[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[3] << endl;
                DistributionFile << machstates_mean_a_e[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[4] << endl;
                DistributionFile << machstates_mean_a_f[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[5] << endl;
                DistributionFile << machstates_mean_a_g[i] << '\t' << "MATCHSTATES" << '\t' << Divergence[6] << endl;
                DistributionFile << mrca_mean_a_a[i] << '\t' << "MRCA" << '\t' << Divergence[0] << endl;
                DistributionFile << mrca_mean_a_b[i] << '\t' << "MRCA" << '\t' << Divergence[1] << endl;
                DistributionFile << mrca_mean_a_c[i] << '\t' << "MRCA" << '\t' << Divergence[2] << endl;
                //            DistributionFile << mrca_mean_a_d[i] << '\t' << "MRCA" << '\t' << Divergence[3] << endl;
                DistributionFile << mrca_mean_a_e[i] << '\t' << "MRCA" << '\t' << Divergence[4] << endl;
                DistributionFile << mrca_mean_a_f[i] << '\t' << "MRCA" << '\t' << Divergence[5] << endl;
                DistributionFile << mrca_mean_a_g[i] << '\t' << "MRCA" << '\t' << Divergence[6] << endl;
                /*
                 DistributionFile << FRQ_mean_a_a[i] << '\t' << "FRQ" << '\t' << Divergence[0] << endl;
                 DistributionFile << FRQ_mean_a_b[i] << '\t' << "FRQ" << '\t' << Divergence[1] << endl;
                 DistributionFile << FRQ_mean_a_c[i] << '\t' << "FRQ" << '\t' << Divergence[2] << endl;
                 //DistributionFile << FRQ_mean_a_d[i] << '\t' << "FRQ" << '\t' << Divergence[3] << endl;
                 DistributionFile << FRQ_mean_a_e[i] << '\t' << "FRQ" << '\t' << Divergence[4] << endl;
                 DistributionFile << FRQ_mean_a_f[i] << '\t' << "FRQ" << '\t' << Divergence[5] << endl;
                 DistributionFile << FRQ_mean_a_g[i] << '\t' << "FRQ" << '\t' << Divergence[6] << endl;
                 */
                DistributionFile << nei_mean_a_a[i] << '\t' << "NEI" << '\t' << Divergence[0] << endl;
                DistributionFile << nei_mean_a_b[i] << '\t' << "NEI" << '\t' << Divergence[1] << endl;
                DistributionFile << nei_mean_a_c[i] << '\t' << "NEI" << '\t' << Divergence[2] << endl;
                //DistributionFile << nei_mean_a_d[i] << '\t' << "NEI" << '\t' << Divergence[3] << endl;
                DistributionFile << nei_mean_a_e[i] << '\t' << "NEI" << '\t' << Divergence[4] << endl;
                DistributionFile << nei_mean_a_f[i] << '\t' << "NEI" << '\t' << Divergence[5] << endl;
                DistributionFile << nei_mean_a_g[i] << '\t' << "NEI" << '\t' << Divergence[6] << endl;
                DistributionFile << PBC_mean_a_a[i] << '\t' << "PBC" << '\t' << Divergence[0] << endl;
                DistributionFile << PBC_mean_a_b[i] << '\t' << "PBC" << '\t' << Divergence[1] << endl;
                DistributionFile << PBC_mean_a_c[i] << '\t' << "PBC" << '\t' << Divergence[2] << endl;
                //DistributionFile << PBC_mean_a_d[i] << '\t' << "PBC" << '\t' << Divergence[3] << endl;
                DistributionFile << PBC_mean_a_e[i] << '\t' << "PBC" << '\t' << Divergence[4] << endl;
                DistributionFile << PBC_mean_a_f[i] << '\t' << "PBC" << '\t' << Divergence[5] << endl;
                DistributionFile << PBC_mean_a_g[i] << '\t' << "PBC" << '\t' << Divergence[6] << endl;
                DistributionFile << MIN_mean_a_a[i] << '\t' << "MIN" << '\t' << Divergence[0] << endl;
                DistributionFile << MIN_mean_a_b[i] << '\t' << "MIN" << '\t' << Divergence[1] << endl;
                DistributionFile << MIN_mean_a_c[i] << '\t' << "MIN" << '\t' << Divergence[2] << endl;
                //DistributionFile << MIN_mean_a_d[i] << '\t' << "MIN" << '\t' << Divergence[3] << endl;
                DistributionFile << MIN_mean_a_e[i] << '\t' << "MIN" << '\t' << Divergence[4] << endl;
                DistributionFile << MIN_mean_a_f[i] << '\t' << "MIN" << '\t' << Divergence[5] << endl;
                DistributionFile << MIN_mean_a_g[i] << '\t' << "MIN" << '\t' << Divergence[6] << endl;
                DistributionFile << ISP_mean_a_a[i] << '\t' << "PP" << '\t' << Divergence[0] << endl;
                DistributionFile << ISP_mean_a_b[i] << '\t' << "PP" << '\t' << Divergence[1] << endl;
                DistributionFile << ISP_mean_a_c[i] << '\t' << "PP" << '\t' << Divergence[2] << endl;
                //DistributionFile << ISP_mean_a_d[i] << '\t' << "PP" << '\t' << Divergence[3] << endl;
                DistributionFile << ISP_mean_a_e[i] << '\t' << "PP" << '\t' << Divergence[4] << endl;
                DistributionFile << ISP_mean_a_f[i] << '\t' << "PP" << '\t' << Divergence[5] << endl;
                DistributionFile << ISP_mean_a_g[i] << '\t' << "PP" << '\t' << Divergence[6] << endl;
            }
        }
        DistributionFile.close();
    } // end: if (!hybrid_analysis)

    
    else {  //If hybrid simulations
        
        cout << " Analysis of ";
        if (equal_hybrid) cout << "equal";
        else cout << "unequal";
        cout << " hybrid simulations" << endl;
        cout << " Analysing datasets:" << endl;
        for (z=0; z<number_replicates; z++) {
            cout << "  -> Performing replicate: " << (z+1) << "\r";
            cout.flush();
            
            /*** Read File ***/
            
            vector<string> tokens;
            infile.open(FileList[z].c_str());
            stringstream ss;
            
            nexus_data.InitializeMember();
            if (!equal_hybrid) nexus_data.GetNTaxa((org_data.ReturnNOrg()*4)+2);
            if (equal_hybrid) nexus_data.GetNTaxa((org_data.ReturnNOrg()+1)*4);
            nexus_data.InitDataMatrix();
            org_data.InitDataMatrix();
            int orgnumber=0;
            string temp2;
            
            while (!infile.eof()) {
                getline(infile, temp, '\n');
                if (temp=="") continue;
                string buf; // buffer string
                ss.str (temp);
                while (ss >> buf)
                    tokens.push_back(buf);
                regex rr("(.+_.+)");
                if (regex_match(tokens[0], rr) ) {
                    temp2 = tokens[0].substr(0, tokens[0].find_first_of("_"));
                    //cout << "temp2 = " << temp2;
                    if (temp2=="4") {
                        if (equal_hybrid || (tokens[0]!="4_3" && tokens[0]!="4_4")) {
                            nexus_data.ReadTaxa(tokens[0]);
                            nexus_data.AddCharacter(orgnumber, tokens[2]);
                            orgnumber++;
                            tokens[0]="3_"+tokens[0];
                        }
                    }
                    //cout << "  |  tokens[0] = " << tokens[0] << endl;
                    nexus_data.ReadTaxa(tokens[0]);
                    nexus_data.AddCharacter(orgnumber, tokens[2]);
                    orgnumber++;
                }
                tokens.clear();
                ss.clear();
            }
            
            temp = nexus_data.ReturnSequence(0);
            nexus_data.GetNChars(temp.size());
            org_data.GetNChar(temp.size());
            
            string char1="";
            string char2="";
            string char3="";
            string char4="";
            string char5="";
            string char6="";
            string char7="";
            string char8="";
            string astring, astring_1, astring_2, astring_3, astring_4,astring_5, astring_6, astring_7,astring_8;
            for (i=0; i<org_data.ReturnNOrg(); i++) {
                if (i == 3 || i==5 || i==6 || i==7) {
                    continue;
                }
                astring = To_Uppercase(org_data.ReturnOrganism(i));
                astring_1 = To_Uppercase(org_data.ReturnOrganism(i)) + "1";
                astring_2 = To_Uppercase(org_data.ReturnOrganism(i)) + "2";
                astring_3 = To_Uppercase(org_data.ReturnOrganism(i)) + "3";
                astring_4 = To_Uppercase(org_data.ReturnOrganism(i)) + "4";
                astring_5 = To_Uppercase(org_data.ReturnOrganism(i)) + "4_1";
                astring_6 = To_Uppercase(org_data.ReturnOrganism(i)) + "4_2";
                astring_7 = To_Uppercase(org_data.ReturnOrganism(i)) + "4_3";
                astring_8 = To_Uppercase(org_data.ReturnOrganism(i)) + "4_4";
                if (get_position_nexus(astring_3) != -1) /* simulation of polyploids */
                {
                    for (j=0; j<nexus_data.ReturnNChars(); j++) {
                        char1 = nexus_data.ReturnChar( get_position_nexus(astring_1), j);
                        char2 = nexus_data.ReturnChar( get_position_nexus(astring_2), j);
                        char3 = nexus_data.ReturnChar( get_position_nexus(astring_3), j);
                        char4 = nexus_data.ReturnChar( get_position_nexus(astring_4), j);
                        temp = char_consensus(char1[0],char2[0]);
                        temp = char_consensus(temp[0],char3[0]);
                        temp = char_consensus(temp[0],char4[0]);
                        if (get_position_nexus(astring_5) == -1) {
                            org_data.AddCharacter(i,temp);
                            continue;
                        }
                        char5 = nexus_data.ReturnChar( get_position_nexus(astring_5), j);
                        char6 = nexus_data.ReturnChar( get_position_nexus(astring_6), j);
                        temp = char_consensus(temp[0],char5[0]);
                        temp = char_consensus(temp[0],char6[0]);
                        if (equal_hybrid) {  // if hybrid has constitution 2:1...
                            char7 = nexus_data.ReturnChar( get_position_nexus(astring_7), j);
                            char8 = nexus_data.ReturnChar( get_position_nexus(astring_8), j);
                            temp = char_consensus(temp[0],char7[0]);
                            temp = char_consensus(temp[0],char8[0]);
                        }
                        org_data.AddCharacter(i,temp);
                    }
                }
                else
                {
                    for (j=0; j<nexus_data.ReturnNChars(); j++)
                    {
                        char1 = nexus_data.ReturnChar( get_position_nexus(astring_1), j);
                        char2 = nexus_data.ReturnChar( get_position_nexus(astring_2), j);
                        temp = char_consensus(char1[0],char2[0]);
                        org_data.AddCharacter(i,temp);
                    }
                }
            }
            
            
            //For methods with FRQ matrix
            AttributeSequencesToOrganisms();
            org_data.InitFRQMatrix();
            org_data.BuildFRQMatrix();
            //For distance based methods
            nexus_data.InitMatrix();
            nexus_data.ComputeDistanceMatrix();
            
            
            /***  Calculate Genpofad distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("genpofad");
            
            // A vs H
            genpofad_mean_a_h.push_back(org_data.ReturnDist(0,2));
            
            // C vs H
            genpofad_mean_h_b.push_back(org_data.ReturnDist(2,4));
            
            
            
            /***  Calculate mrca distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("mrca");
            
            // A vs H
            mrca_mean_a_h.push_back(org_data.ReturnDist(0,2));
            
            // H vs C
            mrca_mean_h_b.push_back(org_data.ReturnDist(2,4));
            
            
            
            /***  Calculate matchstates distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("matchstates");
            
            // A vs H
            machstates_mean_a_h.push_back(org_data.ReturnDist(0,2));
            
            // H vs B
            machstates_mean_h_b.push_back(org_data.ReturnDist(2,4));
            
            
            
            /***  Calculate 2ISP distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("2ISP");
            
            // A vs H
            ISP_mean_a_h.push_back(org_data.ReturnDist(0,2));
            
            // H vs B
            ISP_mean_h_b.push_back(org_data.ReturnDist(2,4));
            
            
            /***  Calculate Nei's genetic distances  ***/
            
            org_data.InitDistMatrix();
            org_data.CalculateDistances("nei");
            
            // A vs H
            nei_mean_a_h.push_back(org_data.ReturnDist(0,2));
            
            // H vs B
            nei_mean_h_b.push_back(org_data.ReturnDist(2,4));
            
            
            /***  Calculate PBC distances  ***/
            
            // A vs H
            PBC_mean_a_h.push_back(get_PBC_distance(0,2));
            
            // H vs B
            PBC_mean_h_b.push_back(get_PBC_distance(2,4));
            
            
            /***  Calculate MIN distances  ***/
            
            // A vs H
            MIN_mean_a_h.push_back(get_MIN_distance(0,2));
            
            // H vs B
            MIN_mean_h_b.push_back(get_MIN_distance(2,4));
            
            
            infile.close();
        }
        
        //Write distance values to a file in column format for drawing histograms
        
        cout << endl;
        cout << " Writing results in file \"" << output_file << "\"" << endl;
        
        DistributionFile.open(output_file.c_str());
        while (!DistributionFile) {
            cout << "Error openning distribution file!" << endl;
            exit(0);
        }

        DistributionFile << "index" << '\t' << "method" << endl;
        
        for (i=0;i<number_replicates;i++) {
            for (j=0;j<7;j++) {
                DistributionFile << (genpofad_mean_a_h[i]/(genpofad_mean_a_h[i]+genpofad_mean_h_b[i])) << '\t' << "GENPOFAD" << endl;
                DistributionFile << (machstates_mean_a_h[i]/(machstates_mean_a_h[i]+machstates_mean_h_b[i])) << '\t' << "MATCHSTATES" << endl;
                DistributionFile << (mrca_mean_a_h[i]/(mrca_mean_a_h[i]+mrca_mean_h_b[i])) << '\t' << "MRCA" << endl;
                DistributionFile << (ISP_mean_a_h[i]/(ISP_mean_a_h[i]+ISP_mean_h_b[i])) << '\t' << "PP" << endl;
                //DistributionFile << (FRQ_mean_a_h[i]/(FRQ_mean_a_h[i]+FRQ_mean_h_b[i])) << '\t' << "FRQ" << endl;
                DistributionFile << (nei_mean_a_h[i]/(nei_mean_a_h[i]+nei_mean_h_b[i])) << '\t' << "NEI" << endl;
                DistributionFile << (PBC_mean_a_h[i]/(PBC_mean_a_h[i]+PBC_mean_h_b[i])) << '\t' << "PBC" << endl;
                DistributionFile << (MIN_mean_a_h[i]/(MIN_mean_a_h[i]+MIN_mean_h_b[i])) << '\t' << "MIN" << endl;
            }
        }
        DistributionFile.close();

    } //end: if(hybrid_analysis)
    
    cout << endl;
    return 0;

}





/*********************************************

Function : char_consensus

Returns the consensus base of two nucleic acids

Variables: - two nucleic acids


*********************************************/


string char_consensus(char char1, char char2)
{

    string firstchar;
    string secondchar;

    firstchar = toupper(char1);
    secondchar = toupper(char2);

    // If both characters are identicals
    if (firstchar == secondchar) return firstchar;

    // Handling 0-1 characters
    if ( (firstchar == "0"  && secondchar == "1") || (firstchar == "1"  && secondchar == "0") ) return "?";
    if (firstchar == "0"  && secondchar == "0") return "0";
    if (firstchar == "1"  && secondchar == "1") return "1";

    // Handling of missing characters
    if (firstchar == "?" ) return secondchar;
    if (secondchar == "?" ) return firstchar;

    // Handling of gaps -> can only treat them as missing data!
    if (firstchar == "-" ) return secondchar;
    if (secondchar == "-" ) return firstchar;

    if ( (firstchar == "A"  && secondchar == "G") || (firstchar == "G"  && secondchar == "A") ) return "R";
    if ( (firstchar == "A"  && secondchar == "C") || (firstchar == "C"  && secondchar == "A") ) return "M";
    if ( (firstchar == "A"  && secondchar == "T") || (firstchar == "T"  && secondchar == "A") ) return "W";
    if ( (firstchar == "G"  && secondchar == "C") || (firstchar == "C"  && secondchar == "G") ) return "S";
    if ( (firstchar == "G"  && secondchar == "T") || (firstchar == "T"  && secondchar == "G") ) return "K";
    if ( (firstchar == "C"  && secondchar == "T") || (firstchar == "T"  && secondchar == "C") ) return "Y";

    if ( (firstchar == "A"  && secondchar == "R") || (firstchar == "R"  && secondchar == "A") ) return "R";
    if ( (firstchar == "C"  && secondchar == "R") || (firstchar == "R"  && secondchar == "C") ) return "V";
    if ( (firstchar == "G"  && secondchar == "R") || (firstchar == "R"  && secondchar == "G") ) return "R";
    if ( (firstchar == "T"  && secondchar == "R") || (firstchar == "R"  && secondchar == "T") ) return "D";

    if ( (firstchar == "A"  && secondchar == "M") || (firstchar == "M"  && secondchar == "A") ) return "M";
    if ( (firstchar == "C"  && secondchar == "M") || (firstchar == "M"  && secondchar == "C") ) return "M";
    if ( (firstchar == "G"  && secondchar == "M") || (firstchar == "M"  && secondchar == "G") ) return "V";
    if ( (firstchar == "T"  && secondchar == "M") || (firstchar == "M"  && secondchar == "T") ) return "H";

    if ( (firstchar == "A"  && secondchar == "W") || (firstchar == "W"  && secondchar == "A") ) return "W";
    if ( (firstchar == "C"  && secondchar == "W") || (firstchar == "W"  && secondchar == "C") ) return "H";
    if ( (firstchar == "G"  && secondchar == "W") || (firstchar == "W"  && secondchar == "G") ) return "D";
    if ( (firstchar == "T"  && secondchar == "W") || (firstchar == "W"  && secondchar == "T") ) return "W";

    if ( (firstchar == "A"  && secondchar == "S") || (firstchar == "S"  && secondchar == "A") ) return "V";
    if ( (firstchar == "C"  && secondchar == "S") || (firstchar == "S"  && secondchar == "C") ) return "S";
    if ( (firstchar == "G"  && secondchar == "S") || (firstchar == "S"  && secondchar == "G") ) return "S";
    if ( (firstchar == "T"  && secondchar == "S") || (firstchar == "S"  && secondchar == "T") ) return "B";

    if ( (firstchar == "A"  && secondchar == "K") || (firstchar == "K"  && secondchar == "A") ) return "D";
    if ( (firstchar == "C"  && secondchar == "K") || (firstchar == "K"  && secondchar == "C") ) return "B";
    if ( (firstchar == "G"  && secondchar == "K") || (firstchar == "K"  && secondchar == "G") ) return "K";
    if ( (firstchar == "T"  && secondchar == "K") || (firstchar == "K"  && secondchar == "T") ) return "K";

    if ( (firstchar == "A"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "A") ) return "H";
    if ( (firstchar == "C"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "C") ) return "Y";
    if ( (firstchar == "G"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "G") ) return "B";
    if ( (firstchar == "T"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "T") ) return "Y";

    if ( (firstchar == "A"  && secondchar == "B") || (firstchar == "B"  && secondchar == "A") ) return "N";
    if ( (firstchar == "C"  && secondchar == "B") || (firstchar == "B"  && secondchar == "C") ) return "B";
    if ( (firstchar == "G"  && secondchar == "B") || (firstchar == "B"  && secondchar == "G") ) return "B";
    if ( (firstchar == "T"  && secondchar == "B") || (firstchar == "B"  && secondchar == "T") ) return "B";

    if ( (firstchar == "A"  && secondchar == "D") || (firstchar == "D"  && secondchar == "A") ) return "D";
    if ( (firstchar == "C"  && secondchar == "D") || (firstchar == "D"  && secondchar == "C") ) return "N";
    if ( (firstchar == "G"  && secondchar == "D") || (firstchar == "D"  && secondchar == "G") ) return "D";
    if ( (firstchar == "T"  && secondchar == "D") || (firstchar == "D"  && secondchar == "T") ) return "D";

    if ( (firstchar == "A"  && secondchar == "H") || (firstchar == "H"  && secondchar == "A") ) return "H";
    if ( (firstchar == "C"  && secondchar == "H") || (firstchar == "H"  && secondchar == "C") ) return "H";
    if ( (firstchar == "G"  && secondchar == "H") || (firstchar == "H"  && secondchar == "G") ) return "N";
    if ( (firstchar == "T"  && secondchar == "H") || (firstchar == "H"  && secondchar == "T") ) return "H";

    if ( (firstchar == "A"  && secondchar == "V") || (firstchar == "V"  && secondchar == "A") ) return "V";
    if ( (firstchar == "C"  && secondchar == "V") || (firstchar == "V"  && secondchar == "C") ) return "V";
    if ( (firstchar == "G"  && secondchar == "V") || (firstchar == "V"  && secondchar == "G") ) return "V";
    if ( (firstchar == "T"  && secondchar == "V") || (firstchar == "V"  && secondchar == "T") ) return "N";

    if ( (firstchar == "A"  && secondchar == "N") || (firstchar == "N"  && secondchar == "A") ) return "N";
    if ( (firstchar == "C"  && secondchar == "N") || (firstchar == "N"  && secondchar == "C") ) return "N";
    if ( (firstchar == "G"  && secondchar == "N") || (firstchar == "N"  && secondchar == "G") ) return "N";
    if ( (firstchar == "T"  && secondchar == "N") || (firstchar == "N"  && secondchar == "T") ) return "N";

    if ( (firstchar == "R"  && secondchar == "N") || (firstchar == "N"  && secondchar == "R") ) return "N";
    if ( (firstchar == "Y"  && secondchar == "N") || (firstchar == "N"  && secondchar == "Y") ) return "N";
    if ( (firstchar == "K"  && secondchar == "N") || (firstchar == "N"  && secondchar == "K") ) return "N";
    if ( (firstchar == "M"  && secondchar == "N") || (firstchar == "N"  && secondchar == "M") ) return "N";
    if ( (firstchar == "S"  && secondchar == "N") || (firstchar == "N"  && secondchar == "S") ) return "N";
    if ( (firstchar == "W"  && secondchar == "N") || (firstchar == "N"  && secondchar == "W") ) return "N";
    if ( (firstchar == "B"  && secondchar == "N") || (firstchar == "N"  && secondchar == "B") ) return "N";
    if ( (firstchar == "D"  && secondchar == "N") || (firstchar == "N"  && secondchar == "D") ) return "N";
    if ( (firstchar == "H"  && secondchar == "N") || (firstchar == "N"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "N") || (firstchar == "N"  && secondchar == "V") ) return "N";

    if ( (firstchar == "Y"  && secondchar == "R") || (firstchar == "R"  && secondchar == "Y") ) return "N";
    if ( (firstchar == "K"  && secondchar == "R") || (firstchar == "R"  && secondchar == "K") ) return "D";
    if ( (firstchar == "M"  && secondchar == "R") || (firstchar == "R"  && secondchar == "M") ) return "V";
    if ( (firstchar == "S"  && secondchar == "R") || (firstchar == "R"  && secondchar == "S") ) return "V";
    if ( (firstchar == "W"  && secondchar == "R") || (firstchar == "R"  && secondchar == "W") ) return "D";
    if ( (firstchar == "B"  && secondchar == "R") || (firstchar == "R"  && secondchar == "B") ) return "N";
    if ( (firstchar == "D"  && secondchar == "R") || (firstchar == "R"  && secondchar == "D") ) return "D";
    if ( (firstchar == "H"  && secondchar == "R") || (firstchar == "R"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "R") || (firstchar == "R"  && secondchar == "V") ) return "V";

    if ( (firstchar == "K"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "K") ) return "B";
    if ( (firstchar == "M"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "M") ) return "H";
    if ( (firstchar == "S"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "S") ) return "B";
    if ( (firstchar == "W"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "W") ) return "H";
    if ( (firstchar == "B"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "B") ) return "B";
    if ( (firstchar == "D"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "D") ) return "N";
    if ( (firstchar == "H"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "H") ) return "H";
    if ( (firstchar == "V"  && secondchar == "Y") || (firstchar == "Y"  && secondchar == "V") ) return "N";

    if ( (firstchar == "M"  && secondchar == "K") || (firstchar == "K"  && secondchar == "M") ) return "N";
    if ( (firstchar == "S"  && secondchar == "K") || (firstchar == "K"  && secondchar == "S") ) return "B";
    if ( (firstchar == "W"  && secondchar == "K") || (firstchar == "K"  && secondchar == "W") ) return "D";
    if ( (firstchar == "B"  && secondchar == "K") || (firstchar == "K"  && secondchar == "B") ) return "B";
    if ( (firstchar == "D"  && secondchar == "K") || (firstchar == "K"  && secondchar == "D") ) return "D";
    if ( (firstchar == "H"  && secondchar == "K") || (firstchar == "K"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "K") || (firstchar == "K"  && secondchar == "V") ) return "N";

    if ( (firstchar == "S"  && secondchar == "M") || (firstchar == "M"  && secondchar == "S") ) return "V";
    if ( (firstchar == "W"  && secondchar == "M") || (firstchar == "M"  && secondchar == "W") ) return "H";
    if ( (firstchar == "B"  && secondchar == "M") || (firstchar == "M"  && secondchar == "B") ) return "N";
    if ( (firstchar == "D"  && secondchar == "M") || (firstchar == "M"  && secondchar == "D") ) return "N";
    if ( (firstchar == "H"  && secondchar == "M") || (firstchar == "M"  && secondchar == "H") ) return "H";
    if ( (firstchar == "V"  && secondchar == "M") || (firstchar == "M"  && secondchar == "V") ) return "V";

    if ( (firstchar == "W"  && secondchar == "S") || (firstchar == "S"  && secondchar == "W") ) return "N";
    if ( (firstchar == "B"  && secondchar == "S") || (firstchar == "S"  && secondchar == "B") ) return "B";
    if ( (firstchar == "D"  && secondchar == "S") || (firstchar == "S"  && secondchar == "D") ) return "N";
    if ( (firstchar == "H"  && secondchar == "S") || (firstchar == "S"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "S") || (firstchar == "S"  && secondchar == "V") ) return "V";

    if ( (firstchar == "B"  && secondchar == "W") || (firstchar == "W"  && secondchar == "B") ) return "N";
    if ( (firstchar == "D"  && secondchar == "W") || (firstchar == "W"  && secondchar == "D") ) return "D";
    if ( (firstchar == "H"  && secondchar == "W") || (firstchar == "W"  && secondchar == "H") ) return "H";
    if ( (firstchar == "V"  && secondchar == "W") || (firstchar == "W"  && secondchar == "V") ) return "N";

    if ( (firstchar == "D"  && secondchar == "B") || (firstchar == "B"  && secondchar == "D") ) return "N";
    if ( (firstchar == "H"  && secondchar == "B") || (firstchar == "B"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "B") || (firstchar == "B"  && secondchar == "V") ) return "N";

    if ( (firstchar == "H"  && secondchar == "D") || (firstchar == "D"  && secondchar == "H") ) return "N";
    if ( (firstchar == "V"  && secondchar == "D") || (firstchar == "D"  && secondchar == "V") ) return "N";

    if ( (firstchar == "V"  && secondchar == "H") || (firstchar == "H"  && secondchar == "V") ) return "N";

    else
        {
        cout << endl << endl << "comparing " << firstchar << " and " << secondchar << endl;
        cerr << "Unable to compute the consensus of characters";
        exit(1);
        }

    return 0;
}

/**************************************************************************************************************************

Function : get_position_nexus

The function get_position returns the position "i" in a haplotype vector (vector[i]) that has been read from the nexus file

variables : 
	- the name of the string to search for

**************************************************************************************************************************/

int get_position_nexus(string the_string)
{
    int i;
    vector<string> a_vector;
    for (i=0; i < nexus_data.ReturnNTax(); i++)
    {
        a_vector.push_back(nexus_data.ReturnTaxa(i));
    }

    #ifdef DEBUG
    cout << "Vec" << '\t';
    #endif

    int position;
    vector<string>::iterator an_iterator;
    an_iterator = find(a_vector.begin(), a_vector.end(), the_string);
    if (an_iterator == a_vector.end())
    {
        return -1;
    }
    else
    {
        a_vector.erase(an_iterator, a_vector.end());
        position = a_vector.size();

        #ifdef DEBUG
        cout << "Pos" << '\t';
        #endif

        return position;
    }
}

/**************************************************************************************************************************
 
 Function : AttributeSequencesToOrganisms
 
 Read the taxa names and attributes them to an organism from the organism file. Organisms need to start the name of the
 allele and is separated from the rest of the allele name by a '_'. ex: Rnitida761_A, Rnitida761_B, etc.
 
 variables :
 - dataset number being treated
 
 **************************************************************************************************************************/

void AttributeSequencesToOrganisms(void)
{
    int i,x;
    string temp,an_organism;
    org_data.InitializeAllelesInOrganisms();
    for (i=0; i < nexus_data.ReturnNTax(); i++) {
        an_organism="";
        temp = nexus_data.ReturnTaxa(i);
        an_organism+=temp[0];
        an_organism+=temp[1];
        if (org_data.AddAllele(an_organism,nexus_data.ReturnTaxa(i)) == 1) {
            cout << " Organism " << an_organism << " is not present in organisms list." << endl;
            exit(1);
        }
        //cout << "added " << nexus_data.ReturnTaxa(i) << " in " << an_organism << endl;
    }
}


//------------------------------------------
// get_PBC_distance
//------------------------------------------

double get_PBC_distance(int org1, int org2)
{
	int x,j;
    //If we compare the individual with itself
    if (org1 == org2) return (0);
    //If one fo the organism is absent from the matrix
    if (org_data.ReturnNbAllelesforOrganism(org1)==0 || org_data.ReturnNbAllelesforOrganism(org2)==0) return(-999);
    
    double distance=0.;
    double tempdist1,tempdist2;
    for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
        tempdist2=9999999.0;
        for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
        distance+=tempdist2;
    }
    for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
        tempdist2=9999999.0;
        for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
        distance+=tempdist2;
    }    
    distance = distance/(org_data.ReturnNbAllelesforOrganism(org1)+org_data.ReturnNbAllelesforOrganism(org2));
    return distance;
}


//------------------------------------------
// get_MIN_distance
//------------------------------------------

double get_MIN_distance(int org1, int org2)
{
	int x,j;
    //If we compare the individual with itself
    if (org1 == org2) return (0);
    //If one fo the organism is absent from the matrix
    if (org_data.ReturnNbAllelesforOrganism(org1)==0 || org_data.ReturnNbAllelesforOrganism(org2)==0) return(-999);
    
    double distance=0.;
    double tempdist1,tempdist2;
    
    tempdist2=9999999.0;
    for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
        for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
    }
    distance=tempdist2;
    return(distance);
}



/************************************************************************************

Function : read_organisms

Description: The function reads the organisms from a file and places them in a vector of organisms

Arguments:   - the name of a file containning the organisms 

************************************************************************************/


void read_organisms(string name)    //When using the batch mode
{
    string stemp;
    char *BigBuffer;

    ifstream inorganisms(name.c_str());

    while (!inorganisms)   //Check if the input file can be found
        {
        cout << endl << "Error with input organisms file: " << name << endl;
        cout << "Make sure the file exist and that it contains a valid organisms block" << endl;
        cerr << "exiting";
        exit(1);
        }

    BigBuffer = new char[MAX_BUFFER_SIZE];
    if (!BigBuffer)
    {
        cerr << "Could not allocate file buffer";
        exit(1);
    }

    for (i=0; !(inorganisms.eof()); i++)
        {
        if (i >= MAX_BUFFER_SIZE-1)  //have to save room for terminating null
            {
            cerr << "Nexus file exceeds 500k maximum";
            exit(1);
            }
        inorganisms.get(BigBuffer[i]);
        }        

    inorganisms.close();
    BigBuffer[i]='\0';
    bufPtr=BigBuffer;

    if ( bufPtr != NULL )
        {
        while ((aTokenPtr=nextToken()) != "") //This has to be checked
            {
            if (aTokenPtr == "BEGIN")
                {
                stemp = nextToken();   /* get the block name and store in 'stemp'*/
                if (stemp == "")
                    {
                        cerr << "Error reading block name";
                        exit(1);
                    }
                if ((aTokenPtr=nextToken()) == ";") /* pop the terminating semicolon */
                    {
                    if (stemp == "ORGANISMS")
                        {
                        do            /* need to put in error checking in case no DIMENSIONS statement */
                            {
                            aTokenPtr=nextToken();
                            if (aTokenPtr == "DIMENSIONS")  
                                doOrgDimensions();
                            if (aTokenPtr == "ORGLABELS")
                                doOrgLabels();
                            }
                        while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
                        aTokenPtr=nextToken();
                        if (aTokenPtr != ";") {
                            cerr << "Error reading block name";
                            exit(1);
                        }
                        }
                    else  /* token is not a recognized block */
                        {  
                        doUnrecognizedBlock();
                        }
                    }
                stemp.clear();				
                }
            }
        }

    if (info_data.IsVerbose())
        {
        cout << endl << "File containning the organisms: " << name << endl;
        cout << "Number of organisms read: " << org_data.ReturnNOrg() << endl;
        cout << "Content of file: " << endl;
        for (i=0; i < org_data.ReturnNOrg(); i++)
            {
            cout << org_data.ReturnOrganism(i) << '\t';
            }
        cout << endl << endl;
        }

    info_data.IsOrganismsFile();

    delete [] BigBuffer;
}


/***********************************

Function : doOrgDimensions

***********************************/

void doOrgDimensions(void)
{
    while ( (aTokenPtr=nextToken()) != ";")
        {
        if (parse_assignment2("NORG"))
            {
            org_data.GetNOrg( atoi(LocalToken.c_str()) );
            }
        }
    return;
}


/***********************************

Function : doOrgLabels

***********************************/


void doOrgLabels(void)
{
    int size;
    while ( (aTokenPtr=nextToken()) != ";")
        {
        org_data.ReadOrganism(aTokenPtr);
        }
    size = org_data.NumberOrgLabels();       //Number of organisms that were read
    if ( size < org_data.ReturnNOrg()) {
        cerr << "Too few organisms labels";
        exit(1);
    }
    else if ( size > org_data.ReturnNOrg()) {
        cerr << "Too many organisms labels";
        exit(1);
    }
    org_data.OrganismsRead();                //Flag to indicate the organisms were read
    return;
}

/******************************************************

Function: To_Uppercase()

Takes a string as input and outputs an uppercase string

******************************************************/
/*
string To_Uppercase(string input_string)
{
    int i;
    char character;
    string output_string, temp_string;
    temp_string = input_string;
    int size_of_string = temp_string.size();

    for (i=0; i < size_of_string; i++)
    {
        character = temp_string[i];
        output_string += toupper(character);
    }
    return output_string;
}
*/
/*******************************

Function : doUnrecognizedBlock

*******************************/

void doUnrecognizedBlock(void)
{
	do 				
		{
		aTokenPtr=nextToken();
		}  while ( (aTokenPtr != "END")  &&
						(aTokenPtr != "ENDBLOCK" ) );
	aTokenPtr=nextToken();  /* Pop the terminating semicolon */
	if (aTokenPtr != ";") {
        cerr << "Block not terminated by a semicolon";
        exit(1);
    }
	return;
}





/****************************************************************/
int parse_assignment2(string target)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

{
	if (aTokenPtr == target)
			{
			aTokenPtr=nextToken();
			if (aTokenPtr != "=") {
                cerr << "Bad assignment statement= " << aTokenPtr << endl;
                exit(1);
				}
			aTokenPtr=nextToken();
			LocalToken = aTokenPtr;
			return 1;
			}
	else return 0;
}


/***************************************************

Function: nextToken

***************************************************/

/*
	Function originally written by Michael J. Sanderson
	original name: nextToken2.c
	included in the program r8s
*/


/*	
	Gets the next token from input stream 'fpointer', and copies it onto the global
	buffer pointed to by 'aTokenPtr'.  If there is NO next token, we copy a null
	string onto that buffer.  That's a signal for the main caller routine...

	If the global variable gNewLine=1 then the newline characters, '\n' and '\r'
	ARE returned as individual tokens,  when encountered.  The normal state is
	gNewLine=0,  which treats these as white space delimiters too.  The only time
	NEXUS file needs to think about newlines is when reading interleaved matrices!

*/

string nextToken(void)
	{
	string TempPtr;		// a pointer to manipulate aTokenPtr
	char c;
	
	if  ((c=*bufPtr++) == '\0') NULL_RETURN;
	
	while (( isNEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */
		{
		while ( isNEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
			{
			c=*bufPtr++;
			if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
			}
			    
		if (c=='[')		/* skip the comment and land on next 'c' after comment */
			{
			while (c !=']')
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				}
			c=*bufPtr++;	/* get next char after ']' */
			if (c=='\0') NULL_RETURN;
			}
		}	    
    
	if (c=='\'')		/* deal with single-quoted tokens */
		{
		TempPtr.push_back(toupper(c));
		while (  (c=*bufPtr++) != '\'')
			{
			if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			TempPtr.push_back(toupper(c));	/* this is a valid character in the word, add to token */
			}
		TempPtr.push_back(toupper(c));	/* add the terminating quote too */
		return(TempPtr);
		}	/* return everything between single quotes, including the quotes, as a token*/		

	TempPtr.push_back(toupper(c));		/* char is either punctuation or part of word, so add it to token */
	    
	if (!isNEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
								    or Token size exceeded */
	    {
	    for (;;)
		    {
		    c=*bufPtr++;
		    if (  isNEXUSpunct(c) || isNEXUSwhiteSpace(c) ||  (c == '\0') )
			    {
			    --bufPtr; /* word is terminated by some c that is not part of word;
							     push c back into stream and deal with it on
							    next call to this function; meantime, break out, 
							    and return this token*/
			    break;
			    }
		    TempPtr.push_back(toupper(c));	/* this is a valid character in the word, add to token */
		    }
	    }
	return(TempPtr);
	}

