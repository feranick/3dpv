// SGA VERSION (Without OPENGL)
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <sstream>
//#include <direct.h> // for getcwd 
//#include <iomanip.h> 
#include <fstream>
#include <assert.h>
#include <list>     // don't think i'm using list anymore
#include <map>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>

#include <unistd.h>

// Genetric algorithm (stuff that WAS included in "userDefinables.cpp")
// make sure to put stuff in the right directories, or indicate correct path
#include "ga.hpp"
#include "random.hpp"
#include "globalSetup.hpp"

// other function libraries
#include "spa.h"
#include "geometry.h"

// my headers
#include "Cell_3.h" // contains, Cell_3 class, 
#include "main.h"

using namespace std;

// ascii codes for certain keys
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

int Ggen_show=0;  // which generation's individual to show in opengl window
double Gtime_show=0; // at what time interval to place the "sun" to view the structure in opengl
			// this is based on the same day-month-year "trajectory" that was used in calculations

// for making glutinit function happy
/*int intforglut=0;
char **charforglut;*/
//test-comment
//##################################################################
// CONSTANTS
//##################################################################
#define PI           3.1415926535897932384626433832795028841971
enum Shapes {MODEL, CUBOID, TETRAHEDRON, TRIFLAT, RECTFLAT};    // "MODEL" is loaded from file
enum DataSets {USER, BERKELEY, PALMDALE, MIT}; // default data set examples
enum CrossType {UNIFORM, ONEPOINT, TWOPOINT};
const int Rays_numberper_iteration = 1000;

// INPUT_examples is an array for default sets of simulation parameters. The first index is the set, second is the parameter.
 // [2][6] (index) (longitude, latitude, avgtemp, avgpressure, elevation, time from GMT)
double INPUT_examples[][6] = {
		{        0,       0,       0,       0,  0,  0},  // [0][] = don't use default examples
		{-122.2581,-37.8709, 14.22, 970.5, 86, -8},   // [1][] = berkeley values (Birge Hall)
		{-118.16259384155273,34.58417227544041, 18.2871, 71.6794, 2600, -7},
		{-71.09248, 42.35994, 10.8889, 1013.25, 22, -5} //cambridge MA values
		};
DataSets input_i = MIT;      // which input set to use (USER(0) for not any defaults)
double globalInputLatitude=30;	  // this will be set to the appropriate value automatically when simulation begins

////// FROM GA userDefinables ////////////////
// these are constructs used by the GA source codes, but for functions used somewhere in main.cpp
#define BLANK_STR " \t\n\r"
GlobalSetup *globalSetup;	// (see globalSetup.cpp globalSetup.hpp)
Random myRandom;		// a class of functions for random numbers (see random.cpp)
/////////////////////////////////////////////

//###################################################################
// INPUT VARIABLES
//###################################################################
// constraints not dependent on genetic algorithm
// these are typically set by reading in the main input file in "sim" operation mode
// see 3DPV documentation for format of input file and variables used
int numStdInputs = 18;
int numGAInputs = 9;
double numSigFigs = pow(10.0,6.0);
double xBound = 68;  // ground-plane footprint in x (east-west). length is in meters
double yBound = 68.0;// same for y-direction (north-south)
double Footprint_xmin = 0.0; 
double Footprint_ymin = 0.0;
double zBound = 68;     // ie, maximum height permitted for structure
double Structure_zmin = 0;   // might as well keep as zero (ground-height), but could change for interesting floating structure
int startYear = 2008; 	// four-digit year
int startMonth = 6;  	// two-digit month
int startDay = 21;        // two-digit day
int endYear = 2010;
int endMonth = 10;
int endDay = 31;
bool Delete_sundata=true; // whether to append or remake (for comparing sun trajectories more easily)
bool Done_sundata=false;  // since each GA only over 1 same day, only plot this data ONCE to avoid large useless filesize of repeated data
bool Write_tsampledata=false; // if true, records sample points from 1st individual 1st time then switches to false
bool Write_vertdata=true;

// for shadows
enum RAYMETHOD {PERAREA, PERCELL};  // the PERCELL method is used now. changing this to PERAREA would require many code changes below....
RAYMETHOD ray_method=PERCELL;
int rays_per_unit_area=2;  // this is not used unless you use defunct PERAREA method
int raysPerCell=10;       // # rays to check per each cell. note, again that input value for this (maininput file) will be squared

// when loading initial population, these models will become the individuals in the pop.
// if i do this, input file names are here no longer const char but just char)
char Input_modelALL[100][100];
string modelSpecDir;
int mode = 1;
int numInputVars = 28;
///////// GA, etc variables that ARE set by input file or command line, so value here shouldn't matter/////////////
/*const */int numCells = 48;	// this is just to initialize number of solar cell panels. Variable is set in input file
int GA_population=12;       // size of population in genetic algorithm, ie # of individuals. Also set in input file
/*const */int GA_maxgenerations=4;   // maximum number of generations (eventually may be 100+). Also set in input file
int GA_numDistinctModels=2; 		// see input file 3DPV documentation
double GA_replace_ratio = 0.5;		// set by input file, also. see input file 3DPV documentation
int GA_select_tourneysize = 			6;
double GA_mutation_probability = 		0.01;
double GA_cross_probability = 			0.8;
int GA_willload_population=1;      // true, else if 0 then random initialization
int doSingleReflections = 1;    // if true, then any reflected parts of rays get a "single chance" to hit another cell
				// assumption is that on any FURTHER reflection it would be insignificant power gain.
////////////////////////////////////////////////////////////////////////

char Main_Dir[255/*200*/];	// main directory. see "sim" part of main() function 
char Clone_Dir[255/*200*/];	// this is made a copy of unique directory and may be changed/copied repeatedly when writing some files
char Unique_Dir[255/*200*/];	// this is set to the unique output directory in the main() function
char Input_CMD[255/*200*/]; // formerly "*Input_CMD" but this caused a segfault I think

/// Here are all the base-names for (most of?) the output files
// these are used when writing and reading files.
const char *Input_GA = 					"inputGAfile.txt";
const char *OutMain_name=				"MAINoutput.txt";
const char *OutModel_prefix=				"Model_gen";
const char *OuthAngle_prefix=				"HorizAngles_grid";
const char *OutpAngle_prefix=				"PolarAngles_grid";
const char *OuthAngleRefl_prefix=			"HorizAnglesREFL_grid";
const char *OutpAngleRefl_prefix=			"PolarAnglesREFL_grid";
const char *OutallAngle_prefix=				"AllAngles_grid";
const char *OutSpherical_prefix=			"AllSphere_grid";
const char *OutMidDist_prefix=				"ToCenterDistances_gen";
const char *OutCountMidDists_prefix=			"CountMidDists_gen";
const char *OutPairCorrelation_prefix=			"PairCorrelation_gen";
const char *OutBest_Indivs=				"outputAllBestGen";
const char *OutBest_Obj=				"outputBestFit";
const char *OutAppBest_Obj=				"outAppBestFit";  // 4/29/09. like previous but appended each generation 
const char *OutAll_Indivs=				"outputAllIndiv";
const char *OutSun_Path=				"outputSunPath";
const char *OutPowerPerCell=				"outputCellPow";	// OBSOLETE FILE OUTPUT
const char *OutPowerMesh=				"outPowerMesh_";	// OBSOLETE FILE OUTPUT
const char *OutPowerPerSubcell=				"outSubcellEnergy_best";
const char *OutGridCoordinates=				"outSubcellCoords_best";
const char *OutCoordANDEnergy=				"outSubcellAllCoordsEnergy"; // combination of all coordinates data with all
											// iterations' energies.
const char *OutTriGridCoords=				"outSubcellTriCoords_best"; // new on 02/11/2010 for opengl visuals

const char *OutBest_Power=				"outputBestPower";

///////// a number of other GA variables that are not set in the input file///////////////////////
// so if you need to change them, this is the place to do it, though most have been optimized
// by trial and error for the 3DPV stuff
int GA_numberof_objectives = 1;	// SGA with one objective ===> just POWER, no POWER/MATERIAL. do not change this
int GA_mass_extinction_rate = 			50;      //# of generations before wiping out all but one, periodic
double GA_extinction_ratio = 			0.99999; // should always be getting rid of all but one
int GA_use_extinction = 			0;       // 0 =false, 1= true do use extinction with above parameters
						// NOTE, there is no extinction GA operation at this time, so these variables are unused
int GA_cross_type = 				TWOPOINT;//ONEPOINT;//UNIFORM; // swap two solutions genes in segment of chromosome
double GA_genewise_swap_probability = 	0.01; 
const char *GA_initial_population_file = "initialPopulation.txt";
int GA_willsave_solutions = 			1;       // true
const char *Output_GA_sol = 			"solutions1.txt";// this is NOT USED ANYMORE 
///////////////***************************/////////////////////////////////////////

//////// Don't need to/shouldn't change these, they are initialized/////////////
int GA_gen_current = 0; 	// just the number ID of current generation. "genID" in the GA object
bool getSubCellData = false;	// used as true when "grid" commnd line option is used, knows to collect subcell data
int numCoordinates=0;
double *modelPointer; // global pointer to first position in array of model coordinates
// shouldn't be needing the allDecVarPointer for loading populations after all
//double *allDecVarPointer; // global pointer to first position in array allModelDecisionVars 
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// various inputs for the sun-position-algorithm (SPA)
double longitude;  // earth-location longitude
double latitude;   // earth-location latitude (note the SPA has some odd +/- sign quirks here, I think, so be careful)
double avgTemp;	// not too important, default set in SPA function below	
double avgPressure; // not too important, default set in SPA function below	
double elevation; // not too important to change?, default set in SPA function below
double timeZone; // default set in SPA function below	
spa_data spa; /*create instance of data structure for sun position calculations*/

// dependent on genetic algorithm (these are just default values for now)
int Cells_baseshape = TRIFLAT;//MODEL;
int Cells_radius = 2; // NOT USED ANYMORE 
double gRatio = (1+sqrt(5))/2; // golden ratio for some of the polyhedra coordinates
//int gotcenter=0; // get general center of the structure true or false
//double midcommon[3] = {0, 0, 0};
bool globalCellsOneSided = false;

double delta_time = 0.2;  // time between sun-repositioning (in fractional hours)

// contain pointers to the cell class objects
map<double, double> DayPowerValues; // key = fractional time, value = power at this time
double totalSurfaceArea=0;
// NOTE: area restriction is NOT USED ANYMORE, requires re-modification if to be used again
int allow_area_restrict=0; // 1 => don't restrict area of cells. 1 = do restrict area of the cells
double max_percell_area=14; // maximum allowed surface area for a triangle (one-side)
double min_percell_area=12;  // minimum allowed surface area for a triangle (one-side)
double max_perdim_dist=5;
double min_perdim_dist=3;

// Power calculation variables/constants
const double solar_constant = 1000;  // insolation  // this "solar_constant" refers to insolation at perp surface to rays (W/m^2)
double cellEfficiency = 0.06; // assume plastic cells are 6% efficient, just a constant factor here
const double Ws_to_kWhr = 2.777778*pow(10.0,-7.0); // ie 1kWhr/(3600s*1000W), conversion factor
bool globalUseKWH=false;	// units, note if this is false there are still some things output in KWH. probably shouldn't change

// Material properties
const double index_ref_air = 1.000293;
const double index_ref_plastic = 100000.0;    // values are somewhere between 1.46 to 1.55
double relative_permeability = 1;/*(1+)/(1+)*/ // u1/u2 = u0(1+Xm1)/u0(1+Xm2)

////////////// SOME INITIALIZATION OF GLOBAL STORAGE VARIABLES FOR OUTPUTS AND STUFF///////////////////////////////
double globalDataParams[1000000]; // changed 3/9, then again 3/17 it now holds all coordinates for ONLY 1 generation at a time
//int	globalParamsCounter=0; // keeps track of next open index for globalDataParams, but unused as of 3/17
//double globalDataFitness[1000010];
double globalDataObjective0[2000000];
double globalDataObjective1[100001];
int globalObjCounter=0; 

double globalParamsAvgBest[160000000]; // of all global data params choose best individual from each generation and store it here
double globalObj0AvgBest[1000010]; //note, can knock a zero or two off this ---3/17
double globalObj1AvgBest[1000010]; // note, can knock a zero or two off this ---3/17, havent done it though

double globalDataTEMP[10000000];

double qwikStoreCoords[16][10000]; 	// the "qwik" stuff is for an efficiency thing I added that gets around a kind of
double qwikStoreEnergies[16];		// annoyance in the GA. If it so happens that there are two structures produced
double qwikStoreReflEnergies[16];	// which are "genetically-identical" and within a few generations of each other
int qwikIndex=0;			// then instead of calculating the daily-power again and wasting time, the known
int qwikIndMax=16;			// value will just be attributed to it
int qwikCoordMax=10000;

int globalDistrPolarAngle[360];		// used in "grid" mode. stores distribution of cells with particular polar angle orientations
int globalDistrHorizAngle[360];		// same for horizontal angle (ie, theta in spherical coordinates with z axis skyward)

int AngleBins_numberof=60;//15; 
int HorizBins_numberof=60;//15; 
double globalCountDistrPAangleBins[60];
double globalCountDistrHAangleBins[60]; // changed from 'int' to 'double' on 6/04/09 for weighting values
double global2dDistrBothAngleBins[60][60];
double globalCountDistrPAangleReflBins[60];
double globalCountDistrHAangleReflBins[60];
//double globalDistrMidDistance[10000000];
int globalCountDistrMidDistance[100000];
int globalMaxSoFar=0;

// added 3/16/09 for power SGA in particular, save percell power of all of 1 gen, then reuse array after save each gen to file
double globalDataCellPower_OneGuy[1000]; // assumes no more than 1000 cells...more than enough space
double globalDataCellPower_Gen[10000];

// 4/6/09, shouldn't be using globalDataCellPower arrays above anymore,,,,,instead:
double globalDataSubCellPower_OneGuy[5000000]; // requires 160cell*16subcell*100dayDivisions ~256,000 slots
double globalDataSubCellPower_Reflect[5000000]; //new 6/25/09. may be useful for analysis of something
double globalDataSubCellPower_Gen[8000000]; // requires slots for all data above times 8 guys in a population
double globalDataSubCellPower_Sum[10000]; 
double globalDataSubCellPower_ReflSum[10000];
double globalBestPower_thisGen=0; // keep track during evaluations of each generation, which is the best power
int globalBestGuy_thisGen=0;  // keep track of the individual who has the above-stored best power in generation so far (0, 1, 2...7)
int globalNumDayDivisions=0;

///////////// added 8/19/09 to separate triangles....//////
// this is a somewhat cumbersome quick-fix. but essentially, all the
// cells are double-sided and these variables contain the information for the "other" side
// if the above versions of these arrays contain the info for the "First" side
double SECONDglobalDataSubCellPower_OneGuy[5000000]; 
double SECONDglobalDataSubCellPower_Reflect[5000000]; 
double SECONDglobalDataSubCellPower_Gen[8000000]; 
double SECONDglobalDataSubCellPower_Sum[10000]; 
double SECONDglobalDataSubCellPower_ReflSum[10000];
//////////////////////////////////////////////////////////

//added on 7/29/09 for subcell stuff 
// "subcells" are the triangles in the plane of a given triangle solar cell
// which each get one ray trace. A good picture to have in mind I think is a "Triforce" decomposition
// where each sub-triangle in the whole has one ray trace to the sun, and thus may be shadowed or not.
// thus the number of sub cells is related to the input number of ray traces (squared)
// the "gridpoints" are the midpoints of each of these subcells, which is the position that the raytrace actually shoots from
double globalDataSubCoords[6000000]; // say 200 cells * 10,000 gdpt * 3 coord's per gdpt = 6,000,000 slots
int globalCountForRecordCoord=0; // don't change, initialized to count from 0.

// more global variables, don't change counters to start at 0, or "havedone" variables for intializing processes just once
int globalIndividualCounter=0; // for the above two arrays' data, increments in globalEval but resets to zero once GA_population value
bool globalHaveDonePreEval=0;

const int maxNumGrids = 90000;

double globalTriGrid[2*maxNumGrids]; // max number of grid points is 64000/2 = 32000...
double globalTriMids[2*maxNumGrids];	// this is really the ultimate set of triangle 2d gridpoints!!! midpoints of subdivided triangles
double globalTriVerts[6*maxNumGrids]; // NEWLY ADDED on 02/11/2010 to help with new OpenGL visualization video
				// the tri verts are actually the subcell vertices whose midpoints are the raytrace "gridpoints"
int use_trigrid=1; // =0 for no, do randomized sample. =1 for yes, use a uniform grid of triangle points
		// DON'T set this to 0. the randomized stuff is OLD and probably doesn't even work anymore. inaccurate, if anything...

double globalMostRecentEnergyReturned=0; // new in august (2009), for last globalEval quick power return value

// new for the calculation of roof stuff
// the roof operation modes of "solar3d" are pretty under-developed, and need work if they
// are to be expanded upon. but basically, they test the power over a day for a bunch of identical
// structures positioned on some flat "roof" in some close-configuration. thus, the calculation may take very long for a lot of copies. 
bool globalDoRoof=false;
double globalRoofx=0;
double globalRoofy=0;

//////////////////////////////////////////
// StoreBestIndividual_Gen() is alternative to below function
// added 3/17/09 because it should require a lot less memory, write to disk every generation or so instead of at end
void StoreBestIndividual_Gen(int ggen)
{
	//double Sum0=0;  //NSGA
	//double Sum1=0;  //NSGA
	//double stake[GA_population]; //NSGA
	//double max_stake=0; // NSGA
	//printf("in storebest ggen is %i\n", ggen);
	double max_power=0;//SGA
	int best_dude=-1;
	//int cskip = GA_population*9*numCells; // constant integer multiplier, unused as of 3/17 too
	int bskip = 9*numCells;

	//int gg=0; // generations
	int ii=0; // individuals
	int cc=0; // cells
	int vv=0; // variables
	
	// For single objective, its easy to get the best individual....
	best_dude=-1;
	max_power=0;

	for (ii=0; ii < GA_population; ii++)
	{
		if (globalDataObjective0[ggen*GA_population + ii] > max_power)
		{
			max_power = globalDataObjective0[ggen*GA_population + ii];
			best_dude = ii;
		}
	}
	for (cc=0; cc<bskip; cc++)
	{
		for (vv=0; vv < 9; vv++) {
			// commented out 3/17/09 to use much smaller globalDataParams array
			//globalParamsAvgBest[gg*bskip + cc*9 + vv] = globalDataParams[gg*cskip + best_dude*bskip + cc*9 + vv];
			globalParamsAvgBest[ggen*bskip + cc*9 + vv] = globalDataParams[best_dude*bskip + cc*9 + vv];
			globalObj0AvgBest[ggen] = globalDataObjective0[ggen*GA_population + best_dude];			
		}
	}

//JIN: for the sake of uniformity, let's output all amounts of energy at the end
//	// 4/30/09// each generation this appends the best fitness. (previosuly only had file printed at end of all GA)
//	strcpy(Clone_Dir, Unique_Dir);
//	char *tmp_strb = strcat(Clone_Dir, OutAppBest_Obj);
//	FILE *outFileObj = fopen (tmp_strb, "a");
//	fprintf(outFileObj, "%i %f\n", ggen, globalObj0AvgBest[ggen]);
//	fclose(outFileObj);
	/////////////////////////////////////////
}


// CALCULATE THE PAIR CORRELATION FUNCTION FOR AN INDIVIDUAL OF CELLS (pass in generation and look at bestAvg individual)
// NOTE THAT it has probably been since early 2009 that I did anything with this paircorrelation
// function, and I never put it to good use. May be an interesting idea if reworked. I forget how well it worked.
void GetPairCorrelationFunction(int thisgen)
{
	double coords_me[9];	// current guy we are looking at
	//double coords_other[9];  // shorthand storage of other guy to compare distance to
	double centers[numCells][2]; // stores x, y, z of the centers of the cells
	double distances[numCells][numCells]; // store distance between cell i and cell j [i][j]
	double temp_d=0;
	int temp_bin=0;
	
	int cc=0;
	int vv=0;
	double dr = 2.0; // this interval will be constant
	double norm_add = 1/(((double)numCells)-1); // partially normalized increment (divide by number of neighboring cells considered)
	
	// number of distribution bins is the largest possible distance (box corner-to-corner) divided by the interval dr
	int num_dist_bins = (int)(sqrt( pow(xBound-Footprint_xmin, 2)+
					pow(yBound-Footprint_ymin, 2)+
					pow(zBound-Structure_zmin, 2))/dr); printf("number of distribution bins is %i\n", num_dist_bins);
	
	double pair_function[numCells][num_dist_bins];	// the main array to hold counts of how many guys in each bin (spherical shell) from each cell cc
	double total_pair_function[num_dist_bins]; // averaged? or just summed? pair_function[][] over all cells, cc, used as references !!
	
	double r = 0.01; // the radius from current cell being considered
	
	double cell_number_density = (double)numCells / ((xBound-Footprint_xmin)*(yBound-Footprint_ymin)*(zBound-Structure_zmin));
	printf( "Adding factor: %f and cell number desnity: %f\n", norm_add, cell_number_density);

	// first get the centers of all the cell triangles
	for (cc=0; cc < numCells; cc++)
	{
		for (vv=0; vv < 9; vv++)
		{
			coords_me[vv] = globalParamsAvgBest[thisgen*numCells*9 + cc*9 + vv];
		}
		centers[cc][0] = (coords_me[0] + coords_me[3] + coords_me[6])/3;
		centers[cc][1] = (coords_me[1] + coords_me[4] + coords_me[7])/3;
		centers[cc][2] = (coords_me[2] + coords_me[5] + coords_me[8])/3;
	}
	printf ("Got past centers\n");
	// then calculate the distances between all the cells
	int iii;
	int jjj;
	for (iii=0; iii < numCells; iii++)
	{
		for (jjj=0; jjj < numCells; jjj++) // only look at higher indices, to not double-calculate
		{
			if (jjj>iii)
			{
				distances[iii][jjj] = 	sqrt(pow(centers[iii][0]-centers[jjj][0],2)+
							pow(centers[iii][1]-centers[jjj][1],2)+
							pow(centers[iii][2]-centers[jjj][2],2));
				//printf("jjj> iii distance: %f\n", distances[iii][jjj]);
			}
			else
			{
				if (jjj<iii)
				{
					distances[iii][jjj] = distances[jjj][iii]; 
					//printf("jjj<iii distance: %f\n", distances[jjj][iii]);
				} // eg, same distance 0->3 as 3->0
			}
			if (jjj==iii) {distances[iii][jjj]=0;} // same guy
		}	
	}
	printf ("Got past getting all distances\n");
	// second, select cell #i and a value dr and see how many others lie in r to r + dr
	for (cc=0; cc < numCells; cc++)
	{
		for (int c2=0; c2 < numCells; c2++)
		{
			if (c2!=cc)	// don't count itself in the 0 to 0+dr bin
			{
				temp_d = distances[cc][c2];
				temp_bin = (int)(temp_d/dr);
				// third, divide each count by N, the total number of considered neighbors (numCells-1) (done above already)
				pair_function[cc][temp_bin] += norm_add; // add 1 to the bin for this distance
				// dividing each +1 by N # of cells is same as summing then dividing by N # of cells, of course			
			}
		}

			
		for (int abin=0; abin < num_dist_bins; abin++)
		{
			// fourth, divide each count by 4*pi*r^2*dr to make up for that at larger r there is more space to check in spherical shell....
			// NOTE: this may turn out wrong unless I somehow take into account that on the boundary edges there are obviously no cells outside so
			//	checking not a sphere
			// but some other kind of volume
			//FOR NOW, 3/10/09 I will just divide by the full volume of the sphere...but the data probably won't look right
			double volume_factor = ( 4*PI*(r+(double)abin*dr)*(r+(double)abin*dr)*dr );
			//printf("volume factor: %f\n", volume_factor);
			pair_function[cc][abin] /= volume_factor; // this is just 4*pi*r^2*dr
			
			// fifth, divide by the number density (here this should just be my total number of cells divided by boundary volume, a box)
			pair_function[cc][abin] /= cell_number_density; // same for all => no effect on shape of curve? just normalize numbers
		}
	}

	// sum them all
	for (int abin=0; abin < num_dist_bins; abin++)	
	{
		for (int cc = 0; cc < numCells; cc++)
		{
			total_pair_function[abin] += pair_function[cc][abin];
		}
	}

	// write this PAIR CORRELATION STUFF TO FILE
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat = strcat(Clone_Dir, OutPairCorrelation_prefix);
	//strcat(output_cat, itoa(whichgen));
	char istring[200];
	sprintf(istring,"%i",thisgen);	
	strcat(output_cat, istring);
	//strcat(output_cat, "_cells.txt");
	//strcat(output_cat, itoa(numCells);

	FILE *OUTfile = fopen(output_cat, "w");
	
	for (int pcf=0; pcf<num_dist_bins; pcf++)
	{	
		fprintf(OUTfile, "%f\n", pair_function[0][pcf]/*total_pair_function[pcf]*/);
	}

	fclose(OUTfile);
}

//#####################################################################
// MATH FUNCTIONS
//#####################################################################
double radtodeg(double radians)
{
    return (180.0/PI)*radians;
}
double degtorad(double degrees)
{
    return (PI/180.0)*degrees;
}

//########################################################################################
// ANALYSIS & DATA GATHERING FOR DISTRIBUTIONS
//########################################################################################

void WriteCountDistanceDistrFile(int whichgen)
{
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat = strcat(Clone_Dir, OutCountMidDists_prefix);
	char istring[200];
	sprintf(istring,"%i",whichgen);	
	strcat(output_cat, istring);

	FILE *OUTfile = fopen(output_cat, "w");
	
	for (int cff=0; cff<10000/*globalMaxSoFar*/; cff++)
	{	
		if (globalCountDistrMidDistance[cff]!=0) {
			fprintf(OUTfile, "%i\t%i\n", cff, globalCountDistrMidDistance[cff]);
		}
	}

	fclose(OUTfile);
}

//#############################################################################
// GENERAL OPENGL FUNCTIONS
//#############################################################################

// WriteTextModelFile writes an output (.txt) with model coordinates for the
// most fit/best-performing structure of a particular generation number (whichgen)
void WriteTextModelFile(int whichgen)
{
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat = strcat(Clone_Dir, OutModel_prefix);
	char istring[200];
	sprintf(istring,"%i",whichgen+0);	
	strcat(output_cat, istring);
	strcat(output_cat, "_cells.txt");

	FILE *OUTfile = fopen(output_cat, "w");
	
	for (int cff=0; cff<numCells; cff++)
	{	
		for (int cdd=0; cdd<9; cdd++)
		{
			fprintf(OUTfile, "%f", globalParamsAvgBest[whichgen*numCells*9+cff*9+cdd]);
			if (cdd<8) {fprintf(OUTfile, ",");}
			if (cdd==8) {fprintf(OUTfile, "\n");}
		}
	}

	fclose(OUTfile);
}

// ReadBestObjFile. YOU CAN "find" on this function, but apparently I don't use it anymore
// in this code. It loaded some best-power values from a previous run...

//JIN: not used currently
void ReadBestObjFile ()  // takes pointer to first element in new coordinate array
{	
	//double outputs_stuff[14];
	
	strcpy(Clone_Dir, Unique_Dir);
	FILE * bestFile = fopen ( strcat(Clone_Dir,  OutBest_Obj), "rb");


	long lSize;
	char * buffer;
	size_t result;

	//if (bestFile==NULL) {fputs ("pre-view loading : File error, main output file\n",stderr); exit (1);}
	
	// obtain file size:
	fseek (bestFile , 0 , SEEK_END);
	lSize = ftell (bestFile);
	rewind (bestFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,bestFile);
//	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}   // check gives integer comparison warning

	// the whole file is now loaded in the memory buffer.

	// terminate
	fclose (bestFile);
	
	int tstore0=0;
	//int tstore1=0;  // NSGA only
	int cycle=0;
	
	// SECOND PASS: store data as array of doubles
	char *toks;
	toks = strtok(buffer, " ,\n,\t");            // break file contents into tokens delimited by commas and whitespaces
	while (toks!=NULL)	
	{
		if (cycle==0)
		{
			// do nothing
			cycle=1;
		}
		else {
			if (cycle==1)
			{
				//cycle=2; //NSGA
				cycle=0; // SGA
				globalObj0AvgBest[tstore0]=atof(toks);
				tstore0++;
			}
			/*else { //NSGA
				if (cycle==2)
				{
					cycle = 0;
					globalObj1AvgBest[tstore1]=atof(toks);
					tstore1++;
				}
			}*/
		}

		//outputs_stuff[tstore] = atof(toks);          // convert char string to double
		toks = strtok (NULL, " ,\n,\t");         // continue on to next token
	}
	
	free (buffer);
}

//#########################################################################################
// FILE LOADING FUNCTIONS
//#########################################################################################

// eg, input_standardGA
// this function is called in main() directly to load in the main input file
// best to look at the 3DPV documentaion for more information
void ReadMainInputFile (string mainInputFile, int mode)  // takes pointer to first element in new coordinate array
{	
	ifstream ifs(mainInputFile.c_str());

	char stdInputs[numStdInputs][100];

	for(int i = 0; i < numStdInputs; i++ ) {
		ifs.getline(stdInputs[i], 100);
	}

	globalCellsOneSided	= atoi(stdInputs[0]);
	numCells			= atoi(stdInputs[1]);
	raysPerCell			= atoi(stdInputs[2]);
	xBound				= atof(stdInputs[3]);
	yBound				= atof(stdInputs[4]);
	zBound				= atof(stdInputs[5]);
	startMonth			= atoi(stdInputs[6]);
	startDay			= atoi(stdInputs[7]);
	startYear			= atoi(stdInputs[8]);
	endMonth			= atoi(stdInputs[9]);
	endDay				= atoi(stdInputs[10]);
	endYear				= atoi(stdInputs[11]);
	longitude 			= atof(stdInputs[12]);
	latitude    		= atof(stdInputs[13]);
	avgTemp     		= atof(stdInputs[14]);
	avgPressure 		= atof(stdInputs[15]);
	elevation   		= atof(stdInputs[16]);
	timeZone	  		= atoi(stdInputs[17]);

	if(mode == 0) { //GA
		char GAInputs[numGAInputs][100];

		for(int i = 0; i < numGAInputs; i++ ) {
			ifs.getline(GAInputs[i], 100);
		}

		GA_population			= atoi(GAInputs[0]);
		GA_maxgenerations		= atoi(GAInputs[1]);
		GA_replace_ratio		= atof(GAInputs[2]);
		GA_select_tourneysize	= atoi(GAInputs[3]);
		GA_cross_probability	= atof(GAInputs[4]);
		GA_mutation_probability	= atof(GAInputs[5]);
		GA_willload_population	= atoi(GAInputs[6]);
		allow_area_restrict		= atoi(GAInputs[7]);
		GA_numDistinctModels	= atoi(GAInputs[8]);

		if(GA_willload_population) {
			char modelDirectory[100];

			for(int i = 0; i < GA_numDistinctModels; i++) {
				ifs.getline(modelDirectory, 100);
				strcpy(Input_modelALL[i], modelDirectory);
			}
		}
	} else if(mode == 1) { //fixed shape
		ifs.getline(Input_modelALL[0], 100);
	} else {
		ifs.getline(Input_modelALL[0], 100);
		string s(Input_modelALL[0]);

		modelSpecDir = s.substr(0, (s.length()-4)) + "_c.txt";
	}

	ifs.close();
}

//########################
// WRITE THE OUTPUT FILE THAT SPECIFIES THE ORIGINAL INPUT AS WELL AS OTHER THINGS
//#########################
// WriteMainOutputFile() is called to store the input file for use later by other modes
// of "solar3d" operation, like "grid"
void WriteMainOutputFile()  // this will be appendable, but this first entry is a copy of the values in the cmd-line-entered input file
{
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat = strcat(Clone_Dir, OutMain_name);

	FILE *OUTfile = fopen(output_cat, "w");
	
	//fprintf(OUTfile, "COMMAND LINE INPUTS:\n");
	fprintf(OUTfile, "number_generations\t%i\n", GA_maxgenerations);
	fprintf(OUTfile, "population_replace_ratio\t%f\n", GA_replace_ratio);
	fprintf(OUTfile, "selection_tournament_size\t%i\n", GA_select_tourneysize);
	fprintf(OUTfile, "crossover_probability\t%f\n", GA_cross_probability);
	fprintf(OUTfile, "mutation_probabilty\t%f\n", GA_mutation_probability);
	fprintf(OUTfile, "number_solarcells\t%i\n", numCells);
	fprintf(OUTfile, "raytraces_per_cell\t%i\n", raysPerCell);
	fprintf(OUTfile, "maximum_xposition\t%f\n", xBound);
	fprintf(OUTfile, "maximum_yposition\t%f\n", yBound);
	fprintf(OUTfile, "maximum_zposition\t%f\n", zBound);
	fprintf(OUTfile, "sim_month\t%i\n", startMonth);
	fprintf(OUTfile, "sim_day\t%i\n", startDay);
	fprintf(OUTfile, "sim_year\t%i\n", startYear);
	fprintf(OUTfile, "population_size\t%i\n",GA_population);
	// stuff below here is printed for reference, but need not be loaded back in program when using "view" option
	if (GA_willload_population==1) {
		for (int mdd=0; mdd < GA_numDistinctModels; mdd++)
		{
			fprintf(OUTfile, "model %i: %s\n", mdd, Input_modelALL[mdd]);
		}
	}
	//fprintf(OUTfile, "---------------------------------------------\n");

	fclose(OUTfile); 
}
//#####################################################################
// READ IN CELL VERTICES
//#####################################################################
// called in main() if loading a population. called just once to count how many vertices are in
// the models. then, ReadStructureVertices_Retrieve() is called once for each model specified in the input file
int ReadStructureVertices_Count (int model_ind) 
{
  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;
  strcpy(Clone_Dir, Main_Dir);
	//printf("will count in: %s\n", strcat(Clone_Dir, Input_modelALL[model_ind]));
	//printf("%s\n", Clone_Dir);
	pFile = fopen (/*Clone_Dir*/strcat(Clone_Dir, Input_modelALL[model_ind]), /*"rb"*/"r");

  if (pFile==NULL) {fputs ("File error, for model file used to count vertices\n",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
//  if (result != lSize) {fputs ("Reading error",stderr); exit (3);} // check gives integer comparison warning

  // the whole file is now loaded in the memory buffer.

  // terminate
  fclose (pFile);
  
  int tsize=0;
  
  // FIRST PASS: count the data (increment tsize for each token)
  char *toks;
  toks = strtok(buffer, " ,\n");            // break file contents into tokens delimited by commas and whitespaces
  while (toks!=NULL)
  {
	  toks = strtok (NULL, " ,\n");         // go to next token
	  ++tsize;                              // add to token (ie coordinate) count
  }
  //cout << tsize << endl;                    // eg, if 480 vertices in model this should be 480x3 = 1440
  
  free (buffer);
  return tsize;                             // return number of tokens (coordinates) in file
}

// ReadStructureVertices_Retrieve(): see comment for previous function
// this function parses a (.txt) model file of format described in the 3DPV documentation
// and then stores the coordinates of the triangular solar cells in an array pointed to in arguments (<coordinates>)
void ReadStructureVertices_Retrieve (double *coordinates, int model_ind)  // takes pointer to first element in new coordinate array
{	
	FILE * pFile;
	long lSize;
	char * buffer;
	size_t result;
	strcpy(Clone_Dir, Main_Dir); printf("%s\n",Input_modelALL[model_ind]);

	pFile = fopen (strcat(Clone_Dir, Input_modelALL[model_ind]), "rb");

	if (pFile==NULL) {fputs ("File error, for model file",stderr); exit (1);}
	
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,pFile);
//	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}   // check gives integer comparison warning

	// the whole file is now loaded in the memory buffer.

	// terminate
	fclose (pFile);
	
	int tstore=0;
	
	// SECOND PASS: store data as array of doubles
	char *toks;
	toks = strtok(buffer, " ,\n");            // break file contents into tokens delimited by commas and whitespaces
	while (toks!=NULL)	
	{
		coordinates[tstore] = atof(toks);          // convert char string to double
		toks = strtok (NULL, " ,\n");         // continue on to next token
		++tstore;                             // next array index for following iteration
	}
	
	free (buffer);
}

//#######################################################################################
// READ A PREVIOUSLY-CREATED  OUTPUT FILE FOR "GRID" OPTION OF PROGRAM
//#######################################################################################
// formerly used just for defunct "view" option of solar3d program
// this reads the input parameters written to an output directory of some previously-run
// "sim" Genetic algorithm run of solar3d.
void ReadMainOutputFile (char *output_mainfile)  // takes pointer to first element in new coordinate array
{	
	FILE * pFile;
	long lSize;
	char * buffer;
	size_t result;

	pFile = fopen ( output_mainfile, "rb" );
	if (pFile==NULL) {fputs ("pre-view loading : File error, main output file\n",stderr); exit (1);}
	
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,pFile);
//	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}   // check gives integer comparison warning

	// the whole file is now loaded in the memory buffer.

	// terminate
	fclose (pFile);
	
	int tstore=0;
	
	// SECOND PASS: store data as array of doubles
	char *toks;
	toks = strtok(buffer, " ,\n,\t");            // break file contents into tokens delimited by commas and whitespaces
	while (toks!=NULL)	
	{
		switch(tstore)
		{
			case 0: break; // the default do-nothing value
			case 1: GA_maxgenerations=atoi(toks); tstore=0; break;
			case 2: GA_replace_ratio=atof(toks); tstore=0; break;
			case 3: GA_select_tourneysize=atoi(toks); tstore=0; break;
			case 4: GA_cross_probability=atof(toks); tstore=0; break;
			case 5: GA_mutation_probability=atof(toks); tstore=0; break;
			case 6: numCells=atoi(toks); tstore=0; break;
			case 7: raysPerCell=atoi(toks); tstore=0; break;
			case 8: xBound=atof(toks); tstore=0; break;
			case 9: yBound=atof(toks); tstore=0; break;
			case 10: zBound=atof(toks); tstore=0; break;
			case 11: startMonth=atoi(toks); tstore=0; break;
			case 12: startDay=atoi(toks); tstore=0; break;
			case 13: startYear=atoi(toks); tstore=0; break;
			case 14: GA_population=atoi(toks); tstore=0; break;
			default: tstore=0; break;
		}
		if (strcmp("number_generations",toks)==0) {tstore=1;}    // set identity for value on next token
		if (strcmp("population_replace_ratio",toks)==0) {tstore=2;}
		if (strcmp("selection_tournament_size",toks)==0) {tstore=3;}
		if (strcmp("crossover_probability",toks)==0) {tstore=4;}
		if (strcmp("mutation_probability",toks)==0) {tstore=5;}
		if (strcmp("number_solarcells",toks)==0) {tstore=6;}
		if (strcmp("raytraces_per_cell",toks)==0) {tstore=7;}
		if (strcmp("maximum_xposition",toks)==0) {tstore=8;}
		if (strcmp("maximum_yposition",toks)==0) {tstore=9;}
		if (strcmp("maximum_zposition",toks)==0) {tstore=10;}
		if (strcmp("sim_month",toks)==0) {tstore=11;}
		if (strcmp("sim_day",toks)==0) {tstore=12;}
		if (strcmp("sim_year",toks)==0) {tstore=13;}
		if (strcmp("population_size",toks)==0) {tstore=14;}
		
		toks = strtok (NULL, " ,\n,\t");         // continue on to next token
	}
	
	free (buffer);
}

void ReadRoofPositions(double in_pos[], int nextline_num, int i_am_best)
{
	strcpy(Clone_Dir, Unique_Dir);
	FILE * posFile = fopen ( strcat(Clone_Dir,  "RoofTestPositions"), "r");

	long lSize;
	char * buffer;
	size_t result;
	//strcpy(Clone_Dir, Main_Dir);
	
	//printf("Check Roof a\n");

	// obtain file size:
	fseek (posFile , 0 , SEEK_END);
	lSize = ftell (posFile);
	rewind (posFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,posFile);
//	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}   // check gives integer comparison warning

	// the whole file is now loaded in the memory buffer.
	//printf("Check Roof b\n");
	// terminate
	fclose (posFile);
	
	int temp_store_a=0; // other simpler method
	int temp_store_b=0;
	
	// SECOND PASS: store data as array of doubles
	char *toks;
	toks = strtok(buffer, " ,\n,\t");            // break file contents into tokens delimited by commas and whitespaces
	
	while (toks!=NULL)	
	{
		if (temp_store_a == i_am_best)
		{
			in_pos[temp_store_b] = atof(toks);// char string to double
		}
		//printf("Check Roof %i %i\n", temp_store_a, temp_store_b );
		
		++temp_store_b;
		if (temp_store_b == nextline_num)
		{
			temp_store_a++;
			temp_store_b=0;
		}
		toks = strtok (NULL, " ,\n,\t");        	 // continue on to next token
	}	
	free (buffer);
	// note: since cellpower file contains all individuals of a generation, I have to sum up all the datas for each individual in the
	// final generation to see which power is highest because I did not save any kind of index that.. or compare number to what is 		// found in the  bestFit file
}

//#####################################################################
// WRITE INPUT FILE FOR GENETIC ALGORTIHM
//#####################################################################
// WriteGAInputFile() is called in main() when the GA "sim" mode of solar3d is being
// performed. This writes automatically an input file for the Genetic algorithm code
// based on inputs written both in the command-line input file and some global variables at top of main.cpp
// See 3DPV documentation and/or Illigal Genetic Algorithms toolkit 2007 documentation
// for more information on the types of inputs and possibilities with the GA code
void WriteGAInputFile(int num_cells, int num_params)
{
	// argument numParameters is number of "decision variables" to be passed to GA (coordinates for now?)
	
	int num_dec = num_cells*num_params; printf("%i\n", num_dec);
	
	strcpy(Clone_Dir, Unique_Dir);
	char *Input_cat = strcat(Clone_Dir, Input_GA);
	FILE *GAfile = fopen(Input_cat, "w");
	if (GA_numberof_objectives==1) {
		fprintf(GAfile, "SGA # type\n");                                      // GA type SGA or NSGA
	}
	else {
		fprintf(GAfile, "NSGA # type\n");                                      // GA type SGA or NSGA
	}
	fprintf(GAfile, "%i\n", num_dec);                        // # of "decision variables"
	for (int dv=0; dv < num_cells; ++dv)
	{	
		// variable ranges: min, max : for the 9-coordinate-per-cell picture
		for (int trip=0; trip < 3; trip++) {
			fprintf(GAfile, "double %f %f\n", /*0.0*/Footprint_xmin, xBound);
			fprintf(GAfile, "double %f %f\n", /*0.0*/Footprint_ymin, yBound);
			fprintf(GAfile, "double %f %f\n", Structure_zmin, zBound);
		}
	}
	if (GA_numberof_objectives==2)
	{
		fprintf(GAfile, "2\n");                                        // # of objectives
		fprintf(GAfile, "Max\n");                                      // objective type
		fprintf(GAfile, "Max\n");
	}
	if (GA_numberof_objectives==1)
	{
		fprintf(GAfile, "1\n");                                        // # of objectives
		fprintf(GAfile, "Max\n");                                      // objective type
	}
	fprintf(GAfile, "0\n");                                        // # of constraints
	// put one line here for each constraint weight (eg 1.0)
	fprintf(GAfile, "%i\n", GA_population);                        // # in poopulation
	fprintf(GAfile, "%i\n", GA_maxgenerations);                    // # of generations to cycle through
	fprintf(GAfile, "%f\n", GA_replace_ratio);                     // proportion of each population to replace?
	
	fprintf(GAfile, "NoNiching\n");                                // type for niching of solutions
	//fprintf(GAfile, "TournamentWOR %f\n", GA_select_tourneysize);  // type and param's for SELECTION operator
	fprintf(GAfile, "TournamentWOR %i\n", GA_select_tourneysize);  // type and param's for SELECTION operator
	//fprintf(GAfile, "RouletteWheel\n");
	//fprintf(GAfile, "SUS\n");										// this is preferred over the roulette method
	
	fprintf(GAfile, "%f\n", GA_cross_probability);                 // probability of a cross between two parents
	switch (GA_cross_type) {
	case ONEPOINT: 
		fprintf(GAfile, "OnePoint\n");
		break;
	case TWOPOINT:
		fprintf(GAfile, "TwoPoint\n");                                 // crossover, twopoint i think should be chunks swapped
		break;
	case UNIFORM:
		fprintf(GAfile, "Uniform %f\n", GA_genewise_swap_probability); // type of cross, and probability per gene to swap 
		break;
	default: break;
	}
	
	fprintf(GAfile, "%f\n", GA_mutation_probability);              // probability of mutation
	//fprintf(GAfile, "Polynomial 10\n");                            // type of mutation method and parameter values
	fprintf(GAfile, "Genewise ");
	for (int msig=0; msig<9*num_cells; msig++)
	{
		fprintf(GAfile, "1 "); 
		//printf("%i\n", msig); // test
	}
	fprintf(GAfile, "\n");
	fprintf(GAfile, "NoScaling\n");                                // set scaling type 
	fprintf(GAfile, "NoConstraints\n"/*"Tournament\n"*/);         // constraint-handling method (could also choose "Penalty Linear"
	fprintf(GAfile, "NoLocalSearch\n");                            // local serach method (could also choose SimplexSearch [parameters]
	
	fprintf(GAfile, "0\n");                                        // number of stopping criteria
	strcpy(Clone_Dir, Unique_Dir);	
	if (GA_willload_population==0) {							   // whether to load population file
		// no, do not load from file
		fprintf(GAfile, "0\n");
	}
	else {
		// yes, load from file
		fprintf(GAfile, "1 %s 1\n", strcat(Clone_Dir, GA_initial_population_file) );       // load file?, which file?, evaluate population
	}
	
	// stopping criteria, for now I will not add any beyond the previous statement of maximum generations
	strcpy(Clone_Dir, Unique_Dir);	
	char *tmp_sol = strcat(Clone_Dir, Output_GA_sol);
	if (GA_willsave_solutions) {
		fprintf(GAfile, "%i %s\n", GA_willsave_solutions, tmp_sol);
	}
	else {
		fprintf(GAfile, "%i\n", GA_willsave_solutions);
	}
	fprintf(GAfile, "#END");
	
	fclose(GAfile);                                                // done 
}

//#####################################################################
// SETUP CELL VERTICES 
//#####################################################################
// SetupCellVertices() is one of the first 3DPV functions written. it was formerly very long because we were still thinking about
// how to approach the problem and set constraints, and so there would be these "3d solar cells" composing the overall 
// structure. These would have a shape like a cube or tetrahderon or ellipse. However, eventually I thought it more general
// to make no 3d shape assumption and just make all the cells of an arbitrary-triangle. 
// Presently, this function just initializes a single solar cell geometrically, though it is rather cumbersome in the array format...
void SetupCellVertices(Cell_3 *cel, int shape_type, int cellindex)  // note cell_index mainly for MODEL
{	
	int v_index_triflat[2][3] = {{0,1,2},{0,2,1}};  // stores integer indices for 3 verts of each triangle
	//int v_index_tetrahedron[4][3] = {{0,1,2},{0,3,1},{0,2,3},{3,2,1}};
	//int v_index_rectflat[4][3] = {{0,3,1},{0,2,3},{3,0,1},{3,2,0}};
	//int v_index_cuboid[12][3] = {{0,1,6},{0,6,5},{1,2,6},{2,7,6},{2,3,4},{3,4,7},{3,0,5},{3,5,4},{5,6,4},{6,7,4},{1,0,3},{2,1,3}};
	Triangle *temptri;
	double pn[3]; 
	double newvertex1[3];
	double newvertex2[3];
	double newvertex3[3];
	//int cellindex = (Cell_3::cell_counter-1);  // for MODEL type to obtain correct vertices
	// note: I don't know why cellindex was like this, I guess it was for old GenerateCell stuff
	// so when SetupVertices was called the cell_counter was not maxed yet
	
	switch(shape_type) { // the use of a switch statement is no longer necessary
	case MODEL:
		// always made of flat triangles
		cel->numVertices = 3;
		cel->numTriangles = 2; // like TRIFLAT, it has a top and bottom triangle
//		cel->numEdges = 3;




		for (int nt=0; nt < cel->numTriangles; nt++)
		{
			cel->triangleList[nt] = new Triangle;
		}
//		for (int ne=0; ne < cel->numEdges; ne++)
//		{
//			cel->edgeList[ne] = new Edge;
//		}
		cel->vertices[0] = new Position;
		cel->vertices[0]->x = modelPointer[cellindex*9+0]; // get coordinate from global pointer to array modelCoordinates
		cel->vertices[0]->y = modelPointer[cellindex*9+1];
		cel->vertices[0]->z = modelPointer[cellindex*9+2];
		
		cel->vertices[1] = new Position;
		cel->vertices[1]->x = modelPointer[cellindex*9+3];
		cel->vertices[1]->y = modelPointer[cellindex*9+4];
		cel->vertices[1]->z = modelPointer[cellindex*9+5];
		
		cel->vertices[2] = new Position;
		cel->vertices[2]->x = modelPointer[cellindex*9+6];
		cel->vertices[2]->y = modelPointer[cellindex*9+7];
		cel->vertices[2]->z = modelPointer[cellindex*9+8];
		
		// for former complicated shapes there could be many more triangles. but a flat cell has two sides-->numTriangles=2
		for (int m=0; m < (cel->numTriangles); ++m)
		{
			temptri = cel->triangleList[m];
			// deativate STUFF ONLY VALID FOR CLOSED SHAPES!!!!!!! for the most part
			//if (m==1) {temptri->is_deactivated=true;} else {temptri->is_deactivated=false;}
			temptri->index[0] = v_index_triflat[m][0];
			temptri->index[1] = v_index_triflat[m][1];
			temptri->index[2] = v_index_triflat[m][2];

			//cal the normal to a particular face
			pn[0]=0; pn[1]=0; pn[2]=0;
			newvertex1[0] = (cel->vertices[temptri->index[0]])->x;
			newvertex1[1] = (cel->vertices[temptri->index[0]])->y;
			newvertex1[2] = (cel->vertices[temptri->index[0]])->z;
			newvertex2[0] = (cel->vertices[temptri->index[1]])->x;
			newvertex2[1] = (cel->vertices[temptri->index[1]])->y;
			newvertex2[2] = (cel->vertices[temptri->index[1]])->z;
			newvertex3[0] = (cel->vertices[temptri->index[2]])->x;
			newvertex3[1] = (cel->vertices[temptri->index[2]])->y;
			newvertex3[2] = (cel->vertices[temptri->index[2]])->z;
			plane_exp_normal_3d(newvertex1, newvertex2, newvertex3, pn);
			temptri->Normal[0]=pn[0]; temptri->Normal[1]=pn[1]; temptri->Normal[2]=pn[2];
			
			temptri->vertList[0] = newvertex1[0]; temptri->vertList[1] = newvertex1[1]; temptri->vertList[2]=newvertex1[2]; 
			temptri->vertList[3] = newvertex2[0]; temptri->vertList[4] = newvertex2[1]; temptri->vertList[5]=newvertex2[2];
			temptri->vertList[6] = newvertex3[0]; temptri->vertList[7] = newvertex3[1]; temptri->vertList[8]=newvertex3[2];
					
			temptri->Area = triangle_area_3d(temptri->vertList); // call to geometry.cpp function to get triangles area
			
			if (m==0) 
			{
				totalSurfaceArea+=temptri->Area;  // get the total area
				temptri->is_positive=true;
			} // finally, add area to the total structure area (m==0 => postive norm side only)
			if (m==1) {temptri->is_positive=false;}	
		}
		break;

	default: break;
	}
}

//##################################################################################
// FIND AND SPECIFY EDGES SHARED BY TRIANGLES THAT COMPOSE the solar cells' geometry
//##################################################################################
// this was more necessary for complicated 3d base-shapes of solar cells (ie, not just triangles)
// but it is actually still used, called right after SetupCellVertices for a solar cell.
// It is probably fine to just leave it as it is, though now (2010) it's been a long time since I went
// through this function. The reason for the edge approach was when I was still deciding
// how to calculate shadows cast by solar cells on others. The present method is just ray-tracing
// from a large number of discrete points, but another former method was to calculate a shadow-volume
// projected down from a 3d object in the direction of the sunlight. Then a particular other
// cell would be checked to lie within this. 
//int BuildEdges(Cell_3 *cel)
//{
//
//	int triangleCount = cel->numTriangles;
//
//	int edgeCount = 0;
//	int ind=0;
//	Edge *edge;
//
//	for (int a = 0; a < triangleCount; a++)
//	{
//		Triangle *triangle = cel->triangleList[a];
//
//		int i1 = triangle->index[0];
//		int i2 = triangle->index[1];
//		int i3 = triangle->index[2];
//
//		// cyclic order
//		if (i1 < i2)
//		{
////			cout << "1,2 " << ind << endl;
//			edge = cel->edgeList[ind];
//			edge->vertexIndex[0] = i1;
//			edge->vertexIndex[1] = i2;
//			edge->triangleIndex[0] = a;
//			edge->triangleIndex[1]= -1;  // leave open spot
//			edgeCount++;
//			ind++;
//		}
//		if (i2 < i3)
//		{
////			cout << "2,3 " << ind << endl;
//			edge = cel->edgeList[ind];
//			edge->vertexIndex[0] = i2;
//			edge->vertexIndex[1] = i3;
//			edge->triangleIndex[0] = a;
//			edge->triangleIndex[1]= -1;  // leave open spot
//			edgeCount++;
//			ind++;
//		}
//		if (i3 < i1)
//		{
////			cout << "3,1 " << ind << endl;
//			edge = cel->edgeList[ind];
//			edge->vertexIndex[0] = i3;
//			edge->vertexIndex[1] = i1;
//			edge->triangleIndex[0] = a;
//			edge->triangleIndex[1]= -1;  // leave open spot
//			edgeCount++;
//			ind++;
//		}
//	}
//
//	// second pass: match triangles to edges
////	triangle = triangleArray;
//	for (int a = 0; a < triangleCount; a++)
//	{
//		Triangle *triangle = cel->triangleList[a];
////		cout << "Looking at triangle # " << a << endl;
//
//		int i1 = triangle->index[0];
//		int i2 = triangle->index[1];
//		int i3 = triangle->index[2];
//
//		// these matching triangles to edges must be anti-cyclic order
//		if (i1 > i2)
//		{
//			for (int b = 0; b < edgeCount; b++)
//			{
//				edge = cel->edgeList[b];
//				if ((edge->vertexIndex[0]==i2) && (edge->vertexIndex[1]==i1) && (edge->triangleIndex[1]==-1))
//				{
//					edge->triangleIndex[1] = a;
//					break; // don't search the rest of the edges after a match
//				}
//				//edge++;
//			}
//		}
//
//		if (i2 > i3)
//		{
//			for (int b = 0; b < edgeCount; b++)
//			{
//				edge = cel->edgeList[b];
//				if ((edge->vertexIndex[0]==i3) && (edge->vertexIndex[1]==i2) && (edge->triangleIndex[1]==-1))
//				{
//					edge->triangleIndex[1] = a;
//					break; // don't search the rest of the edges after a match
//				}
//				//edge++;
//			}
//		}
//
//		if (i3 > i1)
//		{
//			//edge = *(cel->edgeList);
//			for (int b = 0; b < edgeCount; b++)
//			{
//				edge = cel->edgeList[b];
//				if ((edge->vertexIndex[0]==i1) && (edge->vertexIndex[1]==i3) && (edge->triangleIndex[1]==-1))
//				{
//					edge->triangleIndex[1] = a;
//					break; // don't search the rest of the edges after a match
//				}
//				//edge++;
//			}
//		}
//
//		//triangle++;
//	}
//	return (edgeCount);
//}

//###################################################################################
// CREATE STRUCTURE (FROM CELL CONFIGURATION), and the destroy function to go with it
//###################################################################################
void BuildStructure(Cell_3 structure_array[], int structure_component)
{
//	Cell_3::cell_shape = Cells_baseshape; // NOT USED ANYMORE
//	Cell_3::cell_radius = Cells_radius;   // NOT USED ANYMORE

	int true_num=0;
	if (Cells_baseshape==MODEL)
	{
		true_num = numCoordinates/9;  // this is dangerous to divide since both are integers, but 9 should be a factor always
	}
	else
	{
		true_num=numCells;
	}
	
	for (int i=0; i<true_num; i++)
	{
		SetupCellVertices( &structure_array[i], structure_component, i );  // pass address of the ith element
//		BuildEdges( &structure_array[i] );         // pass address of the ith element
	}

	//midpoint. I think this may be defunct presently...
//	if (gotcenter==0)
//	{
//		gotcenter=1;
//		midcommon[0]/= (numCells*3);
//		midcommon[1]/= (numCells*3);
//		midcommon[2]/= (numCells*3);
//	}
}

//##########################################################################
// SUN CURRENT POSITION (returns cartesian point with scaled-down distance)
//##########################################################################
// SunCurrentPosition() is called at each time step in newCalculatePowerDayCycle() to
// set up inputs and call the SunPos algorithm in spa.cpp. the return value is a position 
Position* SunCurrentPosition(double ctime, bool find_init = 0, int year = 0, int month = 0, int day = 0)
{
	// define spherical coordinates
	double theta = 0;
	double phi = 0;
	double r = 0;
	
	// ctime is a fractional hour --> standard hour, minute, second values
	// conversion
	double hour = floor(ctime);
	double decimal = ctime - hour;
	double minute = decimal * 60;
	decimal = minute - floor(minute);
	minute = floor(minute);
	double second = floor(decimal*60);
	
	int errorcode;
	// first time, initialize all spa data structure variables
	if (find_init)
	{
		// enter initial inputs into data structure
		spa.year          = year;
		spa.month         = month;
		spa.day           = day;
		spa.hour          = (int)hour;             //ADDED 3/16/09 CONVERSIONS (INT) TO GET RID OF WARNINGS// 0 to 24
		spa.minute        = (int)minute;           // 0 to 59
		spa.second        = (int)second;           // 0 to 59
		spa.delta_t       = 32.184 + 33 -0.4; // see header file for link
		spa.timezone      = timeZone;   // negative hours west of GMT
		spa.longitude     = longitude;
		spa.latitude      = latitude;
		spa.elevation     = elevation;
		spa.pressure      = avgPressure;
		spa.temperature   = avgTemp;
		spa.slope         = 0;                   // surface slope from horizontal (-360 to 360)
		spa.azm_rotation  = 180;                 // surface azimuth rotation (like compass; 0 = south; 180=north i think)
		spa.atmos_refract = 0.5667;              // (-5 to 5 degrees)
		spa.function      = SPA_ZA_RTS;          // output zenith, azimuth, and time of sunrise/sunset (in fractional hour)
		errorcode     = spa_calculate(&spa); // carry out "sun position algorithm with inputs given for spa
		// if errorcode != 0 may want to crash program
	}
	else
	{
		// update spa variables which change during "day" cycle
		spa.hour          = (int)hour;                // 0 to 24
		spa.minute        = (int)minute;              // 0 to 59
		spa.second        = (int)second;              // 0 to 59
//		spa.azm_rotation  = 180;                 // surface azimuth rotation (like compass; 0 = south; 180=north i think)
		spa.function      = SPA_ZA_INC;          // output zenith, azimuth, and incidence angle
		errorcode     = spa_calculate(&spa); // carry out "sun position algorithm with inputs given for spa
	}
	
	// transform this data into spherical coordinates, scale distance to sun to 100m from the cell configuration
	theta = degtorad(spa.zenith);     
	phi = degtorad(90 - spa.azimuth); // choosing x = east, y = north, z = up
	r =1.49*pow(10.0,8.0) /*10000*/;                        // 100m --> 10,000cm = 10,000 units
	// NOTE: with the JUNE 09 change of parallel rays, the distance really doesn't matter now, as long as it
		// is a few times larger than the scale of the 3d structure
	//printf("%f, --> %f\n",spa.azimuth,phi); 
	//printf("theta= %f\n",spa.zenith);
	// return address of Position storing cartesian coordinate position of the sun
	Position sun_position;
	sun_position.x = r*sin(theta)*cos(phi);  // sin and cos require arg. in radians
	sun_position.y = r*sin(theta)*sin(phi);
	sun_position.z = r*cos(theta);
	/*return &sun_position;*/ // don't want to return this LOCALLY SCOPED THING
	Position *target_result;
	target_result = new Position;
	target_result->x = sun_position.x;
	target_result->y = sun_position.y;
	target_result->z = sun_position.z;
	return target_result;
}

////////////////////////////////////////////////////

// CalculateTransmittedPower() is called in newCalculatePowerDayCycle() to determine
// with some basic electromagnetism how much light is reflected. It returns the (double) transmitted
// power (which is taken to be "absorbed by the solar cell", and reflected power information
// is a value pointed to here. inputs are discussed below
// NOTE: the way I've been testing zero-relfectance consistencies is to set Tcoeff=1.
// NOTE: this transmitted power calculation for REFLECTED rays is done explicitly in newCalculatePowerDayCycle().
// This is a bit redundant, maybe but resulted from developing the code in stages...never got to cleaning it up.
double CalculateTransmittedPower(double net_area, double inc_angle, double *r_power_flux, double efficiency, double nRatio)
{
	// argument0 = the area of triangle face that is not shadowed, ie, in view of current "rays" (from sun)
	// argument1 = incidence angle is angle between face normal and "sun ray"
	// argument2 = pointer to a new reflected ray
	// calculates the transmitted (and hence reflected) power on a single triangle face of a solar cell

//reviewed on Dec 11, Fresnel eqn for parallel part, Bryan's one was wrong. Marco
//reviewed by MARCO on Jun 15th 2011
        double transm_angle = asin(nRatio*sin(inc_angle));	// snell's law...
	//Tcoeff=1; // for testing
	// NEW, for unpolarized light
       // Implement Fresnel Equations as in Wikipedia, using both incident and transmitted angle, and cosine instead of sine and/or tangent               
        double RP_check = (nRatio*cos(transm_angle)-cos(inc_angle))/(nRatio*cos(transm_angle)+cos(inc_angle));
        RP_check = RP_check*RP_check;
        double RS = (nRatio*cos(inc_angle)-cos(transm_angle))/(nRatio*cos(inc_angle)+cos(transm_angle));
        RS = RS*RS;
        double Rcoeff = (RP_check+RS)/2; // take the average for unpolarized light
        double Tcoeff = 1-Rcoeff;// new Tcoefficient
	///////////////////
// changed by Marco, using equation for correction of Air Mass. 11 Dec 2010
     // Need to see whether adding spa.zenith brute-force in this way will work
        double amtheta = degtorad(spa.zenith);
        double amcorrection = pow(1/cos(amtheta),0.678); //define exponent for am correction 
	double solar_variable = 1488*pow(0.7,amcorrection);  // then replace 2-> cos(spa.zenith!)
        //double absorbed_power = Tcoeff*solar_constant*net_area;
	// how much transmitted power X insolation X area of cell X efficiency
	double entire_power = (solar_variable*cos(inc_angle))*net_area;
	double absorbed_power = Tcoeff*entire_power*efficiency;//*Ws_to_kWhr;
	if (globalUseKWH==true) {absorbed_power*=Ws_to_kWhr;}
	double reflected_power_flux = Rcoeff*entire_power/*(net_area)*/; // want the power flux (W/m^2), not just the power.
	*r_power_flux = reflected_power_flux; //DID NOT MULTIPLY BY CELL EFFICIENCY YET!!! THAT IS DONE UPON other cell absorbing ray

	return absorbed_power;
//	return cos(inc_angle)*net_area;
}

///TRIANGLE RAY TRACE POINTS!!!
void GetSubTrianglesFromGrid(int subs, int subtri_ind[])
{
	 // subtri_ind[] has size = number of triangles needed is sub_num^2. x3 for storing each point index, 
	
	int row=0; // which row of points in trianlge (0th row is most points)
	int num_rows = subs+1; // just number of points in biggest row.
	int pts_in_row[num_rows];
	int tcnt=0; // count triangles
	int skipper=0; 
	for (row=0; row < num_rows; row++) // store num points per rows
	{
		pts_in_row[row] = num_rows-row; // so max number points for 0th row (base of triangle)
	}
	for (row=0; row<(num_rows-1); row++)
	{
		// first get the triangles with lower side base in row
		for (int pp=0; pp< pts_in_row[row]-1; pp++) // record one triangle for the first N-1 points in row with flat lower side
		{
			subtri_ind[tcnt*3 + 0] = (skipper) + pp;
			subtri_ind[tcnt*3 + 1] = (skipper) + pp+1;
			subtri_ind[tcnt*3 + 2] = (skipper) + pts_in_row[row]+pp; // gets the right point from next row
			tcnt++; // increment to next triangle
		}
		
		// second, get triangles with upper side base in the row
		// the 2nd to last row will skip this because 1 is not less than 2-1 so it never enters loop
		for (int pp=1; pp< pts_in_row[row]-1; pp++) // record one triangle for the first N-1 points in row with flat lower side
		{
			subtri_ind[tcnt*3 + 0] = (skipper) + pp;
			subtri_ind[tcnt*3 + 1] = (skipper) + pts_in_row[row]+pp; // gets right point from next row
			subtri_ind[tcnt*3 + 2] = (skipper) + pts_in_row[row]+pp-1; // gets the right point from next row
			tcnt++; // increment to next triangle
		}
		skipper += pts_in_row[row]; // skip the index ahead for next row by the number of points in the current row.
	}
}

// Get2DTriangleGrid() is called only once, before any simulations
// this function simply locates on a 1x1 unit square --> (half) triangle where to 
// put the gridpoints from where rays will be traced. then this mapping of
// points will get transformed onto each of the arbitrarily-shaped triangular cells 
// in the power calculation code.
void Get2DTriangleGrid()
{
	// the basic 2D triangle
	double tri2d[6] = {0.0,0.0, 1.0,0.0, 0.0,1.0};
	int grid_num_out=0;
	int subdiv_num= raysPerCell;
	
	// note sub_num not really equivalent to raytraces #, but i will just use it , in input file this raytrace # should be
	// approximately 4 for 15 rays, 5 for 21 rays and 6 for 28 rays...etc...

	// (triangle array of points[2*3], # subdivisions, # max grid points, *actual returned # grid points, returned array of grids
	triangle_gridpoints_2d ( tri2d, subdiv_num, maxNumGrids, &grid_num_out, globalTriGrid );

	// set the new number of raytraces to the actual number of gridpoints returned
	//raysPerCell = grid_num_out;
	//printf("Raytraces for %i gridpoints per triangle\n", raysPerCell);

	double tdiv_coord[subdiv_num*subdiv_num*6]; // number of triangles need is sub_num^2. x6 for 6 coordinates per 2d triangle
	int tdiv_ind[subdiv_num*subdiv_num*3];	// 3 indices (verts) per triangle times number of triangles
	//double tdiv_midpoint[subdiv*subdiv*2];	// 2 coordinates (x,y) per triangles times number of triangles


	GetSubTrianglesFromGrid(subdiv_num, tdiv_ind);

	// now finally store the coordinates based on indices of all the sub triangles
	for (int ts=0; ts < (subdiv_num*subdiv_num); ts++)
	{
		for (int tv=0; tv<3; tv++)
		{
			tdiv_coord[ts*6 + tv*2 + 0] = globalTriGrid[ 2*(tdiv_ind[ts*3 + tv]) + 0];
			tdiv_coord[ts*6 + tv*2 + 1] = globalTriGrid[ 2*(tdiv_ind[ts*3 + tv]) + 1];

			// here global part added 02/11/2010 for openGL new visuals, just make a copy of local variable
			globalTriVerts[ts*6 + tv*2 + 0] = tdiv_coord[ts*6 + tv*2 + 0];
			globalTriVerts[ts*6 + tv*2 + 1] = tdiv_coord[ts*6 + tv*2 + 1];
		}
	}

	// once all the triangles are known, I need the last step: the midpoint of each 2d triangle
	for (int ts=0; ts < (subdiv_num*subdiv_num); ts++)
	{
		for (int tv=0; tv<3; tv++)
		{
			globalTriMids[ts*2 + 0] += tdiv_coord[ts*6 + tv*2 + 0];
			globalTriMids[ts*2 + 1] += tdiv_coord[ts*6 + tv*2 + 1];
		} 
		globalTriMids[ts*2 + 0] /= 3;	// getting average x value of the points
		globalTriMids[ts*2 + 1] /= 3;	// getting average y value of the points
	}
	
	raysPerCell = subdiv_num*subdiv_num;
	printf("Raytraces for %i gridpoints per triangle\n", raysPerCell);
}	

// this is the part of the gridmesh onto triangle which is actually called to 
// transform a specific  point by using the pregenerated gridmesh points
// which all have equal areas assigned to them (at least in the 2D unit plane presumably...)
void OrderedTrianglePoint3D(double tri_v[], double out_v[], int pt_index, bool willwrite)
{
	//double sample1 = globalTriGrid[pt_index*2 + 0];
	//double sample2 = globalTriGrid[pt_index*2 + 1];
	double sample1 = globalTriMids[pt_index*2 + 0]; // commented out the above on 4/3/09
	double sample2 = globalTriMids[pt_index*2 + 1];	// replaced it with this
	
	// pt_index is counter to know which grid point to get from precalculated set 
	out_v[0] = tri_v[0] + (tri_v[3]-tri_v[0])*sample1 + (tri_v[6]-tri_v[0])*sample2;
	out_v[1] = tri_v[1] + (tri_v[4]-tri_v[1])*sample1 + (tri_v[7]-tri_v[1])*sample2;
	out_v[2] = tri_v[2] + (tri_v[5]-tri_v[2])*sample1 + (tri_v[8]-tri_v[2])*sample2;

	if (willwrite==true) {
		strcpy(Clone_Dir, Unique_Dir);
		char *output_catgc = strcat(Clone_Dir, OutGridCoordinates);
		
 		FILE *fileSAMPLES = fopen(output_catgc,"a");
		fprintf(fileSAMPLES, "%f\t%f\t%f\n", out_v[0], out_v[1], out_v[2]);
		fclose (fileSAMPLES);

		// added 7/25/09 for writing subcell "master" file ////////////////////////////////
		globalDataSubCoords[ globalCountForRecordCoord ] = out_v[0];
		globalDataSubCoords[ globalCountForRecordCoord + 1 ] = out_v[1];
		globalDataSubCoords[ globalCountForRecordCoord + 2 ] = out_v[2];
		globalCountForRecordCoord+=3; // for next time be on the next available slot
		/////////////////////////////////////////////////////////////////////////////////
	}
	
	///// HERE ADDED 02/11/2010 for new OpenGL visuals
	// pt_index is counter to know which grid point to get from precalculated set 
	// this actually gets the 9 coordinates of the gridpoint triangles instead of just the midpoint
	// these are only used so far for making the movie showing intensity of power on each cell over the day.
	double trisampleA0 = globalTriVerts[pt_index*6 + 0*2 + 0];
	double trisampleA1 = globalTriVerts[pt_index*6 + 0*2 + 1];
	double trisampleB0 = globalTriVerts[pt_index*6 + 1*2 + 0];
	double trisampleB1 = globalTriVerts[pt_index*6 + 1*2 + 1];
	double trisampleC0 = globalTriVerts[pt_index*6 + 2*2 + 0];
	double trisampleC1 = globalTriVerts[pt_index*6 + 2*2 + 1];

	double svt_v[9];
	// first vertex
	svt_v[0] = tri_v[0] + (tri_v[3]-tri_v[0])*trisampleA0 + (tri_v[6]-tri_v[0])*trisampleA1;
	svt_v[1] = tri_v[1] + (tri_v[4]-tri_v[1])*trisampleA0 + (tri_v[7]-tri_v[1])*trisampleA1;
	svt_v[2] = tri_v[2] + (tri_v[5]-tri_v[2])*trisampleA0 + (tri_v[8]-tri_v[2])*trisampleA1;
	// second vertex
	svt_v[3] = tri_v[0] + (tri_v[3]-tri_v[0])*trisampleB0 + (tri_v[6]-tri_v[0])*trisampleB1;
	svt_v[4] = tri_v[1] + (tri_v[4]-tri_v[1])*trisampleB0 + (tri_v[7]-tri_v[1])*trisampleB1;
	svt_v[5] = tri_v[2] + (tri_v[5]-tri_v[2])*trisampleB0 + (tri_v[8]-tri_v[2])*trisampleB1;
	// third vertex
	svt_v[6] = tri_v[0] + (tri_v[3]-tri_v[0])*trisampleC0 + (tri_v[6]-tri_v[0])*trisampleC1;
	svt_v[7] = tri_v[1] + (tri_v[4]-tri_v[1])*trisampleC0 + (tri_v[7]-tri_v[1])*trisampleC1;
	svt_v[8] = tri_v[2] + (tri_v[5]-tri_v[2])*trisampleC0 + (tri_v[8]-tri_v[2])*trisampleC1;

	// write to file "outSubcellTriCoords_best", a list of all subcell coordinates (not midpoints but actual triangles)
	if (willwrite==true) {
		strcpy(Clone_Dir, Unique_Dir);
		char *output_catgc = strcat(Clone_Dir, OutTriGridCoords);
		
 		FILE *fileSAMPLES = fopen(output_catgc,"a");
		fprintf(fileSAMPLES, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 	svt_v[0], svt_v[1], svt_v[2], 
										svt_v[3], svt_v[4], svt_v[5], 											svt_v[6], svt_v[7], svt_v[8]);
		fclose (fileSAMPLES);
	}
	////////////////
}

// RANDOMTRIANGLEPOINT3D
// takes array of 9 doubles as coordinates of the 3d triangle, and also pointer to first element for 3-dim array
// creates a random point on a 3d solar cell from which a ray would be traced
// Note, this is an older function not being used anymore, because the random point approach is unreliable for
// obvious reasons, especially with not enough gridpoints
void RandomTrianglePoint3D(double tri_v[], double out_v[], bool willwrite)
{
	double sample1=0;
	double sample2=0;
	
	// sample a random point in the unit square [0,1]
	do {
		sample1 = myRandom.random01();
		sample2 = myRandom.random01();
	} while ((sample1+sample2) > 1);    // want triangle not whole parallelogram
	
	// from wolfram mathworld: "triangle point picking"
	// correct mapping using triangle edge vectors
	out_v[0] = tri_v[0] + (tri_v[3]-tri_v[0])*sample1 + (tri_v[6]-tri_v[0])*sample2;
	out_v[1] = tri_v[1] + (tri_v[4]-tri_v[1])*sample1 + (tri_v[7]-tri_v[1])*sample2;
	out_v[2] = tri_v[2] + (tri_v[5]-tri_v[2])*sample1 + (tri_v[8]-tri_v[2])*sample2;
	
	if (willwrite==true) {
 		FILE *fileSAMPLES = fopen("/home/bryan/SGA_Solar/samplexyz","a");
		fprintf(fileSAMPLES, "%f\t%f\t%f\n", out_v[0], out_v[1], out_v[2]);
		fclose (fileSAMPLES);
	}
}

//assuming year1 <= year2, month1 <= month2.
bool isBefore(int year1, int month1, int day1, int year2, int month2, int day2) {
	if(year1 < year2) return true;
	if(year1 > year2) return false;

	if(month1 < month2) return true;
	if(month1 > month2) return false;

	if(day1 <= day2) return true;
	else return false;
}

bool isLeapYear(int year){
	if(year % 400 == 0) {
		return true;
	} else if(year % 100 == 0) {
		return false;
	} else if(year % 4 == 0) {
		return true;
	} else {
		return false;
	}
}

void nextDay(int &year, int &month, int &day) {
	switch(month) {
	case 1:
		if(day == 31) {
			month = 2;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 2:
		if(day == 29) {
			month = 3;
			day = 1;
		} else if(day == 28) {
			if(isLeapYear(year)) {
				day += 1;
			} else {
				month = 3;
				day = 1;
			}
		} else {
			day += 1;
		}

		break;
	case 3:
		if(day == 31) {
			month = 4;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 4:
		if(day == 30) {
			month = 5;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 5:
		if(day == 31) {
			month = 6;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 6:
		if(day == 30) {
			month = 7;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 7:
		if(day == 31) {
			month = 8;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 8:
		if(day == 31) {
			month = 9;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 9:
		if(day == 30) {
			month = 10;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 10:
		if(day == 31) {
			month = 11;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 11:
		if(day == 30) {
			month = 12;
			day = 1;
		} else {
			day += 1;
		}

		break;
	case 12:
		if(day == 31) {
			month = 1;
			day = 1;
			year += 1;
		} else {
			day += 1;
		}

		break;
	}
}

//##################################################################################3
// new function to calculate shadowing and power
//##################################################################################
// newCalculatePowerDayCycle() is really the main part of this code. It serves as the
// "fitness function" calculation for the genetic algorithm by determining how much
// power a given structure will collect in one day. In this problem, the fitness is very
// straightforward, it is just the power (or equivalently total energy in a day).
double newCalculatePowerDayCycle(Cell_3 cells_Array[], double xgenes[], double efficiencies[], double reflIndices[],
								 bool outputStructurePower = false, bool outputCellPower = false, string outDirectory = "")
{
	int curDate = startDay;
	int curMonth = startMonth;
	int curYear = startYear;
	double totalEnergy = 0;
	double totalEnergyFromRefl = 0;
	
	double time;
	clock_t start,finish;
	start = clock();

	ofstream structPowerFile, cellPowerFile;
	string structPowerDir, cellPowerDir;

	if(outputStructurePower) {
		structPowerDir = "struct_power";
	}

	if(outputCellPower) {
		cellPowerDir = "cell_power";
		//JIN: make directory
	}

	if ((Done_sundata==false) && (Delete_sundata==true)) {
		//remove ("/home/brosbam/Desktop/solar_results/optimization/sunXYZ");
		strcpy(Clone_Dir, Unique_Dir);

		char *output_cat = strcat(Clone_Dir, OutSun_Path);

		FILE *fileSUN= fopen(output_cat, "w");

		fprintf(fileSUN, "%f %f %f\n", 0.0, 0.0, 0.0); // origin
		fclose(fileSUN);
	}


	for (int rep1=0; rep1<5000001; rep1++) {
		globalDataSubCellPower_OneGuy[rep1]=0;
		SECONDglobalDataSubCellPower_OneGuy[rep1]=0;
		globalDataSubCellPower_Reflect[rep1]=0;
	}

	// first reset these cell powers to zero, since they are incremented, and not set, in below calculations
	// this one is OBSOLETE
	for (int nnn=0; nnn < numCells; nnn++) {
		globalDataCellPower_OneGuy[nnn] =0;
	}



	int qwik_skip=-1; // -1 means do not skip day energy calculation. positive integer value means same as that index
	for (int qi0=0; qi0<qwikIndMax; qi0++) {
		for (int qi1=0; qi1<numCells*9; qi1++) {
			if (xgenes[qi1] != qwikStoreCoords[qi0][qi1])
			{
				qi1=numCells*9;
				qwik_skip=-1;
			}
			else
			{
				qwik_skip=qi0;
			}
		}
		if (qwik_skip != -1 ) // if a complete match has been found, then no need to go further.
		{
			totalEnergy = qwikStoreEnergies[qwik_skip];		// set energy to be same as stored match
			totalEnergyFromRefl = qwikStoreReflEnergies[qwik_skip];	// set reflected energy to be same as stored match
			qi0=qwikIndMax;						// then force quit the FOR loop by setting max
		}
	}


	if (qwik_skip == -1) {

		while(isBefore(curYear, curMonth, curDate, endYear, endMonth, endDay)){ //JIN: looping over days
			if(outputStructurePower) {
				stringstream ss;
				ss << outDirectory << "/" << structPowerDir << "/" << curYear << "_" << curMonth << "_" << curDate << ".txt";
				structPowerFile.open(ss.str().data());
				//structPowerFile << setiosflags(ios::fixed) << setprecision(5);
			}

			int num_times=0; // counts the number of iterations that day is broken into
			int which_iteration=0; // like the above one but defined on 4/7/09 and incremented at end of iteration...not start

			SunCurrentPosition(6.0, 1, curYear, curMonth, curDate);
			double start_time = spa.sunrise + 0.1;  // add ~9 minutes to ensure sun is "totally risen"
			double end_time = spa.sunset - 0.1;     // subtract ~9 minutes to ensure

			double sun_pos[3];
		//	double sunDist=0;
			double curCellPower=0;
			double curCellReflPower=0;

			//Vector shadowRay1;
			//Vector shadowRay2;

			double lightDOTnorm1=0;
		//	double lightDOTnorm2=0;
		//	double pn[3];

		//	double difference[3];
			double flatArea=0;
		//	double ptDOTplane1=0

		//	double ptDOTplane3=0;
			Cell_3 *curr_cell;

			//double total_time = end_time-start_time; // unused variable
		//	Edge *curr_edge;
			Triangle *curr_triangle;

			map<int, double> CellSunDists2;

			double incidenceAngle=0;
			double totalDayEnergy=0;        // holds sum (integrated/averaged?? for each time interval)
			double structPower=0;    // holds sum of all power for just the current time interval (pos. of sun)
			double totalDayEnergyFromRefl=0;
			double totalInstantPower_fromRefl=0;

			for (double curTime=start_time; curTime<=end_time; curTime+=delta_time){ //JIN: looping over time intervals in one day
				if(outputCellPower) {
					stringstream ss;
					ss << outDirectory << "/" << cellPowerDir << "/" << curYear << "_" << curMonth << "_" << curDate << "_" << curTime << ".txt";
					cellPowerFile.open(ss.str().data());
					//cellPowerFile << setiosflags(ios::fixed) << setprecision(5);
				}

				// record sun position for 3d plot
				num_times++;
				int count_checkers=0; // for testing, check how many of the total shadowing cells are checked

				Position *tempsun;
				tempsun = SunCurrentPosition(curTime);      // update sun position for this time interval
				sun_pos[0]=tempsun->x; sun_pos[1]=tempsun->y; sun_pos[2]=tempsun->z;

				structPower=0;                                // reset the sum of all powers of cells over time interval
				totalInstantPower_fromRefl=0; //reset

				if (Done_sundata==false) {
					strcpy(Clone_Dir, Unique_Dir);

					char *output_cat = strcat(Clone_Dir, OutSun_Path);

					FILE *fileSUN= fopen(output_cat, "a");
					fprintf(fileSUN, "%f %f %f\n", sun_pos[0], sun_pos[1], sun_pos[2]);
					fclose(fileSUN);
				}

				double shadow_vert_dists[3];
				double shadow_least_dist=0;
				int shadow_which_least=0;

				double tri_vt[9]; // 9 coordinates of triangle
				double tri_s_vt[9];
				double rand_vt[3];
				Position *fpt_ptr;
				Triangle *curr_s_triangle;

				// univ_dir[] will be a direction vector of light rays from the sun to approximate the area source
				// a far distance away.
				double refer_pt[3];
				if (globalDoRoof==false)
				{
					refer_pt[0] = xBound/2; // these three set-variable lines could be outside this time loop
					refer_pt[1] = yBound/2;
					refer_pt[2] = zBound/2;
				}
				else
				{
					refer_pt[0] = globalRoofx/2;
					refer_pt[1] = globalRoofy/2;
					refer_pt[2] = zBound/2;
				}

				double refer_dist = points_dist_3d(sun_pos, refer_pt);
				double univ_dir[3];	// direction of parallel rays from sun
				univ_dir[0] = (refer_pt[0]-sun_pos[0])/refer_dist;
				univ_dir[1] = (refer_pt[1]-sun_pos[1])/refer_dist;
				univ_dir[2] = (refer_pt[2]-sun_pos[2])/refer_dist;


				/////////////////////////////////////////////////////////////////////////////////////

				// in a first pass through the all the cells in the current structure we are investigating
				// we do some basic dot product stuff to determine the direction relative to the sun and which sides
				// of the two-sided cells are illuminated.
				for (int cc=0; cc<numCells; cc++)
				{
					//int sPlaneIndex=0;
					curr_cell = &cells_Array[cc];
		//			curr_cell->numShadowingEdges = 0; //reset

					// in key equal to distance from sun, put value as index of cell (so distances are sorted)
		//			double cell_pos[3] = {curr_cell->position.x, curr_cell->position.y, curr_cell->position.z};
		//			sunDist = points_dist_3d(sun_pos, cell_pos);
					//CellSunDists2[cc] = sunDist; // 6/27/09 this could be commented out, not using ANY cellsun dists !!!



					for(int dd=0; dd<curr_cell->numTriangles; dd++) {
						curr_triangle = curr_cell->triangleList[dd];
						lightDOTnorm1 = r8vec_dot(3, univ_dir/*sun_unitv.vec*/, curr_triangle->Normal);
						curr_triangle->SunScalarProduct = lightDOTnorm1;
						if (lightDOTnorm1 < 0) {curr_triangle->is_facing_sun=true;}
						else {curr_triangle->is_facing_sun=false;}
					}


					// look-at-cell's-edges loop
		//			for (int dd=0; dd < (curr_cell->numEdges); dd++)
		//			{
		//				curr_edge = (curr_cell->edgeList[dd]);
		//
		//				curr_edge->is_shadower_edge=false; // reset, to be sure
		//
		//				// dot first triangle's norm with sun-to-cell vector, determine sign +/-
		//				curr_triangle = curr_cell->triangleList[curr_edge->triangleIndex[0]];
		//				lightDOTnorm1 = r8vec_dot(3, univ_dir/*sun_unitv.vec*/, curr_triangle->Normal);
		//				curr_triangle->SunScalarProduct = lightDOTnorm1;
		//				if (lightDOTnorm1 < 0) {curr_triangle->is_facing_sun=true;} else {curr_triangle->is_facing_sun=false;}
		//
		//				// dot other triangle's norm with sun-to-cell vector, determine sign +/-
		//			//	Triangle *other_triangle = curr_triangle;
		//				curr_triangle = curr_cell->triangleList[curr_edge->triangleIndex[1]];
		//				lightDOTnorm2 = r8vec_dot(3, univ_dir/*sun_unitv.vec*/, curr_triangle->Normal);
		//				curr_triangle->SunScalarProduct = lightDOTnorm2;
		//				if (lightDOTnorm2 < 0) {curr_triangle->is_facing_sun=true;} else {curr_triangle->is_facing_sun=false;}
		//			}
				}

				// second pass through all solar cells composing structure
				// this loop is the bulk of the code in this power function.
				// the general order of things is
				// 1) set the (ii)th cell to the "shadowed" cell that will will trace rays from
				// 2) look at the triangle (out of "front" or "back") which was determined previously to be facing sun
				// 3) set-up a mesh of gridpoints on this triangle, each for a ray to be traced in direction "univ_dir" seen above
				// 4) start to loop through all other cells, the "shadowing cells".
				// 5) with a shadowing cell assigned, loop through all of the gridpoints, raytrace to see if shadowing is blocking ray
				// 6) keep track of how many / which raytraces were blocked "tsample_blocked[rep]"
				// 7) when the shadowing for loop is on a new cell, don't need to trace rays seen to be blocked by other shadowing cells
				// 8) at the end of this, there is a kind of "effective area" of the cell now
				// 9) so, call the calculateTransmittedPower function to do some E&M computations with this information
				// 10) the transmitted power function also stores a reflected ray flux, so the next part of the code
				// 	repeats a lot of the previous stuff, though slightly different, to deal with the reflected parts of rays
				// 11) the rest is mainly summing and storing various quantities for later output, like distinguishing between
				// 	power contributions from once-reflected rays etc....
				//
				// NOTE: for adding further reflection steps, in principle the same thing could be done, but I'm sure
				// there is a much better way to do this that doesn't require bulky coding. maybe a generalization
				// to a function that allows an arbitrary number of reflections to be looked at (to some limit, due to computation time)

				for (int curCell=0; curCell<numCells; ++curCell) // JIN: looping through all cells
				{
					Cell_3 *shadowed_cell = &cells_Array[curCell];
					
					double efficiency = efficiencies[curCell];
					double nRatio = index_ref_air/reflIndices[curCell];

					// loop through each of the current victim cell's triangle faces
					for(int ff=0; ff < (shadowed_cell->numTriangles); ++ff)
					{
						curr_triangle = shadowed_cell->triangleList[ff];
						curr_triangle->numVerticesShadowed=0; // can reset here since above the shadowing cell loop

						// only check for shadowing if it is directly illuminated by sun in the first place
						// SECOND CONDITION COMMENTED OUT WHEN I WANT DOUBLE_SIDE COLLECTING CELLS
						if ((curr_triangle->is_facing_sun==true) and ((curr_triangle->is_positive==true) or (globalCellsOneSided==false)) )
						{

							//get the 3 vertex points for the current triangle
							// 5/08/09. which face of cell is important to know. because vertices were
							// defined in different orders, and for saving of subcell power to be
							// consistent, the gridpoints of each cell side must be in same order/position
							if (curr_triangle->is_positive==true)
							{

								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[0]];
								tri_vt[0]=fpt_ptr->x; tri_vt[1]=fpt_ptr->y; tri_vt[2]=fpt_ptr->z;
								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[1]];
								tri_vt[3]=fpt_ptr->x; tri_vt[4]=fpt_ptr->y; tri_vt[5]=fpt_ptr->z;
								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[2]];
								tri_vt[6]=fpt_ptr->x; tri_vt[7]=fpt_ptr->y; tri_vt[8]=fpt_ptr->z;


							}
							else
							{	// for back-facing triangle, the second two indices were flipped
								// direction of the surface normal is already known, so it doesn't
								// matter that the vertices are not CCW

								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[0]];
								tri_vt[0]=fpt_ptr->x; tri_vt[1]=fpt_ptr->y; tri_vt[2]=fpt_ptr->z;
								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[2]];
								tri_vt[3]=fpt_ptr->x; tri_vt[4]=fpt_ptr->y; tri_vt[5]=fpt_ptr->z;
								fpt_ptr = shadowed_cell->vertices[curr_triangle->index[1]];
								tri_vt[6]=fpt_ptr->x; tri_vt[7]=fpt_ptr->y; tri_vt[8]=fpt_ptr->z;


							}

							/*if (Write_vertdata==true)
							{
								FILE *fileVERT = fopen("/home/brosbam/Desktop/solar_results/optimization/sampleXYZ","a");
								fprintf(fileVERT, "%f %f %f\n", tri_vt[0], tri_vt[1], tri_vt[2]);
								fprintf(fileVERT, "%f %f %f\n", tri_vt[3], tri_vt[4], tri_vt[5]);
								fprintf(fileVERT, "%f %f %f\n", tri_vt[6], tri_vt[7], tri_vt[8]);
								fclose (fileVERT);
							}*/

							flatArea = curr_triangle->Area; /// calculated in cell construction since its constant throughout the "day"
							//  note: be careful with this rounding as (int) will be rounding down always


							int max_rayfires;
							if (ray_method==PERCELL) {
								max_rayfires = raysPerCell;
							}
							else {
								max_rayfires = (int)(rays_per_unit_area * flatArea);
							}
							double tsample_points[max_rayfires*3]; // store 3 coordinates for each point to ray-trace from
							int tsample_blocked[max_rayfires];     // 0 if ray is still unblocked, 1 if blocked
							double tsample_sundist[max_rayfires];
							double tsample_find_max=0;
							//double tsample_find_min=pow(10,15); // larger than any possible sun-cell distance

							Cell_3 *shadowing_cell;
							bool ray_in_tri = false;
							double ray_tri_pt[3]; // intersection point of line and triangle, if needed
							double raytri_dist; 	// distance from sun of this point, if needed
							int num_points_shadowed=0; // this will be between 0 and max_ray_fires



							for (int rep=0; rep<max_rayfires; rep++)
							{
								if (use_trigrid==0)
									{RandomTrianglePoint3D(tri_vt, rand_vt, Write_tsampledata);}
								if (use_trigrid==1)
									{OrderedTrianglePoint3D(tri_vt, rand_vt, rep, Write_tsampledata);}
								tsample_points[rep*3 + 0] = rand_vt[0];
								tsample_points[rep*3 + 1] = rand_vt[1];
								tsample_points[rep*3 + 2] = rand_vt[2];

								// should be one of the powers to this method: the memory of blocked points will endure
								// at least during this timestep, so other shadowing cells know not to check it
								tsample_blocked[rep] = 0;
								tsample_sundist[rep] = points_dist_3d(sun_pos, rand_vt);
								if (tsample_sundist[rep] > tsample_find_max) {tsample_find_max = tsample_sundist[rep];}
							//	if (tsample_sundist[rep] < tsample_find_min) {tsample_find_min = tsample_sundist[rep];}
							}
							// currently looping through all cells



							//looping through shadowing cells
							for (int kk=0; kk<numCells; ++kk)
							{
								shadowing_cell = &cells_Array[kk];

								if ( (kk!=curCell) /*&& (CellSunDists2[kk] < tsample_find_max)*/)
								{
									ray_in_tri = 0; // reset
									count_checkers++; // test
									curr_s_triangle = shadowing_cell->triangleList[0];
									if (curr_s_triangle->is_facing_sun == false) // assumes only two triangles, switch if wrong
									{
										curr_s_triangle = shadowing_cell->triangleList[1];
									}

									fpt_ptr = shadowing_cell->vertices[curr_s_triangle->index[0]];
									tri_s_vt[0]=fpt_ptr->x; tri_s_vt[1]=fpt_ptr->y; tri_s_vt[2]=fpt_ptr->z;
									fpt_ptr = shadowing_cell->vertices[curr_s_triangle->index[1]];
									tri_s_vt[3]=fpt_ptr->x; tri_s_vt[4]=fpt_ptr->y; tri_s_vt[5]=fpt_ptr->z;
									fpt_ptr = shadowing_cell->vertices[curr_s_triangle->index[2]];
									tri_s_vt[6]=fpt_ptr->x; tri_s_vt[7]=fpt_ptr->y; tri_s_vt[8]=fpt_ptr->z;

									double svert1[] = {tri_s_vt[0], tri_s_vt[1], tri_s_vt[2]};
									double svert2[] = {tri_s_vt[3], tri_s_vt[4], tri_s_vt[5]};
									double svert3[] = {tri_s_vt[6], tri_s_vt[7], tri_s_vt[8]};
									shadow_vert_dists[0] = points_dist_3d(sun_pos, svert1); // added 6/18/09
									shadow_vert_dists[1] = points_dist_3d(sun_pos, svert2); // added 6/18/09
									shadow_vert_dists[2] = points_dist_3d(sun_pos, svert3); // added 6/18/09

									if ((shadow_vert_dists[0]<shadow_vert_dists[1]) and
										(shadow_vert_dists[0]<shadow_vert_dists[2])) {
										shadow_least_dist=shadow_vert_dists[0];
										shadow_which_least=0;
									}
									if ((shadow_vert_dists[1]<shadow_vert_dists[2]) and
										(shadow_vert_dists[1]<shadow_vert_dists[0])) {
										shadow_least_dist=shadow_vert_dists[1];
										shadow_which_least=1;
									}
									if ((shadow_vert_dists[2]<shadow_vert_dists[0]) and
										(shadow_vert_dists[2]<shadow_vert_dists[1])) {
										shadow_least_dist=shadow_vert_dists[2];
										shadow_which_least=2;
									}

									if ( (shadow_least_dist) <= (tsample_find_max) ) {

										// need a second point on line formed by gridpoint in direction 'univ_dir'
										// 6/27/09: this will replace 'sun_pos' in geometry check below
										double newline_pt[3];

										for (int rep=0; rep<max_rayfires; rep++ )
										{
											ray_in_tri = 0;
											if (tsample_blocked[rep]==0) // do not check a ray already blocked
											{
												double temp_sample[3];
												temp_sample[0]=tsample_points[rep*3+0];
												temp_sample[1]=tsample_points[rep*3+1];
												temp_sample[2]=tsample_points[rep*3+2];

												newline_pt[0]=(temp_sample[0] + univ_dir[0]*10);
												newline_pt[1]=(temp_sample[1] + univ_dir[1]*10);
												newline_pt[2]=(temp_sample[2] + univ_dir[2]*10);
												/*triangle_contains_line_exp_3d(tri_s_vt, sun_pos, temp_sample, &ray_in_tri, ray_tri_pt, curr_s_triangle->Normal);*/
												triangle_contains_line_exp_3d(tri_s_vt, newline_pt, temp_sample, &ray_in_tri, ray_tri_pt, curr_s_triangle->Normal);
												raytri_dist=points_dist_3d(sun_pos, ray_tri_pt);
		//ATTENTION: raytri_dist<=tsample_sundist[rep]!!!

		//										double delta[3] = {ray_tri_pt[0]-temp_sample[0], ray_tri_pt[1]-temp_sample[1], ray_tri_pt[2]-temp_sample[2]};
		//										double dot = delta[0]*univ_dir[0] + delta[1]*univ_dir[1] + delta[2]*univ_dir[2];
		//
		//										if(ray_in_tri==1 and dot > 0) {
		//											num_points_shadowed+=1;
		//											tsample_blocked[rep]=1;
		//										}

												if ((ray_in_tri==1) and (raytri_dist<tsample_sundist[rep])) // true, collision, increment shadowed points. THIS raytri<tsample... IS THE MOST IMPORTANT DISTANCE CHECK!!!
												{
													num_points_shadowed+=1;
													tsample_blocked[rep]=1;
												}
											}
										}
									}
								} // end "if shadowing closer to sun"

		//						if(ii == 3) {
		//							double area;
		//							if (ray_method==PERCELL) {
		//								area = (max_rayfires-num_points_shadowed)*flatArea/(max_rayfires);
		//							}
		//							else {
		//								area = (max_rayfires-num_points_shadowed)*flatArea/(max_rayfires);
		//							}
		//							incidenceAngle = acos(-curr_triangle->SunScalarProduct);
		//							double reflPowerFlux=0;
		//							temp_power = CalculateTransmittedPower(area, incidenceAngle, &reflPowerFlux, nRatio);
		//
		//							//after the kk-th shadowing cell
		////							std::cout << kk << ": " << temp_power << std::endl;
		//						}

							} // end: loop through shadowing cells for given shadowed cell	*/







							//printf("points: %i\n", num_points_shadowed);
							// now calculate the power for this cell face using area, angle of incidence, and efficiency
							if (ray_method==PERCELL) {
								flatArea = (max_rayfires-num_points_shadowed)*flatArea/(max_rayfires);
							}
							else {
								flatArea = (max_rayfires-num_points_shadowed)*flatArea/(max_rayfires);
							}
							incidenceAngle = acos(-curr_triangle->SunScalarProduct);

							double reflPowerFlux=0;
							curCellPower = CalculateTransmittedPower(flatArea, incidenceAngle, &reflPowerFlux, efficiency, nRatio); // added 3/16

		//					std::cout << ii << ": " << temp_power << std::endl;

		//					if(ii == 3) std::cout << temp_power << std::endl;

							if (max_rayfires-num_points_shadowed !=0) {
								reflPowerFlux /= (max_rayfires-num_points_shadowed);}
							structPower += curCellPower; // previously 3/16, no intermed. variable called "temp_power"








							// SINGLE REFLECTION CODE started on 5/27/09////////
							double L_vec[3]; // vec from gridpt to lightsource (incident w/ opp. direction)
							double temp_magnitude=1;
							double temp_dot=0;
							double reflRay[3];
							double AaMm=0.005; // arbitrary multiplier to get second point on the line along ray direction
							double pointOnRay0[3];
							double pointOnRay1[3];
							bool ReflRay_in_tri = 0;
							double ReflRay_tri_pt[3];
							curCellReflPower = 0; 		// reset

							if ( (doSingleReflections==1) and (num_points_shadowed < max_rayfires) )
							{
							for (int ngp=0; ngp<max_rayfires; ngp++) {
								if (tsample_blocked[ngp]==0) { // for any gridpoint subtriangle that is not shadowed

									// first, find the direction of reflected ray from incident ray
									/*L_vec[0] = sun_pos[0] - tsample_points[ngp*3+0];
									L_vec[1] = sun_pos[1] - tsample_points[ngp*3+1];
									L_vec[2] = sun_pos[2] - tsample_points[ngp*3+2];*/
									L_vec[0] = -univ_dir[0]; // changed 6/27/09
									L_vec[1] = -univ_dir[1];
									L_vec[2] = -univ_dir[2];
									temp_magnitude = sqrt(L_vec[0]*L_vec[0] +
												L_vec[1]*L_vec[1] +
												L_vec[2]*L_vec[2]);
									if (temp_magnitude == 0) {temp_magnitude=1;}
									L_vec[0] /=temp_magnitude;
									L_vec[1] /=temp_magnitude;
									L_vec[2] /=temp_magnitude;

									// R = 2(N dot L)N - L. assumption that L and N are unit normalized
									temp_dot = r8vec_dot(3, L_vec, curr_triangle->Normal);
									reflRay[0] = 2*temp_dot*(curr_triangle->Normal[0]) - L_vec[0];
									reflRay[1] = 2*temp_dot*(curr_triangle->Normal[1]) - L_vec[1];
									reflRay[2] = 2*temp_dot*(curr_triangle->Normal[2]) - L_vec[2];
									pointOnRay0[0] = (tsample_points[ngp*3+0]);
									pointOnRay0[1] = (tsample_points[ngp*3+1]);
									pointOnRay0[2] = (tsample_points[ngp*3+2]);
									pointOnRay1[0] = (tsample_points[ngp*3+0] + AaMm*reflRay[0]);
									pointOnRay1[1] = (tsample_points[ngp*3+1] + AaMm*reflRay[1]);
									pointOnRay1[2] = (tsample_points[ngp*3+2] + AaMm*reflRay[2]);

									// second, check if this line intersects any other triangles...
									// of those it does, figure out which intersect the RAY not just the LINE
									// and then which one is closest to gridpoint because first hit is only hit
									Cell_3 *reflget_cell;
									int intersected_counter=0;		// count # of maybe "hit" triangles
									int intersected_ones[numCells]; // store indices of triangles maybe "hit"
									double int_dists0[numCells]; // distances intersect pt's to refl gridpt.
									double int_dists1[numCells]; // distances intersect pt's to point on RAY





									for (int rr=0; rr<numCells; rr++)
									{
										ReflRay_in_tri = 0; // reset, in case
										if (rr!=curCell) { // don't look at cell which reflected the ray...
											reflget_cell = &cells_Array[rr];
											// next line just uses same variable as shadowing method above
											curr_s_triangle = reflget_cell->triangleList[0];
											// shouldn't matter which triangle is used? 0 or 1

											// these variables defined above shadowed cell loop
											fpt_ptr = reflget_cell->vertices[curr_s_triangle->index[0]];
												tri_s_vt[0]=fpt_ptr->x;
												tri_s_vt[1]=fpt_ptr->y;
												tri_s_vt[2]=fpt_ptr->z;
											fpt_ptr = reflget_cell->vertices[curr_s_triangle->index[1]];
												tri_s_vt[3]=fpt_ptr->x;
												tri_s_vt[4]=fpt_ptr->y;
												tri_s_vt[5]=fpt_ptr->z;
											fpt_ptr = reflget_cell->vertices[curr_s_triangle->index[2]];
												tri_s_vt[6]=fpt_ptr->x;
												tri_s_vt[7]=fpt_ptr->y;
												tri_s_vt[8]=fpt_ptr->z;
											// condition #1: line described by reflected vector intersects tri.
											triangle_contains_line_exp_3d(tri_s_vt, pointOnRay0, pointOnRay1, &ReflRay_in_tri, ReflRay_tri_pt, curr_s_triangle->Normal);
											if (ReflRay_in_tri==1) // intersection found
											{
												double delta[3] = {ReflRay_tri_pt[0]-pointOnRay0[0], ReflRay_tri_pt[1]-pointOnRay0[1], ReflRay_tri_pt[2]-pointOnRay0[2]};
												double dot = delta[0]*reflRay[0] + delta[1]*reflRay[1] + delta[2]*reflRay[2];


												int_dists0[rr]=points_dist_3d(pointOnRay0, ReflRay_tri_pt);
		//										int_dists1[rr]=points_dist_3d(pointOnRay1, ReflRay_tri_pt);
												// condition #2: the tri. must be intersected by the RAY
												// not just any part of the LINE
												// therefore its distance to reflecting gridpoint must be
												// greater than distance to pointonray1 (+ along vector)

												if(dot > 0) {
													intersected_ones[intersected_counter]=rr;
													intersected_counter++;
												}

		//										if (int_dists0[rr]>int_dists1[rr]) {
		//											intersected_ones[intersected_counter]=rr;
		//											intersected_counter++;
		//										}
											}
										}
									} // end loop through possible cells to hit

									// condition #3: of all the possible intersections with ray, only one
									// actually gets hit. this one is that of closest distance to reflecting gridpoint
									// this doesn't require a second pass necessarily... but here it is for now
									double best_dist_sofar=999999999;
									int best_tri_sofar=-1;
									for (int intc=0; intc<intersected_counter; intc++)
									{
										if (int_dists0[intersected_ones[intc]]<best_dist_sofar) {
											best_dist_sofar = int_dists0[intersected_ones[intc]];
											best_tri_sofar = intersected_ones[intc];
										}
									}

									// third, calculate incidence angle of this ray with the HIT triangle.
									// multiply cosine of this by the reflPower variable
									// then add this to the totalInstantPower counter variable
									//double refGet_inc_angle=0;
									double final_dot=0;
									//int rti_factor = max_rayfires - num_points_shadowed;
									//double rtd_factor = (double)rti_factor;
									if (best_tri_sofar >=0) // there actually exists triangle hit with refl'd ray
									{
										reflget_cell = &cells_Array[best_tri_sofar]; // set the hit triangle
										curr_s_triangle = reflget_cell->triangleList[0]; //positive side
										final_dot = r8vec_dot(3, reflRay, curr_s_triangle->Normal);

										bool reflection = true;

										//proceed only if the positive side is facing the light
										if(globalCellsOneSided && final_dot >= 0) {
											reflection = false;
										}

										if(reflection) {
											if (final_dot<0) {// "wrong" triangle side normal
												curr_s_triangle = reflget_cell->triangleList[1];
												final_dot = r8vec_dot(3, reflRay, curr_s_triangle->Normal);
											}
											//refGet_inc_angle=acos( final_dot );

											// sum of all power absorbed by ANY cell from reflected rays
											//tmpref_power += (final_dot * reflPower)/(rtd_factor);
											// MARCO
											double inc_angle = acos(final_dot);
											double transm_angle = asin(nRatio*sin(inc_angle));	// snell's law...
											// ///////// Jun 16th 2011 - load best-hit cell (best_tri_sofar) properties
				 	                                                efficiency = efficiencies[best_tri_sofar];
                                                                                        nRatio = index_ref_air/reflIndices[best_tri_sofar];
	
											// Implement Fresnel Equations as in Wikipedia, using both incident and transmitted angle		
	                                                                                double RP_check = (nRatio*cos(transm_angle)-cos(inc_angle))/(nRatio*cos(transm_angle)+cos(inc_angle));
                                                                                        RP_check = RP_check*RP_check;
                                                                                        double RS = (nRatio*cos(inc_angle)-cos(transm_angle))/(nRatio*cos(inc_angle)+cos(transm_angle));
                                                                                        RS = RS*RS;
                                                                                        double Rcoeff = (RP_check+RS)/2; // take the average for unpolarized light
                                                                                        double Tcoeff_ref = 1-Rcoeff;// new Tcoefficient

											///////////////// Jun 16th 2011 - Eliminate cosine factor (final_dot) since best-hit cell hit completely
										        // by the projection of the reflecting cell's gridpoint for a fine enough grid, with an error of 1/sqrt(N), sqrt(N)=gridpoints in input file

											curCellReflPower+=reflPowerFlux*efficiency*Tcoeff_ref/*finaldot(curr_s_triangle->Area)/max_rayfires*/;
											double REFL_NRG=(3600*delta_time)*reflPowerFlux*efficiency*final_dot*Tcoeff_ref/max_rayfires;
											if (getSubCellData==true) {
												for (int otherc=0; otherc<max_rayfires; otherc++) {
													globalDataSubCellPower_OneGuy[which_iteration*numCells*max_rayfires+best_tri_sofar*max_rayfires+otherc]+=REFL_NRG;
													globalDataSubCellPower_Reflect[which_iteration*numCells*max_rayfires+best_tri_sofar*max_rayfires+otherc]+=REFL_NRG;
												}
											}
										}
									}
								}
							} // end for loop through gridpoints

		//					std::cout << ii << ": " << tmpref_power << std::endl;

							structPower += curCellReflPower;
							totalInstantPower_fromRefl += curCellReflPower;
							} // endif of whether to do reflections

							///END OF REFLECTION CODE /////////////////////////////////////////////////////////////





							//// THIS STUFF IS USED IN "GRID" verson of program run to save subcell power to file
							// so total power cell got in instant was temp_power +

							if ( /*(GA_gen_current%1000==0) ||*/ (getSubCellData==true) )
							{
								int ptt_factor = max_rayfires - num_points_shadowed;
								double ptd_factor= (double)ptt_factor;
								double scaled_power = (curCellPower)*delta_time*3600; //actually not power at all now
								if (globalUseKWH==true) {scaled_power/=Ws_to_kWhr;}
								if (ptd_factor!=0) {
									scaled_power /= ptd_factor;
									//printf("%f\n", scaled_power);
								}
								for (int subb=0; subb < max_rayfires; subb++) // save power for each subguy
								{
									if (tsample_blocked[subb]==0) {
										if (ff==0) {
											globalDataSubCellPower_OneGuy[which_iteration*numCells*max_rayfires+curCell*max_rayfires+subb]+=scaled_power;
										}
										else
										{
											SECONDglobalDataSubCellPower_OneGuy[which_iteration*numCells*max_rayfires+curCell*max_rayfires+subb]+=scaled_power;
										}
							// 6/25/09 changed '=scaledpower' to '+=scaledpower' because above there is reflection part added!!
									}
									else {
										if (ff==0) {
											globalDataSubCellPower_OneGuy[which_iteration*numCells*max_rayfires +curCell*max_rayfires+subb]+=0;
										}
										else
										{
											SECONDglobalDataSubCellPower_OneGuy[which_iteration*numCells*max_rayfires +curCell*max_rayfires+subb]+=0;
										}
									}
								}
							}
							////// END OF CODE that save subcell power data ////////////////

						} // end: if triangle faces sun block  (no pun intended)
					} // end: loop through victim's triangle faces

					cellPowerFile << curCell << "\t\t" << (curCellPower + curCellReflPower) << endl;

				} // end: shadow victim cell loop



				structPowerFile << curTime << "\t\t" << structPower << endl;

				cellPowerFile.close();

				//totalInstantPower is the power gather within one time slot
				DayPowerValues[curTime] = structPower;    // save value of dP (infinitesimal power) with key = current time tt

		//		std::cout << tt << ": " << totalInstantPower << std::endl;
		//		std::cout << totalInstantPower << std::endl;

				totalDayEnergy += structPower*delta_time*3600; // delta_time in fractional hours times #seconds in 1 whole hour
				totalDayEnergyFromRefl += totalInstantPower_fromRefl*delta_time*3600; // just to keep track of reflection contribution
				which_iteration++; // next iteration for subcell array
		//		printf("average number for all cells, checked shadowings: %i\n", count_checkers/numCells);
		//		Write_tsampledata=false;
		//		printf("test: %i: %f\n", which_iteration, totalInstantPower);
				Write_vertdata=false;
				Write_tsampledata=false;	// only write the coordinate data for one iteration in day
				delete tempsun;
			}


			structPowerFile.close();


			///////// store new individual in old slot ///////////
			qwikStoreEnergies[qwikIndex] = totalDayEnergy;
			qwikStoreReflEnergies[qwikIndex] = totalDayEnergyFromRefl;
			for (int qi1=0; qi1<numCells*9; qi1++) {
				qwikStoreCoords[qwikIndex][qi1] = xgenes[qi1];
			}
			qwikIndex++; //update for the next time around
			if (qwikIndex==qwikIndMax) {qwikIndex=0;}
			///////////////////////////////////////////////////////

			// test
			/*for (int rep1=0; rep1<500001; rep1++) {
				if (globalDataSubCellPower_Reflect[rep1]!=0)
				{
					printf("%f\n", globalDataSubCellPower_Reflect[rep1]);
				}
			}*/

			globalNumDayDivisions=num_times; // iterations.
			Write_tsampledata=false;

			totalEnergy += totalDayEnergy;///totalSurfaceArea;   // return final value for fitness of individual
			totalEnergyFromRefl += totalDayEnergyFromRefl;

			nextDay(curYear, curMonth, curDate);


		}

	}


	// delete
	for (int ddd=0; ddd<numCells; ddd++)
	{
		Cell_3 *curr_cell = &cells_Array[ddd];
		/*for (int dds=0; dds<curr_cell->numEdges; dds++) {
			delete curr_cell->edgeList[dds]->shadowPlane;
		}*/
//		for (int dde=0; dde<curr_cell->numEdges; dde++) {
//			delete curr_cell->edgeList[dde];
//		}
		for (int ddt=0; ddt<curr_cell->numTriangles; ddt++) {
			delete curr_cell->triangleList[ddt];
		}
		for (int ddv=0; ddv<curr_cell->numVertices; ddv++) {
			delete curr_cell->vertices[ddv];
		}
	}

	Done_sundata=true;  // set true whether or not it was true/false
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;

	if (globalUseKWH==true)
		{printf("time = %f sec, energy = %f kWh\n", time, totalEnergy);}//SGA OUTPUT
	else
		{printf("time = %f sec, energy = %f WattSec : %f kWh, %f WattSec : %f KWh regained from refl\n", time, totalEnergy, totalEnergy/(3600000), totalEnergyFromRefl, totalEnergyFromRefl/(3600000));}//SGA OUTPUT


	return totalEnergy;

}

// stuff at end of optimization
/////############################33
//////// OUTPUT                 ###
//################################
void SaveGenerationOutput () //NOT BEING USED
{
	
}

void DisplayFinalOutput ()
{
	// first output a file with just the best fitness value in each generation (for gnuplot??)
	printf("%s", Unique_Dir);
	strcpy(Clone_Dir, Unique_Dir);	
	char *tmp_str1 = strcat(Clone_Dir, "outputBestFit");
	FILE *outFileObj = fopen (tmp_str1/*"/outputBestFit"*//*"/home/brosbam/Desktop/solar_results/optimization/outputBestFit"*/, "w");
	int gen_factor = numCells*9/*6*/; // ie 0 + all decision vars +1
	for (int jj=0; jj<GA_maxgenerations; jj++) {
		if (GA_numberof_objectives==1) {
			fprintf(outFileObj, "%i %f\n", jj, globalObj0AvgBest[jj]/*globalDataFitness[jj]*/);
		}
		else 
		{
			fprintf(outFileObj, "%i %f %f\n", jj, globalObj0AvgBest[jj], globalObj1AvgBest[jj]);
		}

	}
	fclose (outFileObj);


	// output a single individual from all the generations (not always the best indiv,.. ehhhh)
	strcpy(Clone_Dir, Unique_Dir);
	char *tmp_str0 = strcat(Clone_Dir, OutBest_Indivs);
	FILE *outFileAll = fopen (tmp_str0,"w");
	//int gen_factor = numCells*9;
	for (int jj=0; jj<GA_maxgenerations; jj++) {
		if (GA_numberof_objectives==1) {
			for (int wn=0; wn<numCells; wn++) // added FOR SGA
			{
				for (int wdv=0; wdv<9; wdv++) {
					fprintf(outFileAll, "%f\t", globalParamsAvgBest[jj*gen_factor+wn*9+wdv]);
				}
			}
			fprintf(outFileAll, "%f\n", globalObj0AvgBest[jj]);
		}
		else 
		{
			for (int wn=0; wn<numCells; wn++)
			{
				for (int wdv=0; wdv<9; wdv++) {
					fprintf(outFileAll, "%f\t", globalParamsAvgBest[jj*gen_factor+wn*9+wdv]);
				}
			}
			fprintf(outFileAll, "%f\t%f\n", globalObj0AvgBest[jj], globalObj1AvgBest[jj]);
		}
	}
	fclose (outFileAll);
	
	// second, recreate best indiv. to find vertices, then output vertices as I did before using GA (for povray)
//	Cell_3::cell_counter = 0;             // initialize class static 'counter' variable
	


/*preliminary setup for the best structure of the last generation*/
//	Cell_3 cellOut[numCells];
//	Cell_3 *cp;
//
//	int l_gen=GA_maxgenerations-1;
//	WriteTextModelFile(l_gen); // ADD FOR SGA WHEN OPENGL CODE REMOVED: saves this file for viewing in mm3d
//
//	// model pointer to first position of first decision variable of first cell in best individual of this generation
//	modelPointer = &globalParamsAvgBest[l_gen*gen_factor];
//	totalSurfaceArea=0;
//	BuildStructure(cellOut, MODEL);

	


/*output outputBestPOV*/
//	strcpy(Clone_Dir, Unique_Dir);
//	char *tmp_str2 = strcat(Clone_Dir, "outputBestPOV");
//	FILE *outFilePOV = fopen (tmp_str2 /*"results/outputBestPOV"*//*"/home/brosbam/Desktop/solar_results/optimization/outputBestPOV.pov"*/, "w");
//	fprintf(outFilePOV, "#include \"colors.inc\"\n  background { color  Cyan }\n global_settings { ambient_light rgb<1, 1, 1> }\n camera {location <60, 40, 55> look_at  <12,20, 14>}\n");
//	fprintf(outFilePOV, "plane { y, -0.1\n pigment { checker Green White }}\n");
//	for (int pp=0; pp<numCells; ++pp)
//	{
//		cp=&cellOut[pp];
//		for (int tp=0; tp < cp->numTriangles; ++tp)
//		{
//			Triangle *tri = cp->triangleList[tp];
//
//		fprintf(outFilePOV, "triangle {< %f,%f,%f>, <%f,%f,%f>, <%f,%f,%f>\n pigment{ Red }}\n",
//				tri->vertList[1],
//				tri->vertList[2],
//				tri->vertList[0],
//				tri->vertList[4],
//				tri->vertList[5],
//				tri->vertList[3],
//				tri->vertList[7],
//				tri->vertList[8],
//				tri->vertList[6]);
//		}
//	}
//	fprintf(outFilePOV, "light_source { <-200, 400, 152> color White}\n");
//	fclose(outFilePOV);






	// third, since I have this individual might as well output its instantaneous power vs time as well, for gnuplot
	//CalculatePowerDayCycle(cellOut);  // getting DaypowerValues map for best individual
	
/*output outputBestPower*/
//	newCalculatePowerDayCycle(cellOut, modelPointer); // second argument added 08/04/09
//
//	strcpy(Clone_Dir, Unique_Dir);
//	char *tmp_str3 = strcat(Clone_Dir, OutBest_Power);
//	FILE *outFilePow = fopen (tmp_str3/*"/results/outputBestPower"*//*/home/brosbam/Desktop/solar_results/optimization/outputBestPower"*/, "w");
//	for (std::map<double, double>::iterator jj=DayPowerValues.begin(); jj!=DayPowerValues.end(); ++jj)
//	{
//		fprintf(outFilePow, "%f %f \n", jj->first, jj->second);
//	}
//	fclose(outFilePow);

	
	
	//system("/usr/bin/gnuplot"); // COMMENTED OUT FOR MULTIPLE RUNS AT SAME TIME
	//system("plot \"/home/brosbam/workspace_solar_geom/output.txt\"");
}

//######################################################
// SAVE POWER PLOT FOR ANY GENERATIONEST INDIVIDUAL!!!!!!!!!!!!!!!!!!!!
//#####################################################
// this first one IS OBSOLETE: USED AT THE TIME WHEN OPENGL WAS being used for graphics!1!
//JIN: NOT USED RIGHT NOW
//void DisplayPowerOutput (int l_gen)
//{
//	int gen_factor = numCells*9/*6*/; // ie 0 + all decision vars +1
//
//	// second, recreate best indiv. to find vertices, then output vertices as I did before using GA (for povray)
////	Cell_3::cell_counter = 0;             // initialize class static 'counter' variable
//	Cell_3 cellOut[numCells];
//
//	//int l_gen=GA_maxgenerations-1;
//	// model pointer to first position of first decision variable of first cell in best individual of this generation
//	modelPointer = &globalParamsAvgBest[l_gen*gen_factor];
//	totalSurfaceArea=0;
//	BuildStructure(cellOut, MODEL);
//
//	//CalculatePowerDayCycle(cellOut);  // getting DaypowerValues map for best individual
//	newCalculatePowerDayCycle(cellOut, modelPointer); // 2nd argument added on 08/04/09
//
//	strcpy(Clone_Dir, Unique_Dir);
//	char *tmp_str3 = strcat(Clone_Dir, OutBest_Power);
//	strcat(tmp_str3, "_grid");
//
//	char istring[200];
//	sprintf(istring,"%i", raysPerCell);
//	strcat(tmp_str3, istring);
//
//	FILE *outFilePow = fopen (tmp_str3, "w");
//	for (map<double, double>::iterator jj=DayPowerValues.begin(); jj!=DayPowerValues.end(); ++jj)
//	{
//		fprintf(outFilePow, "%f %f \n", jj->first, jj->second);
//	}
//	fclose(outFilePow);
//}

void DisplayPowerOutputGridVersion() // use the power values from last calculated power in "grid" version of solar3d
{				// and also append the # of gridpoints to the end of filename
	strcpy(Clone_Dir, Unique_Dir);
	char *tmp_str3 = strcat(Clone_Dir, OutBest_Power);
	strcat(tmp_str3, "_grid");
	
	char istring[200];
	sprintf(istring,"%i", raysPerCell);	
	strcat(tmp_str3, istring);	

	FILE *outFilePow = fopen (tmp_str3, "w");
	for (map<double, double>::iterator jj=DayPowerValues.begin(); jj!=DayPowerValues.end(); ++jj)
	{
		fprintf(outFilePow, "%f %f \n", jj->first, jj->second);
	}
	fclose(outFilePow);
}

// added 6-4-09 for use in "grid"
void GetZenithAngleDistribution(double myCoordinates[])
{
	// reset
	for (int dff=0; dff<HorizBins_numberof; dff++) {globalCountDistrHAangleBins[dff]=0;}
	for (int dff=0; dff<HorizBins_numberof; dff++) {globalCountDistrHAangleReflBins[dff]=0;}

	//int thisgen = GA_maxgenerations-1;
	double h_angles[numCells];	// horizontal angles
	double SECONDh_angles[numCells];
	double vref[3] = {1, 0, 0};		// horizontal vector against which to measure angle
	double cellnorm[3] = {0, 1, 0};		// calculated for each cell	
	
	Cell_3 cells_Array[numCells];
	Cell_3 *current_cellh;
	modelPointer = myCoordinates; 		// set the pointer for SetupCell... function
	BuildStructure(cells_Array, MODEL);
	// ANGLE FROM HORIZONTAL!!!
	for (int hhh=0; hhh<numCells; hhh++)
	{		
		current_cellh = &cells_Array[hhh];
		// get surface normal for the current cell	
		vref[0] = 0;
		vref[1] = 0;
		vref[2] = 1;
		cellnorm[0]=current_cellh->triangleList[0]->Normal[0];
		cellnorm[1]=current_cellh->triangleList[0]->Normal[1];
		cellnorm[2]=current_cellh->triangleList[0]->Normal[2];

		h_angles[hhh] = radtodeg( acos(cellnorm[0]*vref[0] + cellnorm[1]*vref[1] + cellnorm[2]*vref[2]) );
		/*if (h_angles[hhh] > 90) {
			h_angles[hhh]=180 - h_angles[hhh]; // since double-sided => then 45 degrees equivalent to 135 degrees..etc
		}*/
		
		//printf("%i %f\n", hhh, h_angles[hhh]);
		SECONDh_angles[hhh] = 180-h_angles[hhh]; //flipside....
		//if (SECONDh_angles[hhh] > 360) {SECONDh_angles[hhh]-=360;}
	}

	// for now, create an approximated discrete distribution (round all angles to nearest degree)
	int ang_sort=0;
	double curr_ang=0.0;
	double deg_per_binh = (double)(180/HorizBins_numberof);
	double temp_weighted_inc = 1;
	
	// new on 10/06/09, 10/6/09
	// print out a list of just the unbinned angles and energies for each cell (one cell per line)
	// then in Matlab, it will be easier to choose a good number of bins and experiment with the data, print plot
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat1 = strcat(Clone_Dir, "UnbinnedZenithAnglesList_grid");
	//strcat(output_cat, itoa(whichgen));
	char istring[200];
	sprintf(istring,"%i", raysPerCell);	
	strcat(output_cat1, istring);

	FILE *OUTfile = fopen(output_cat1, "w");
	
	fprintf(OUTfile, "Angle\t"); fprintf(OUTfile, "Contribution\n");
	for (int hah=0; hah<numCells; hah++)
	{	
		// do for both triangles on cell
		fprintf(OUTfile, "%f\t%f\n", h_angles[hah], globalDataSubCellPower_Sum[hah]);
		fprintf(OUTfile, "%f\t%f\n", SECONDh_angles[hah], SECONDglobalDataSubCellPower_Sum[hah]);
	}

	fclose(OUTfile);
	//////////////////////////////////

	// sort the horizontal and polar angles into bins of certain size before saving to file in next function
	for (int hhh=0; hhh<numCells; hhh++) {
		curr_ang = h_angles[hhh];
		if (curr_ang==180) {
			// include it in the highest bin (since below calculation would give bin=15 which is out of bounds)
			ang_sort=HorizBins_numberof-1; // note, this is not autmomatic, need to change if i change # of bins 
		}
		else
			{ang_sort = (int)(floor(curr_ang/deg_per_binh));}
		temp_weighted_inc=globalDataSubCellPower_Sum[hhh];
		globalCountDistrHAangleBins[ang_sort]+=temp_weighted_inc; // count this angle by tossing it in the appropriate bin with its 
									// appropriate weighted value(angle interval)

		temp_weighted_inc=globalDataSubCellPower_ReflSum[hhh];	// same as above but only counting reflected rays absorbed 
		globalCountDistrHAangleReflBins[ang_sort]+=temp_weighted_inc;

		// for other triangle
		curr_ang = SECONDh_angles[hhh];
		if (curr_ang==180) {
			// include it in the highest bin (since below calculation would give bin=15 which is out of bounds)
			ang_sort=HorizBins_numberof-1; // note, this is not autmomatic, need to change if i change # of bins 
		}
		else
			{ang_sort = (int)(floor(curr_ang/deg_per_binh));}
		temp_weighted_inc=SECONDglobalDataSubCellPower_Sum[hhh];
		globalCountDistrHAangleBins[ang_sort]+=temp_weighted_inc; // count this angle by tossing it in the appropriate bin with its 
									// appropriate weighted value(angle interval)
	}
}

void GetAzimuthAngleDistribution(double myCoordinates[])
{
	// reset
	for (int dff=0; dff<AngleBins_numberof; dff++) {globalCountDistrPAangleBins[dff]=0;}
	for (int dff=0; dff<AngleBins_numberof; dff++) {globalCountDistrPAangleReflBins[dff]=0;}

	//int thisgen = GA_maxgenerations-1;
	double p_angles[numCells];	// horizontal angles
	double SECONDp_angles[numCells];
	double vref[3] = {1, 0, 0};		// horizontal vector against which to measure angle
//	double v_up[3] = {0, 0, 1};		// know which way is up
	double cellnorm[3] = {0, 1, 0};		// calculated for each cell	
	
	Cell_3 cells_Array[numCells];
	Cell_3 *current_cellp;
	modelPointer = myCoordinates; 		// set the pointer for SetupCell... function
	BuildStructure(cells_Array, MODEL);
	// ANGLE FROM EAST?!!
	for (int hhh=0; hhh<numCells; hhh++)
	{		
		current_cellp = &cells_Array[hhh];
		// get surface normal for the current cell	
		cellnorm[0]=current_cellp->triangleList[0]->Normal[0];
		cellnorm[1]=current_cellp->triangleList[0]->Normal[1];
		cellnorm[2]=current_cellp->triangleList[0]->Normal[2];

		/*if (r8vec_dot(3, cellnorm, v_up) < 0) {
			// if this triangle faces down, switch to the other one
			// the up facing triangle most likely gets more sun anyway...besides this consistency
			cellnorm[0]=current_cellp->triangleList[1]->Normal[0];
			cellnorm[1]=current_cellp->triangleList[1]->Normal[1];
			cellnorm[2]=current_cellp->triangleList[1]->Normal[2];
		}*/

		p_angles[hhh] = radtodeg( acos(cellnorm[0]*vref[0] + cellnorm[1]*vref[1] + cellnorm[2]*vref[2]) );
		if (cellnorm[1]<0) {p_angles[hhh] += 2*(180-p_angles[hhh]);} // y-direction, reflect about x if y negative
		//printf("%i %f\n", hhh, p_angles[hhh]);
		if ((cellnorm[1]==0) and (cellnorm[0]==0))
		{
			// vertical facing => has no polar direction
			p_angles[hhh]=12001; // secret code
			SECONDp_angles[hhh]=12001;
		}
		else
		{
			SECONDp_angles[hhh] = p_angles[hhh]+180; //flipside....
			if (SECONDp_angles[hhh] > 360) {SECONDp_angles[hhh]-=360;}
		}
	}	

	// new on 10/06/09, 10/6/09
	// print out a list of just the unbinned angles and energies for each cell (one cell per line)
	// then in Matlab, it will be easier to choose a good number of bins and experiment with the data, print plot
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat1 = strcat(Clone_Dir, "UnbinnedPolarAnglesList_grid");
	//strcat(output_cat, itoa(whichgen));
	char istring[200];
	sprintf(istring,"%i", raysPerCell);	
	strcat(output_cat1, istring);

	FILE *OUTfile = fopen(output_cat1, "w");
	
	fprintf(OUTfile, "Angle\t"); fprintf(OUTfile, "Contribution\n");
	for (int hah=0; hah<numCells; hah++)
	{	
		// do for both triangles on cell
		if (p_angles[hah]!=12001) {
			fprintf(OUTfile, "%f\t%f\n", p_angles[hah], globalDataSubCellPower_Sum[hah]);
			fprintf(OUTfile, "%f\t%f\n", SECONDp_angles[hah], SECONDglobalDataSubCellPower_Sum[hah]);
		}
	}

	fclose(OUTfile);
	//////////////////////////////////

	// for now, create an approximated discrete distribution (round all angles to nearest degree)
	int ang_sort=0;
	double curr_ang=0.0;
	double deg_per_bin = (double)(360/AngleBins_numberof);
	double temp_weighted_inc = 1;
	double spec1_weighted_inc = 1;
	double spec2_weighted_inc = 1;
	
	// sort the horizontal and polar angles into bins of certain size before saving to file in next function
	for (int hhh=0; hhh<numCells; hhh++) {
		curr_ang = p_angles[hhh];		
		if (curr_ang<0) {
			curr_ang +=360; // want these angles to be characterized by their positive definition (eg, -75 -> 285)		
		}
		if (curr_ang==360) {curr_ang=0;}
		if (curr_ang==12001) {
			// if it has no polar angle, just split the energy between all angle bins
			spec1_weighted_inc = (globalDataSubCellPower_Sum[hhh]) / ((double)AngleBins_numberof);
			spec2_weighted_inc = (SECONDglobalDataSubCellPower_Sum[hhh]) / ((double)AngleBins_numberof);
			
			for (int a12 = 0; a12 < AngleBins_numberof; a12++) {
				globalCountDistrPAangleBins[a12] += spec1_weighted_inc + spec2_weighted_inc;
			}
		}
		else
		{
			ang_sort = (int)(floor(curr_ang/deg_per_bin));
			temp_weighted_inc=globalDataSubCellPower_Sum[hhh];
			globalCountDistrPAangleBins[ang_sort]+=temp_weighted_inc; // count this angle by tossing in the approp bin with its 
									// appropriate weighted value(angle interval)

			temp_weighted_inc=globalDataSubCellPower_ReflSum[hhh];
			globalCountDistrPAangleReflBins[ang_sort]+=temp_weighted_inc;

			// flipside, other triangle
			curr_ang = SECONDp_angles[hhh];		
			if (curr_ang<0) {
				curr_ang +=360; // want these angles to be charac'd by their (+) definition (eg, -75 -> 285)		
			}
			if (curr_ang==360) {curr_ang=0;}
			ang_sort = (int)(floor(curr_ang/deg_per_bin));
			temp_weighted_inc=SECONDglobalDataSubCellPower_Sum[hhh];
			globalCountDistrPAangleBins[ang_sort]+=temp_weighted_inc;
		}  					
	}
}

void Get2dAngleDistribution(double myCoordinates[])
{
	for (int xff=0; xff<HorizBins_numberof; xff++) {
		for (int yff=0; yff<AngleBins_numberof; yff++) {	
			global2dDistrBothAngleBins[xff][yff]=0; // reset all sums 
		}
	}

	//int thisgen = GA_maxgenerations-1;
	double h_angles[numCells];	// horizontal angles
	double p_angles[numCells];
	double SECONDh_angles[numCells];	// horizontal angles
	double SECONDp_angles[numCells];
	double vref[3] = {1, 0, 0};		// horizontal vector against which to measure angle
//	double v_up[3] = {0, 0, 1};		// know which way is up
	double cellnorm[3] = {0, 1, 0};		// calculated for each cell	
	
	Cell_3 cells_Array[numCells];
	Cell_3 *current_cellh;
	Cell_3 *current_cellp;
	modelPointer = myCoordinates; 		// set the pointer for SetupCell... function
	BuildStructure(cells_Array, MODEL);
	// ANGLE FROM HORIZONTAL!!!
	for (int hhh=0; hhh<numCells; hhh++)
	{		
		current_cellh = &cells_Array[hhh];
		// get surface normal for the current cell	
		vref[0] = 0;
		vref[1] = 0;
		vref[2] = 1;
		cellnorm[0]=current_cellh->triangleList[0]->Normal[0];
		cellnorm[1]=current_cellh->triangleList[0]->Normal[1];
		cellnorm[2]=current_cellh->triangleList[0]->Normal[2];

		h_angles[hhh] = radtodeg( acos(cellnorm[0]*vref[0] + cellnorm[1]*vref[1] + cellnorm[2]*vref[2]) );
		/*if (h_angles[hhh] > 90) {
			h_angles[hhh]=180 - h_angles[hhh]; // since double-sided => then 45 degrees equivalent to 135 degrees..etc
		}*/
		
		//printf("%i %f\n", hhh, h_angles[hhh]);
		SECONDh_angles[hhh] = 180-h_angles[hhh]; //flipside....
		//if (SECONDh_angles[hhh] > 360) {SECONDh_angles[hhh]-=360;}

		//if (h_angles[hhh]==90) {printf("Got a 90\n");}
	}

	// ANGLE FROM EAST?!!
	for (int hhh=0; hhh<numCells; hhh++)
	{		
		current_cellp = &cells_Array[hhh];
		// get surface normal for the current cell
		vref[0] = 1;
		vref[1] = 0;
		vref[2] = 0;	
		cellnorm[0]=current_cellp->triangleList[0]->Normal[0];
		cellnorm[1]=current_cellp->triangleList[0]->Normal[1];
		cellnorm[2]=current_cellp->triangleList[0]->Normal[2];

		p_angles[hhh] = radtodeg( acos(cellnorm[0]*vref[0] + cellnorm[1]*vref[1] + cellnorm[2]*vref[2]) );
		if (cellnorm[1]<0) {p_angles[hhh] += 2*(180-p_angles[hhh]);} // y-direction, reflect about x if y negative
		//printf("%i %f\n", hhh, p_angles[hhh]);
		if ((cellnorm[1]==0) and (cellnorm[0]==0))
		{
			// vertical facing => has no polar direction
			p_angles[hhh]=12001; // secret code
			SECONDp_angles[hhh]=12001;
		}
		else
		{
			SECONDp_angles[hhh] = p_angles[hhh]+180; //flipside....
			if (SECONDp_angles[hhh] > 360) {SECONDp_angles[hhh]-=360;}
		}
		if (p_angles[hhh]==90) {printf("Got a 90\n");}
	}
	
	int ang_sort_x=0;
	int ang_sort_y=0;
	double curr_ang_x=0.0;
	double curr_ang_y=0.0;
	double deg_per_bin_x = (double)(180/HorizBins_numberof);
	double deg_per_bin_y = (double)(360/AngleBins_numberof);
	double temp_weighted_inc = 1;
	double spec1_weighted_inc = 1;
	double spec2_weighted_inc = 1;

	for (int axy=0; axy<numCells; axy++) {
		curr_ang_x = h_angles[axy];
		curr_ang_y = p_angles[axy];
		
		// take three precautions to be sure all cells are counted
		if (curr_ang_x==180) {ang_sort_x = HorizBins_numberof-1;} // set to the max bin
		if (curr_ang_y<0) {
			curr_ang_y+=360; // want these angles to be characterized by their positive definition (eg, -75 ->285)		
		}
		if (curr_ang_y==360) {curr_ang_y=0;}		

		ang_sort_x = (int)(floor(curr_ang_x/deg_per_bin_x));
		if (curr_ang_y==12001) {
			// if it has no polar angle, just split the energy between all angle bins
			spec1_weighted_inc = (globalDataSubCellPower_Sum[axy]) / ((double)AngleBins_numberof);
			spec2_weighted_inc = (SECONDglobalDataSubCellPower_Sum[axy]) / ((double)AngleBins_numberof);
			
			for (int a12 = 0; a12 < AngleBins_numberof; a12++) {
				global2dDistrBothAngleBins[ang_sort_x][a12] += spec1_weighted_inc + spec2_weighted_inc;
			}
		}
		else
		{
			ang_sort_y = (int)(floor(curr_ang_y/deg_per_bin_y));
		
			temp_weighted_inc=globalDataSubCellPower_Sum[axy] + SECONDglobalDataSubCellPower_Sum[axy];

			global2dDistrBothAngleBins[ang_sort_x][ang_sort_y]+=temp_weighted_inc;
		}
	}

	strcpy(Clone_Dir, Unique_Dir);
	char *output_cat = strcat(Clone_Dir, OutSpherical_prefix);
	char jstring[200];
	sprintf(jstring,"%i", raysPerCell);	
	strcat(output_cat, jstring);
	FILE *OUTfile = fopen(output_cat, "w");
	for (int axy=0; axy<numCells; axy++) {
		// each line in file is a point in spherical polar coordinates
		fprintf(OUTfile, "%f\t", globalDataSubCellPower_Sum[axy] + SECONDglobalDataSubCellPower_Sum[axy]);	// r
		fprintf(OUTfile, "%f\t", p_angles[axy]);			// phi
		fprintf(OUTfile, "%f\n", h_angles[axy]);			// theta
	}
	fclose(OUTfile);
}

void WriteAngleDistrFile(/*int whichgen*/)
{
	// file 1
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat1 = strcat(Clone_Dir, OutpAngle_prefix);
	//strcat(output_cat, itoa(whichgen));
	char istring[200];
	sprintf(istring,"%i", raysPerCell);	
	strcat(output_cat1, istring);

	FILE *OUTfile = fopen(output_cat1, "w");
	
	fprintf(OUTfile, "Angles\t"); fprintf(OUTfile, "Contribution\n");
	for (int cff=0; cff<AngleBins_numberof; cff++)	
	{	
		fprintf(OUTfile, "%f\t%f\n", ((double)(cff*360))/((double)AngleBins_numberof), (globalCountDistrPAangleBins[cff])/(3600*1000) );
	}

	fclose(OUTfile);

	// NOW WRITE THE ONE FOR HORIZONTAL ANGLES
	// file 2
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat = strcat(Clone_Dir, OuthAngle_prefix);
	char jstring[200];
	sprintf(jstring,"%i", raysPerCell/*whichgen*/);	
	strcat(output_cat, jstring);

	FILE *OUTfile2 = fopen(output_cat, "w");

	fprintf(OUTfile2, "Angles\t"); fprintf(OUTfile2, "Contribution\n");
	for (int cff=0; cff<HorizBins_numberof; cff++)	
	{	
		//fprintf(OUTfile2, "%i\t%i\n", cff*90/HorizBins_numberof, globalCountDistrHAangleBins[cff]); // yes x values
		fprintf(OUTfile2, "%f\t%f\n", ((double)(cff*180))/((double)HorizBins_numberof), (globalCountDistrHAangleBins[cff])/(3600*1000) ); 
	}
	fclose(OUTfile2);

	// file 3
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat1r = strcat(Clone_Dir, OutpAngleRefl_prefix);
	//strcat(output_cat, itoa(whichgen));
	char irstring[200];
	sprintf(irstring,"%i", raysPerCell);	
	strcat(output_cat1r, irstring);

	FILE *OUTfiler = fopen(output_cat1r, "w");
	
	for (int cff=0; cff<AngleBins_numberof; cff++)	
	{	
		fprintf(OUTfiler, "%f\n",  (globalCountDistrPAangleReflBins[cff]) );
	}

	fclose(OUTfiler);

	// NOW WRITE THE ONE FOR HORIZONTAL ANGLES
	// file 4
	strcpy(Clone_Dir, Unique_Dir);

	char *output_catr = strcat(Clone_Dir, OuthAngleRefl_prefix);
	char jrstring[200];
	sprintf(jrstring,"%i", raysPerCell/*whichgen*/);	
	strcat(output_catr, jrstring);

	FILE *OUTfile2r = fopen(output_catr, "w");

	for (int cff=0; cff<HorizBins_numberof; cff++)	
	{	
		//fprintf(OUTfile2, "%i\t%i\n", cff*90/HorizBins_numberof, globalCountDistrHAangleBins[cff]); // yes x values
		fprintf(OUTfile2r, "%f\n", globalCountDistrHAangleReflBins[cff]); // no x values
	}
	fclose(OUTfile2r);

	// NOW WRITE THE ONE FOR BOTH ANGLES FOR 3D graph
	strcpy(Clone_Dir, Unique_Dir);

	char *output_cat2 = strcat(Clone_Dir, OutallAngle_prefix);
	char kstring[200];
	sprintf(kstring,"%i", raysPerCell);	
	strcat(output_cat2, kstring);

	FILE *OUTfile3 = fopen(output_cat2, "w");

	//commented out is the GRAPHIS-friendly format

	/*for (int xff=0; xff<HorizBins_numberof; xff++) {
		for (int yff=0; yff<AngleBins_numberof; yff++) {	
			fprintf(OUTfile3, "%f\t", (global2dDistrBothAngleBins[xff][yff])/(3600*1000));
		} // tab when incrementing on same line (x++), but newline for incrementing y++
		fprintf(OUTfile3, "\n");
	}*/

	// newer method for CoPlot
	fprintf(OUTfile3, "x\ty\tz\n");
	for (int xff=0; xff<HorizBins_numberof; xff++) {
		for (int yff=0; yff<AngleBins_numberof; yff++)  {
			if ( (global2dDistrBothAngleBins[xff][yff])/(3600*1000) !=0 ) {
				fprintf(OUTfile3, "%i\t%i\t%f\n", xff*180/HorizBins_numberof,yff*360/AngleBins_numberof, (global2dDistrBothAngleBins[xff][yff])/(3600*1000) );
			}
		}
	}
	fclose(OUTfile3);


	// added 8/21/09 
	// statistics
	
	// 1) (avg angle of power contrib) = sum(angle[i] x power[i])/sum(power[i])
	double hstat_sumAnglePower=0;
	double hstat_sumSqAnglePower=0;
	double hstat_sumPower=0;
	double hstat_avgAngle=0;
	double hstat_avgSqAngle=0;
	double hstat_stdDevAngle=0;
	for (int aii=0; aii<HorizBins_numberof; aii++) 
	{
		double the_s_angle = ((double)(aii*180))/((double)HorizBins_numberof);
		double the_s_power = (globalCountDistrHAangleBins[aii])/(3600*1000);
		hstat_sumAnglePower += the_s_angle * the_s_power;
		hstat_sumSqAnglePower+=the_s_angle*the_s_angle*the_s_power;
		hstat_sumPower += the_s_power; 
	}

	hstat_avgAngle = hstat_sumAnglePower / hstat_sumPower;
	hstat_avgSqAngle = hstat_sumSqAnglePower / hstat_sumPower;

	hstat_stdDevAngle = sqrt ( hstat_avgSqAngle - hstat_avgAngle*hstat_avgAngle );

	strcpy(Clone_Dir, Unique_Dir);
	char *stats_cat = strcat(Clone_Dir, "HAngle_stats_grid");
	char wstring[200];
	sprintf(wstring,"%i", raysPerCell);	
	strcat(stats_cat, wstring);

	FILE *OUTfileHstat = fopen(stats_cat, "w");

	fprintf(OUTfileHstat, "\n************Horizontal******************\n");
	fprintf(OUTfileHstat, "total Power\t%f\n", hstat_sumPower);
	fprintf(OUTfileHstat, "Average Angle\t%f\n", hstat_avgAngle);
	fprintf(OUTfileHstat, "StdDev. Angle\t%f\n", hstat_stdDevAngle);
	fprintf(OUTfileHstat, "****************************************\n");
	fclose(OUTfileHstat);

	
	double pstat_sumAnglePower=0;
	double pstat_sumSqAnglePower=0;
	double pstat_sumPower=0;
	double pstat_avgAngle=0;
	double pstat_avgSqAngle=0;
	double pstat_stdDevAngle=0;
	for (int aii=0; aii<AngleBins_numberof; aii++) 
	{
		double the_s_angle = ((double)(aii*360))/((double)AngleBins_numberof);
		double the_s_power = (globalCountDistrPAangleBins[aii])/(3600*1000);
		pstat_sumAnglePower += the_s_angle * the_s_power;
		pstat_sumSqAnglePower+=the_s_angle*the_s_angle*the_s_power;
		pstat_sumPower += the_s_power; 
	}

	pstat_avgAngle = pstat_sumAnglePower / pstat_sumPower;
	pstat_avgSqAngle = pstat_sumSqAnglePower / pstat_sumPower;

	pstat_stdDevAngle = sqrt ( pstat_avgSqAngle - pstat_avgAngle*pstat_avgAngle );

	strcpy(Clone_Dir, Unique_Dir);
	char *stats2_cat = strcat(Clone_Dir, "PAngle_stats_grid");
	char w2string[200];
	sprintf(w2string,"%i", raysPerCell);	
	strcat(stats2_cat, w2string);

	FILE *OUTfilePstat = fopen(stats2_cat, "w");

	fprintf(OUTfilePstat, "\n**************Polar****************\n");
	fprintf(OUTfilePstat, "total Power\t%f\n", pstat_sumPower);
	fprintf(OUTfilePstat, "Average Angle\t%f\n", pstat_avgAngle);
	fprintf(OUTfilePstat, "StdDev. Angle\t%f\n", pstat_stdDevAngle);
	fprintf(OUTfilePstat, "***********************************\n");
	fclose(OUTfilePstat);
}

//######################################################################
//######################################################################
// PREVIOUSLY IN USERDEFINABLES.CPP of GA CODE
//######################################################################
//######################################################################

// fitness function and constraints
void globalEvaluate(double *x, double *objArray, double *constraintViolation, double *penalty, int *noOfViolations) 
{	
	int ii;
	//FILE *outEvals; // COMMENTED OUT SINCE IM NO LONGER LETTING THE GA SAVE ITS OWN SOLUTIONS
	FILE *outMine;

	for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
		objArray[ii] = 0.0;
	}
	
	// OUR CUSTOM FITNESS FUNCTION !!!!
	if (globalSetup->gaType == SGA || globalSetup->gaType == NSGA  ) {	// 4/29/09 this is a useless 'if' statement....
		// single objective --> objArray[0]
		
		/* first time get coordinates from file unless random initial population
		 * but then the model coordinates will just be given by the x area passed
		 * into this globalEvaluate call */
//		Cell_3::cell_counter = 0;             // initialize class static 'counter' variable
		totalSurfaceArea=0;
		Cell_3 cellArray[globalSetup->noOfDecisionVariables/9];
	
		// commented out above blocky block for this easy one-line thing:
		modelPointer = x; // global model pointer simply points to this array "x" passed in
		BuildStructure(cellArray, MODEL/*TRIFLAT*/);                // store cells in cellArray. compose overall energy-collecting structure

		double efficiencies[numCells], reflIndices[numCells];

		for(int i=0; i<numCells; i++) {
			efficiencies[i] = cellEfficiency;
			reflIndices[i] = index_ref_plastic;
		}

		double day_energee=0;
		
		day_energee=newCalculatePowerDayCycle(cellArray, x, efficiencies, reflIndices);   // calculate power at several instants in the day

		globalMostRecentEnergyReturned = day_energee;

		// SGA verson of the above NSGA objective-saving code
		objArray[0] = day_energee;
		// added "or" part to next conditional statement on 4/29/09
		if ( ( (GA_gen_current%1000==0) && (day_energee  > globalBestPower_thisGen) ) || ( getSubCellData == true))
		{
			globalBestPower_thisGen=day_energee;
			globalBestGuy_thisGen=globalIndividualCounter; // current index
		}
		
		if (globalHaveDonePreEval==1) // SAVE THIS INDIVIDUAL if an 1000th generation, OR JUST INCREMENT COUNTER if not
		{
			int indN = globalNumDayDivisions*numCells*raysPerCell;
			int dayN = numCells*raysPerCell;
			int celN = raysPerCell;
			double temp_subpower=0;
			double SECONDtemp_subpower=0;
			double temp_sreflpower=0;
			
			for (int sreset=0; sreset < numCells; sreset++) {
				globalDataSubCellPower_Sum[sreset] = 0;
				SECONDglobalDataSubCellPower_Sum[sreset] = 0;
				globalDataSubCellPower_ReflSum[sreset]=0;
			}

			if ( /*(GA_gen_current%1000==0) ||*/ (getSubCellData == true) ) // right of "||" on 4/29/09 for "grid" execution
			{
				for (int stt=0; stt < globalNumDayDivisions; stt++)
				{
					for (int scc=0; scc<numCells; scc++)
					{
						for (int srr=0; srr< raysPerCell; srr++)
						{
							temp_subpower = globalDataSubCellPower_OneGuy[stt*dayN + scc*celN + srr];
							SECONDtemp_subpower = SECONDglobalDataSubCellPower_OneGuy[stt*dayN + scc*celN + srr];
							globalDataSubCellPower_Gen[globalIndividualCounter*indN + stt*dayN + scc*celN + srr]=temp_subpower;
							SECONDglobalDataSubCellPower_Gen[globalIndividualCounter*indN + stt*dayN + scc*celN + srr]=SECONDtemp_subpower;
							globalDataSubCellPower_Sum[scc] += temp_subpower;
							SECONDglobalDataSubCellPower_Sum[scc] += SECONDtemp_subpower;

							// added 6/25/09 for reflection analysis...hopefully useful
							temp_sreflpower = globalDataSubCellPower_Reflect[stt*dayN + scc*celN + srr];
							
							/*if (temp_sreflpower!=0 and stt<73) // test
							{
								printf("%f", temp_sreflpower); 
								printf("--->daypart: %i, cell#%i, ray#%i ", stt,scc,srr);
							}*/
							// 7/13/09, noticed I never finished this by adding temp... to the sum
							globalDataSubCellPower_ReflSum[scc] += temp_sreflpower;
						}
					}
				}
			}
			globalIndividualCounter++; // increment after each eval ON EVERY GENERATION EVERY EVAL!!!!
			
			// last individual has just been evaluated
			if ( /*( (globalIndividualCounter==GA_population) && (GA_gen_current%1000==0) ) ||*/ (getSubCellData==true) )
			{
				// write the subtriangle energies once every 1000 generations for the best individual

				//FILE *fileCPOW= fopen(output_catc, "a"); // moved 5/13/09 to down in first for loop
				for (int stt=0; stt < globalNumDayDivisions; stt++)
				{
					// added 5/13/09 to split up each iteration of subcell data over the day
					strcpy(Clone_Dir, Unique_Dir);

					char *output_catc = strcat(Clone_Dir, OutPowerPerSubcell);
					char scstring[200];

					sprintf(scstring,"%i",stt);	// 5/13/09 added
					strcat(output_catc, scstring); 	// 5/13/09 added: append iteration # to the filename
					FILE *fileCPOW= fopen(output_catc, "w");  // 5/13/09 switched this to "w" since one time, not "a"
					/////////////////////////////
					for (int scc=0; scc<numCells; scc++)
					{
						for (int srr=0; srr< raysPerCell; srr++)
						{
							fprintf(fileCPOW, "%f\n", globalDataSubCellPower_Gen[globalBestGuy_thisGen*indN + stt*dayN + scc*celN + srr] + SECONDglobalDataSubCellPower_Gen[globalBestGuy_thisGen*indN + stt*dayN + scc*celN + srr]);
						}
					}
					fclose(fileCPOW);
					//fprintf(fileCPOW, "\n");
				}
			}
		}
	}

	// note, here I deleted a bunch of penalty handling I can look up in userdefinables again if needed

	// added finally in March 9, 3/9/09, to get multiobjective data more reliably/easily/unambiguously...hopefully
	if ((globalHaveDonePreEval==1) /*&& (GA_gen_current%200 == 0)*/)  
	{
		strcpy(Clone_Dir, Unique_Dir);
		char *output_cat = strcat(Clone_Dir, OutAll_Indivs);	
		outMine = fopen(output_cat,"a");
		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			if (GA_gen_current%100==0) {fprintf(outMine, "%f\t", x[ii]);} // print only every 100th gen
			//globalDataParams[globalParamsCounter]=x[ii]; // commented out 3/17, replaced with below line
			globalDataParams[(globalIndividualCounter-1)*numCells*9+ii] = x[ii]; //counter>0 here since ++ above
			//globalParamsCounter++; // unused as of 3/17
		}
		for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
			if (GA_gen_current%100==0) {fprintf(outMine, "%f\t", objArray[ii]);}// print only every 100th gen
			if (ii==0)
				{globalDataObjective0[globalObjCounter] = objArray[0];}	// already distinguishes SGA or NSGA
			if (ii==1)
				{globalDataObjective1[globalObjCounter] = objArray[1];}	
		}
		globalObjCounter++;
		fprintf(outMine, "\n");
		fflush(outMine);
		fclose(outMine);
	}
	
	// ADDED 3/17
	// if last evaluation of 1 generation, then save the best guy to array, discard rest from memory (already written to disk)
	// and do some other stuff.....
	if (globalIndividualCounter==GA_population) {
		globalIndividualCounter=0; // reset
		globalBestPower_thisGen=0; // reset
		globalBestGuy_thisGen=0;  //reset	
		//printf("%i\n", GA_gen_current); // debugging

//JIN: for now
		StoreBestIndividual_Gen(GA_gen_current);

		// 5/31/09 : reset the random number seed every few generations just to make sure diff. runs don't get identical results
		if (GA_gen_current%50==0) {
			myRandom.setSeed(myRandom.makeSeed()); // from random.hpp and random.cpp stuff. reset seed 
		}
	
		// finally, for every 1000th generation print out a Model file after deleting old previous one
		strcpy(Clone_Dir, Unique_Dir);
		char *output_cmodel = strcat(Clone_Dir, OutModel_prefix);
		if (GA_gen_current%100==0) {
			if (GA_gen_current>100) { // shouldn't need to delete anything if this is not at least (2*period)th generation
				/*char pgenstring[200];
				sprintf(pgenstring,"%i", GA_gen_current-2); 
				strcat(output_cmodel, pgenstring);
				strcat(output_cmodel, "_cells.txt");*/
				strcat(output_cmodel, "*.txt");
				char delmodel_command[200];
				sprintf(delmodel_command, "mv ");
				strcat(delmodel_command, output_cmodel);
				strcat(delmodel_command, " ");
				strcpy(Clone_Dir, Unique_Dir);
				strcat(delmodel_command, Unique_Dir);
				strcat(delmodel_command, "models");
				//strcat(delmodel_command, output_cmodel);				
				printf(delmodel_command);
				system(delmodel_command);

				// could have more simply done "rm ....Model_gen*"
			}
			WriteTextModelFile(GA_gen_current); // write new model file after previous one is deleted
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

/*
  Read a non-comment, non-blank line from the specified file stream
  At EOF, it returns NULL.
  Otherwise, it returns the pointer to the first token.
  (The string just read in is altered by strtok().)
*/
static char* readOneLine(char *pcBuf, int iMaxSize, FILE *fStream) {
	
	char *pToken;

	do {
		pToken = NULL;

		*pcBuf = '\0';
		fgets(pcBuf, iMaxSize, fStream);
		if (feof(fStream))
			break;

		// get the first token
		pToken = strtok(pcBuf, BLANK_STR);

		// if there is no first token, it is a blank line.
		// if the first token starts with '#', it is a comment line.
	} while ((pToken == NULL) || (*pToken == '#'));

	return (pToken);
}

int mainGA(/*int argc, char *argv[]*/) {
	
	int ii;
	// the buffer size was what was killing me!!! screw this naive GA reader, finally found it
	const int ciBufSize = 50000;  // formerly 1024, but increased
	char *pToken, caBuf[ciBufSize];
	FILE *fInput, *fOutput;

	/*
	if (argc != 2) {
		printf("Error! Usage is GAtbx inputfile\n");
		exit(1);
	}

	fInput = fopen(argv[1], "r");
	if (fInput == NULL) {
		printf("Error! opening file %s\n", argv[1]);
		exit(1);
	}*/
	strcpy(Clone_Dir, Unique_Dir);
	fInput = fopen( strcat(Clone_Dir, Input_GA) , "r");
	if (fInput == NULL) {
		printf("Error! opening file %s\n", Input_GA);
		exit(1);
	}

	globalSetup = new GlobalSetup;

	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}

	///added 4/2/09 ////for area restriction
	globalSetup->do_area_restriction = allow_area_restrict;
	////////////////////////////////////////////

	//GA type
	if (strcmp("SGA", pToken) == 0) {
		globalSetup->gaType = SGA;
	} else if (strcmp("NSGA", pToken) == 0) {
		globalSetup->gaType = NSGA;
	} else {
		fclose(fInput);
		printf("Unknown parameter! It should be either SGA or NSGA\n");
		exit(1);
	}

	// decision variables

	// read the number of decision variables
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	// the number can't be less than or equal to zero
	if ((globalSetup->noOfDecisionVariables = atoi(pToken)) <= 0) {
		fclose(fInput);
		printf("Error! number of decision variables should be > 0\n");
		exit(1);
	}

	globalSetup->variableTypes = new VariableType[globalSetup->noOfDecisionVariables];
	globalSetup->variableRanges = new double*[(globalSetup->noOfDecisionVariables)];
	for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
		// read a line
		if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}
		// variable type
		if (strcmp("double", pToken) == 0) {
			globalSetup->variableTypes[ii] = Real;
		} else if (strcmp("int", pToken) == 0) {
			globalSetup->variableTypes[ii] = Integer;
		} else {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}

		// variable ranges
		// allocate memory
		globalSetup->variableRanges[ii] = new double[2];

		// lower bound
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}
		globalSetup->variableRanges[ii][0] = atof(pToken);

		// upper bound
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}
		globalSetup->variableRanges[ii][1] = atof(pToken);
		if (globalSetup->variableRanges[ii][1]
				<= globalSetup->variableRanges[ii][0]) {
			fclose(fInput);
			printf(
					"Error! lower bound, %f, must be lower than the upper bound, %f\n",
					globalSetup->variableRanges[ii][0],
					globalSetup->variableRanges[ii][1]);
			exit(1);
		}
	}

	// objectives

	// read the number of objectives
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	// the number can't be less than or equal to zero
	if ((globalSetup->noOfRawObjectives = atoi(pToken)) <= 0) {
		fclose(fInput);
		printf("Error! number of objectives should be > 0\n");
		exit(1);
	}
	globalSetup->finalNoOfObjectives = globalSetup->noOfRawObjectives;
	globalSetup->noOfLinearObjectiveCombinations = 0;

	globalSetup->typeOfOptimizations
			= new OptimType[globalSetup->finalNoOfObjectives];
	for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
		// read a line
		if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}
		// optimization type
		if (strcmp("Min", pToken) == 0) {
			globalSetup->typeOfOptimizations[ii] = Minimization;
		} else if (strcmp("Max", pToken) == 0) {
			globalSetup->typeOfOptimizations[ii] = Maximization;
		} else {
			fclose(fInput);
			printf("Error! optimization type can either be Min or Max\n");
			exit(1);
		}
	}

	// constrained variables

	// read the number of constrained variables
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	// the number can't be less than zero
	if ((globalSetup->noOfRawConstraints = atoi(pToken)) < 0) {
		fclose(fInput);
		printf("Error! number of constraints should be >= 0\n");
		exit(1);
	} else if (globalSetup->noOfRawConstraints == 0) {
		if ((strlen(pToken) != 1) || (*pToken != '0')) {
			fclose(fInput);
			printf("Error! number of constraints should be >= 0\n");
			exit(1);
		}
	}
	globalSetup->noOfLinearConstraintCombinations = 0;
	globalSetup->finalNoOfConstraints = globalSetup->noOfRawConstraints;

	// penalty weights
	globalSetup->penaltyWeights = new double[globalSetup->finalNoOfConstraints];
	for (ii = 0; ii < globalSetup->finalNoOfConstraints; ii++) {
		// read a line
		if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
			fclose(fInput);
			printf("Error in the input file, please refer to the documentation\n");
			exit(1);
		}
		globalSetup->penaltyWeights[ii] = atof(pToken);
	}

	// general parameters
	// population size
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		// Use a default value of 30*ell*log(ell);
		printf("Using default population-sizing thumbrule: n = 30*ell*log(ell)\n");

		// Take care of small problem sizes where log(ell) < 1
		if (globalSetup->noOfDecisionVariables > 2)
			globalSetup->populationSize = (int)(30
					*(globalSetup->noOfDecisionVariables)
					*log((double)(globalSetup->noOfDecisionVariables)));
		else
			globalSetup->populationSize = (int)(30
					*(globalSetup->noOfDecisionVariables));

		//Round it to next nearest tenth number
		if ((globalSetup->populationSize)%10)
			globalSetup->populationSize += (globalSetup->populationSize)%10;
		printf("The population size used is: %d\n", globalSetup->populationSize);
	}
	// the number can't be less than or equal to zero
	else if ((globalSetup->populationSize = atoi(pToken)) <= 0) {
		fclose(fInput);
		printf("The population size must be > 0\n");
		exit(1);
	} else if ((globalSetup->populationSize % 2) != 0) {
		// the number can't be an odd number
		fclose(fInput);
		printf("Error! population size must be an even number\n");
		exit(1);
	}
	// maximum generations
	// the number can't be less than or equal to zero
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		// Use a default value of 6*ell;
		printf("Using default convergence-time thumbrule: tc = 6*ell\n");
		globalSetup->maxGenerations = 6*(globalSetup->noOfDecisionVariables);
		if ((globalSetup->maxGenerations)%10)
			globalSetup->maxGenerations += 10-(globalSetup->maxGenerations)%10;
		printf("The maximum number of generations set is: %d\n",
				globalSetup->maxGenerations);
	} else if ((globalSetup->maxGenerations = atoi(pToken)) <= 0) {
		fclose(fInput);
		printf("Error! maximum number of generations must be > 0\n");
		exit(1);
	}
	//cout << "test: " << globalSetup->maxGenerations << endl;
	// replace proportion
	// the number should be in (0.0, 1.0]
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		// Setting a default value of 0.9
		printf("Using default replacement proportion of 0.9\n");
		globalSetup->replaceProportion = 0.9;
	} else if (((globalSetup->replaceProportion = atof(pToken)) <= 0.0)
			|| (globalSetup->replaceProportion > 1.0)) {
		fclose(fInput);
		printf("Error! proportion of parent population that should be replaced must be > 0 and <= 1\n");
		exit(1);
	}

	// niching (multimodal handling)
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		// Using NoNiching by default
		printf("No niching method is used by default\n");
		globalSetup->nichingType = NoNiching;
	}
	// niching type
	else if (strcmp("NoNiching", pToken) == 0) {
		globalSetup->nichingType = NoNiching;
		globalSetup->nichingParameters = NULL;

	} else if (strcmp("Sharing", pToken) == 0) {
		globalSetup->nichingType = Sharing;
	} else if (strcmp("RTS", pToken) == 0) {
		globalSetup->nichingType = RTS;
	} else if (strcmp("DeterministicCrowding", pToken) == 0) {
		globalSetup->nichingType = DeterministicCrowding;
	} else {
		fclose(fInput);
		printf("Error! valid niching types are: NoNiching, Sharing, RTS, and DeterministicCrowding\n");
		exit(1);
	}
	// check niching type
	if ((globalSetup->gaType == NSGA)
			&& (globalSetup->nichingType != NoNiching)) {
		fclose(fInput);
		printf("Error! valid choice for niching types with NSGA is: NoNiching\n");
		exit(1);
	}
	// read niching parameters
	switch (globalSetup->nichingType) {
	case NoNiching:
	case DeterministicCrowding:
		// no extra parameters
		break;
	case Sharing: {
		globalSetup->nichingParameters = new double[2];
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default sharing radius of 4.24\n");
			((double *)(globalSetup->nichingParameters))[0] = 4.24;
		} else
			((double *)globalSetup->nichingParameters)[0]
					= atof(pToken);
		if (((double*)globalSetup->nichingParameters)[0] <= 0.0) {
			fclose(fInput);
			printf("Error! niching radius must be > 0\n");
			exit(1);
		}
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default scaling of 1 for fitness sharing\n");
			((double*)(globalSetup->nichingParameters))[1] = 1.0;
		} else
			((double*)(globalSetup->nichingParameters))[1] = atof(pToken);
	}
		break;
	case RTS: {
		globalSetup->nichingParameters = new int[1];

		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Setting the window size for RTS using default rule: w = min(ell, n/20)\n");
			if (globalSetup->noOfDecisionVariables
					< (globalSetup->populationSize/20))
				((int*)(globalSetup->nichingParameters))[0]
						= globalSetup->noOfDecisionVariables;
			else
				((int*)(globalSetup->nichingParameters))[0]
						= (globalSetup->populationSize)/20;

			//Adjust for small problem sizes: w = n/20;
			if (globalSetup->noOfDecisionVariables < 20)
				((int*)(globalSetup->nichingParameters))[0]
						= (globalSetup->populationSize)/20;

			//Check if the window size is greater than the population size
			if (((int*)globalSetup->nichingParameters)[0]
					> globalSetup->populationSize)
				((int*)(globalSetup->nichingParameters))[0]
						= globalSetup->populationSize;

			printf("The window size used for RTR is: %d\n", ((int*)globalSetup->nichingParameters)[0]);
		} else
			((int*)(globalSetup->nichingParameters))[0] = atoi(pToken);

		// the window size should be in (0, populationSize]
		if ((((int*)globalSetup->nichingParameters)[0] <= 0) || (((int*)globalSetup->nichingParameters)[0]
				> globalSetup->populationSize)) {
			fclose(fInput);
			printf("Error! window size for RTR should be > 0 and less than the population size\n");
			exit(1);
		}
	}
		break;
	default: {
		fclose(fInput);
		printf("Error! valid choices for niching type are: NoNiching, Sharing, RTS, and DeterministicCrowding\n");
		exit(1);
	}
		break;
	}

	// selection
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using tournament selection w/o replacement as a default selection method\n");
		globalSetup->selectionType = TournamentWOR;
	}
	// selection type
	else if (strcmp("TournamentWOR", pToken) == 0) {
		globalSetup->selectionType = TournamentWOR;
	} else if (strcmp("SUS", pToken) == 0) {
		globalSetup->selectionType = SUS;
	} else if (strcmp("Truncation", pToken) == 0) {
		globalSetup->selectionType = Truncation;
	} else if (strcmp("RouletteWheel", pToken) == 0) {
		globalSetup->selectionType = RouletteWheel;
	} else if (strcmp("TournamentWR", pToken) == 0) {
		globalSetup->selectionType = TournamentWR;
	} else {
		fclose(fInput);
		printf("Error! valid selection methods are: RouletteWheel, SUS, TournamentWOR, TournamentWR, and Truncation\n");
		exit(1);
	}

	// check selection type
	if ((globalSetup->gaType == NSGA) && ((globalSetup->selectionType == SUS)
			|| (globalSetup->selectionType == RouletteWheel))) {
		fclose(fInput);
		printf("Error! with NSGA, valid selection methods are: TournamentWOR, TournamentWR, and Truncation\n");
		exit(1);
	}
	// read selection parameters
	switch (globalSetup->selectionType) {
	case TournamentWOR:
	case Truncation:
	case TournamentWR: {
		globalSetup->selectionParameters = new int[1];

		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default tournament size of 2\n");
			((int*)globalSetup->selectionParameters)[0] = 2;
		} else
			((int*)globalSetup->selectionParameters)[0]
					= atoi(pToken);
	}
		break;
	case SUS:
	case RouletteWheel:
		// no extra parameters
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid selection parameter\n");
		exit(1);
	}
		break;
	}

	// Crossover
	// crossover probability
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using a default crossover probability of 0.9\n");
		globalSetup->xOverProbability = 0.9;
	}
	// the number should be [0.0, 1.0]
	else if (((globalSetup->xOverProbability = atof(pToken)) < 0.0)
			|| (globalSetup->xOverProbability > 1.0)) {
		fclose(fInput);
		printf("Error! crossover probability must be >= 0.0 and <= 1.0\n");
		exit(1);
	}
	// crossover type
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using SBX as the default crossover method. Note that this might be an inappropriate choice if your variables are binary\n");
		globalSetup->xOverType = SBX;
	} else if (strcmp("OnePoint", pToken) == 0) {
		globalSetup->xOverType = OnePoint;
	} else if (strcmp("TwoPoint", pToken) == 0) {
		globalSetup->xOverType = TwoPoint;
	} else if (strcmp("Uniform", pToken) == 0) {
		globalSetup->xOverType = Uniform;
	} else if (strcmp("SBX", pToken) == 0) {
		globalSetup->xOverType = SBX;
	} else {
		fclose(fInput);
		printf("Error! valid crossover types are: OnePoint, TwoPoint, Uniform, and SBX\n");
		exit(1);
	}
	// read crossover parameters
	switch (globalSetup->xOverType) {
	case OnePoint:
	case TwoPoint:
		//no extra parameters
		break;
	case Uniform: {
		globalSetup->xOverParameters = new double[1];
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default genewise swap probability of 0.5\n");
			((double*)globalSetup->xOverParameters)[0] = 0.5;
		} else
			((double*)globalSetup->xOverParameters)[0]
					= atof(pToken);
		if ((((double*)globalSetup->xOverParameters)[0] <= 0.0)||(((double*)globalSetup->xOverParameters)[0] >= 1.0)) {
			fclose(fInput);
			printf("Genewise probability must be > 0.0 and < 1.0\n");
			exit(1);
		}
	}
	case SBX: {
		globalSetup->xOverParameters = new double[2];

		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default genewise swap probability of 0.5\n");
			((double*)globalSetup->xOverParameters)[0] = 0.5;
		} else
			((double*)globalSetup->xOverParameters)[0]
					= atof(pToken);
		if ((((double*)globalSetup->xOverParameters)[0] <= 0.0)||(((double*)globalSetup->xOverParameters)[0] >= 1.0)) {
			fclose(fInput);
			printf("Error! genewise probability must be > 0.0 and < 1.0\n");
			exit(1);
		}
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default polynomial order for SBX: 10\n");
			((double*)globalSetup->xOverParameters)[1] = 10;
		} else
			((double*)globalSetup->xOverParameters)[1]
					= atof(pToken);
		if (((double*)globalSetup->xOverParameters)[1] < 0.0) {
			fclose(fInput);
			printf("Error! genewise probability must be >= 0.0\n");
			exit(1);
		}
	}
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid crossover parameter\n");
		exit(1);
	}
		break;
	}

	// mutation
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using a default mutation probability of 0.1\n");
		globalSetup->mutationProbability = 0.1;
	}
	// mutation probability
	// the number should be [0.0, 1.0]
	else if (((globalSetup->mutationProbability = atof(pToken)) < 0.0)
			|| (globalSetup->mutationProbability > 1.0)) {
		fclose(fInput);
		printf("Error! mutation probability must be >= 0.0 and <= 1.0.\n");
		exit(1);
	}
	// mutation type
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using Polynomial as the default mutation method. Note that this might be an inappropriate choice if your variables are binary\n");
		globalSetup->mutationType = Polynomial;
	} else if (strcmp("Selective", pToken) == 0) {
		globalSetup->mutationType = Selective;
	} else if (strcmp("Genewise", pToken) == 0) {
		globalSetup->mutationType = Genewise;
		//printf("got %s\n",pToken);
	} else if (strcmp("Polynomial", pToken) == 0) {
		globalSetup->mutationType = Polynomial;
	} else {
		fclose(fInput);
		printf("Error! valid mutation types are: Selective, Genewise, and Polynomial\n");
		exit(1);
	}
	// read mutation parameters
	switch (globalSetup->mutationType) {
	case Selective:
		// no extra parameters
		break;
	case Polynomial: {
		globalSetup->mutationParameters = new int[1];
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using a default value for the polynomial probability: 20\n");
			((int*)globalSetup->mutationParameters)[0] = 20;
		} else
			((int*)globalSetup->mutationParameters)[0]
					= atoi(pToken);
		if (((int*)globalSetup->mutationParameters)[0] < 0) {
			fclose(fInput);
			printf("Error! polynomial order for polynomial mutation must be > 0\n");
			exit(1);
		}
	}
		break;
	case Genewise: {
		globalSetup->mutationParameters
				= new double[globalSetup->noOfDecisionVariables];

		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
				printf("Using default std. deviation of 10 percent of the variable range, not %s, for #%i\n", pToken, ii);
				((double*)globalSetup->mutationParameters)[ii]
						= 0.1*(globalSetup->variableRanges[ii][1]
								- globalSetup->variableRanges[ii][0]);
			} else
				((double*)globalSetup->mutationParameters)[ii]
						= atof(pToken);
			if (((double*)globalSetup->mutationParameters)[ii] <= 0.0) {
				fclose(fInput);
				printf(
						"Error! standard deviation for gene %d for genewise mutation must be > 0\n",
						ii);
				exit(1);
			}
		}
	}
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid mutation parameter\n");
		exit(1);
	}
		break;
	}

	// scaling method
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Scaling is not used by default\n");
		globalSetup->scalingMethod = NoScaling;
	}
	// scaling method
	else if (strcmp("NoScaling", pToken) == 0) {
		globalSetup->scalingMethod = NoScaling;
	} else if (strcmp("Ranking", pToken) == 0) {
		globalSetup->scalingMethod = Ranking;
	} else if (strcmp("SigmaScaling", pToken) == 0) {
		globalSetup->scalingMethod = SigmaScaling;
	} else {
		fclose(fInput);
		printf("Error! valid scaling methods are: NoScaling, Ranking, and SigmaScaling, not %s\n",pToken);
		exit(1);
	}
	// read scaling parameters
	switch (globalSetup->scalingMethod) {
	case NoScaling:
	case Ranking:
		// no extra parameters
		break;
	case SigmaScaling: {
		globalSetup->scalingParameters = new double[1];

		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			((double*)globalSetup->scalingParameters)[0] = 1.0;
		} else
			((double*)globalSetup->scalingParameters)[0]
					= atof(pToken);
		if (((double*)globalSetup->scalingParameters)[0] <= 0.0) {
			fclose(fInput);
			printf("Error! scaling parameter for SigmaScaling must be > 0\n");
			exit(1);
		}
	}
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid scaling parameter\n");
		exit(1);
	}
		break;
	}

	// constraint method
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		if (globalSetup->finalNoOfConstraints == 0) {
			printf("Using no constraint handling method by default\n");
			globalSetup->constraintMethod = NoConstraints;
		} else {
			printf("Using tournament selection as the default constraint handling method\n");
			globalSetup->constraintMethod = Tournament;
		}
	}
	// constraint method
	else if (strcmp("NoConstraints", pToken) == 0) {
		globalSetup->constraintMethod = NoConstraints;
	} else if (strcmp("Penalty", pToken) == 0) {
		globalSetup->constraintMethod = Penalty;
	} else if (strcmp("Tournament", pToken) == 0) {
		globalSetup->constraintMethod = Tournament;
	} else {
		fclose(fInput);
		printf("Error! valid constraint handling methods are: NoConstraint, Penalty, and Tournament\n");
		exit(1);
	}
	// check constraint method
	if ((globalSetup->gaType == NSGA) && (globalSetup->constraintMethod
			== Penalty)) {
		fclose(fInput);
		printf("Error! penalty based constraint handling method cannot be used with NSGA\n");
		exit(1);
	}
	if ((globalSetup->finalNoOfConstraints == 0)
			&& (globalSetup->constraintMethod != NoConstraints)) {
		fclose(fInput);
		printf("Error! valid constraint-handling method when there are no constraints is NoConstraints\n");
		exit(1);
	}
	// read penalty function
	switch (globalSetup->constraintMethod) {
	case NoConstraints:
	case Tournament:
		// no extra parameters
		break;
	case Penalty: {
		// penalty function
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using linear penalty be default\n");
			globalSetup->penaltyFunction = Linear;
		} else if (strcmp("Linear", pToken) == 0) {
			globalSetup->penaltyFunction = Linear;
		} else if (strcmp("Quadratic", pToken) == 0) {
			globalSetup->penaltyFunction = Quadratic;
		} else {
			fclose(fInput);
			printf("Error! valid penalty function methods are: Linear and Quadratic\n");
			exit(1);
		}
	}
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid constraint-handling method parameter\n");
		exit(1);
	}
		break;
	}

	// local search method
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using no local search by default\n");
		globalSetup->localSearchMethod = NoLocalSearch;
	}
	// local search method
	else if (strcmp("NoLocalSearch", pToken) == 0) {
		globalSetup->localSearchMethod = NoLocalSearch;
	} else if (strcmp("SimplexSearch", pToken) == 0) {
		globalSetup->localSearchMethod = SimplexSearch;
	} else {
		fclose(fInput);
		printf("Error! valid local search methods are: NoLocalSearch and SimplexSearch\n");
		exit(1);
	}
	// check local search method
	if ((globalSetup->localSearchMethod != NoLocalSearch)
			&& (globalSetup->gaType == NSGA)) {
		fclose(fInput);
		printf("Error! cannot use local search with NSGA\n");
		exit(1);
	}
	switch (globalSetup->localSearchMethod) {
	case NoLocalSearch:
		// no extra parameters
		break;
	case SimplexSearch:
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default tolerance of 0.001\n");
			globalSetup->maxLocalTolerance = 1.0E-3;
		} else
			globalSetup->maxLocalTolerance = atof(pToken);
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Setting maximum local search evaluations to 20\n");
			globalSetup->maxLocalEvaluations = 20;
		} else
			globalSetup->maxLocalEvaluations = atoi(pToken);
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default local penalty parameter of 1.0\n");
			globalSetup->initialLocalPenaltyParameter = 1.0;
		} else
			globalSetup->initialLocalPenaltyParameter = atof(pToken);
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default local update parameter of 2.0\n");
			globalSetup->localUpdateParameter = 2.0;
		} else
			globalSetup->localUpdateParameter = atof(pToken);
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default Lamarckian probability of 0.15\n");
			globalSetup->lamarckianProbability = 0.15;
		} else
			globalSetup->lamarckianProbability = atof(pToken);
		if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
			printf("Using default overall local search probability of 0.5\n");
			globalSetup->localSearchProbability = 0.5;
		} else
			globalSetup->localSearchProbability = atof(pToken);
		break;
	default: {
		fclose(fInput);
		printf("Error! invalid local search parameters\n");
		exit(1);
	}
		break;
	}

	// stopping criteria

	// noOfStoppingCriterias
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using no extra stopping criterias by default.\n");
		globalSetup->noOfStoppingCriterias = 0;
	}
	// the number can't be less than zero
	else if ((globalSetup->noOfStoppingCriterias = atoi(pToken)) < 0) {
		fclose(fInput);
		printf("Error! number of stopping criterias must be > 0\n");
		exit(1);
	} else if (globalSetup->noOfStoppingCriterias == 0) {
		if ((strlen(pToken) != 1) || (*pToken != '0')) {
			fclose(fInput);
			printf("Error! number of stopping criterias must be a number!\n");
			exit(1);
		}
	}

	// allocate memory
	if (globalSetup->noOfStoppingCriterias > 0) {
		// number of generation window
		if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
			printf("The default window size used in stopping critiria is: 5.\n");
			globalSetup->genNumWindow = 5;
		} else
			globalSetup->genNumWindow = atoi(pToken);

		// the number should be greater than zero  
		if (globalSetup->genNumWindow <= 0) {
			fclose(fInput);
			printf("Error! window size used in stopping critieria must be > 0\n");
			exit(1);
		}
		globalSetup->otherStoppingCriteria
				= new StoppingCriterias[globalSetup->noOfStoppingCriterias];
		globalSetup->stoppingParameter
				= new double[globalSetup->noOfStoppingCriterias];

		for (ii = 0; ii < globalSetup->noOfStoppingCriterias; ii++) {
			// read a line
			if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
				fclose(fInput);
				printf("Error in the input file, Please see the documentation.\n");
				exit(1);
			}
			// stopping criterion type
			if (strcmp("NoOfEvaluations", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = NoOfEvaluations;
			} else if (strcmp("FitnessVariance", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = FitnessVariance;
			} else if (strcmp("AverageFitness", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = AverageFitness;
			} else if (strcmp("AverageObjective", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = AverageObjective;
			} else if (strcmp("ChangeInBestFitness", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInBestFitness;
			} else if (strcmp("ChangeInAvgFitness", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInAvgFitness;
			} else if (strcmp("ChangeInFitnessVar", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInFitnessVar;
			} else if (strcmp("ChangeInBestObjective", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInBestObjective;
			} else if (strcmp("ChangeInAvgObjective", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInAvgObjective;
			} else if (strcmp("NoOfFronts", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = NoOfFronts;
			} else if (strcmp("NoOfGuysInFirstFront", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = NoOfGuysInFirstFront;
			} else if (strcmp("ChangeInNoOfFronts", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = ChangeInNoOfFronts;
			} else if (strcmp("BestFitness", pToken) == 0) {
				globalSetup->otherStoppingCriteria[ii] = BestFitness;
			} else {
				fclose(fInput);
				printf("Error! invalid stopping criteria parameter\n");
				exit(1);
			}

			// check criteria
			switch (globalSetup->otherStoppingCriteria[ii]) {
			case NoOfEvaluations:
			case FitnessVariance:
				// 1 - nothing to check
				break;

			case AverageFitness:
			case AverageObjective:
			case ChangeInBestFitness:
			case ChangeInAvgFitness:
			case ChangeInFitnessVar:
			case ChangeInBestObjective:
			case ChangeInAvgObjective:
			case BestFitness:
				// 2 - SGA w/o MM only
				if ((globalSetup->gaType != SGA) || (globalSetup->nichingType
						!= NoNiching)) {
					fclose(fInput);
					printf("Cannot use the following stopping critier with NSGA or when Niching is used:\n");
					printf("\t AverageFitness\n\t AverageObjective\n");
					printf("\t ChangeInBestFitness\n\t ChangeInAvgFitness\n");
					printf("\t ChangeInFitnessVar\n\t ChangeInBestObjective\n");
					printf("\t ChangeInAvgObjective\n\t BestFitness\n");
					printf("Error! invalid choice of stopping criteria");
					exit(1);
				}
				break;
			case NoOfFronts:
			case NoOfGuysInFirstFront:
			case ChangeInNoOfFronts:
				// 3 - NSGA only
				if (globalSetup->gaType != NSGA) {
					fclose(fInput);
					printf("Error! following stopping criteria can be used only with NSGA:\n");
					printf("\t NoOfFronts\n\t NoOfGuysInFirstFront\n\t ChangeInNoOfFronts\n");
					exit(1);
				}
				break;

			default: {
				fclose(fInput);
				printf("Error! invalid stopping criteria\n");
				exit(1);
			}
				break;
			}
			// parameter
			if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
				fclose(fInput);
				printf("Error! invalid stopping criteria parameter\n");
				exit(1);
			}
			((double*)globalSetup->stoppingParameter)[ii]
					= atof(pToken);
		}
	}
	// Initial population generation

	// Load population
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Using random initialization by default.\n");
		globalSetup->loadPopulation = false;
		globalSetup->evaluateAgain = true;
	}
	// Check if the load population is 0 or 1
	else {//
		// load the population from a file
		if (atoi(pToken) == 1) {
			globalSetup->loadPopulation = true;

			globalSetup->populationFileName = new char[256];
			// Read the name of the file to load the population from
			if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
				fclose(fInput);
				printf("Error! invalid file name to load the population from, see documentation for more details.\n");
				exit(1);
			}
			strcpy(globalSetup->populationFileName, pToken);

			// Evaluate the population or not
			if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
				printf("Evaluating the initial population by default.\n");
				globalSetup->evaluateAgain = true;
			} else if (atoi(pToken) == 0) {
				globalSetup->evaluateAgain = false;
			} else
				globalSetup->evaluateAgain = true;
		}
		// Use random initialization
		else {
			globalSetup->loadPopulation = false;
			globalSetup->evaluateAgain = true;
		}
	}

	// Save population
	if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
		fclose(fInput);
		printf("Error in the input file, please refer to the documentation\n");
		exit(1);
	}
	if (strcmp("default", pToken) == 0) {
		printf("Saving the evaluated individuals to a file by default.\n");
		globalSetup->savePopulation = true;
		globalSetup->saveEvalSolutions = new char[256];
		printf("Enter the filename you want to save the population to.\n");
		fflush(stdout);
		scanf("%s", globalSetup->saveEvalSolutions);
	}
	// Check if the save population is 0 or 1
	else {
		fflush(stdout);
		// Save the population to a file
		if (atoi(pToken) == 1) {
			globalSetup->savePopulation = true;
			// Filename to save the evaluated solutions to.
			globalSetup->saveEvalSolutions = new char[256];
			if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
				fclose(fInput);
				printf("Error in the input file, please refer to the documentation\n");
				exit(1);
			}
			strcpy(globalSetup->saveEvalSolutions, pToken);
		}
	}
	fclose(fInput);
	fflush(stdout);

	if (globalSetup->savePopulation) {
		fOutput = fopen(globalSetup->saveEvalSolutions, "w");
		fprintf(fOutput, "%% ");
		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			fprintf(fOutput, "var #%d\t", ii);
		}
		for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
			fprintf(fOutput, "obj #%d\t", ii);
		}
		if (globalSetup->finalNoOfConstraints > 0) {
			for (ii = 0; ii < globalSetup->finalNoOfConstraints; ii++) {
				fprintf(fOutput, "const #%d\t", ii);
			}
			fprintf(fOutput, "penalty");
		}
		fprintf(fOutput, "\n");
		fclose(fOutput);
	}
	strcpy(Clone_Dir, Unique_Dir);
	char *output_cat1 = strcat(Clone_Dir, OutAll_Indivs);
	fOutput = fopen(output_cat1, "w"); // create this file newly, in case it already exists cuz later calls are APPENDS
	if (globalSetup->gaType == SGA) {
		printf("Instantiating genetic algorithm, single objective...\n");
	}
	printf("Evaluating initial population, but not saving this data\n");
	int *gen_current_pointer=&GA_gen_current;
	GA ga;
	globalHaveDonePreEval=1; // only after those initial evaluations, start saving data
	printf("\n");
	if (globalSetup->gaType == SGA) {
		while (ga.generate(globalDataParams, numCells*9/*6*/, gen_current_pointer))
			;
	} else {
		while (ga.nsgaGenerate(globalDataParams, numCells*9/*6*/, globalDataObjective0, globalDataObjective1, gen_current_pointer))
			;
	}
	//StoreBestIndividuals(); // commented out 3/17 to use a smaller array, and save these Best not at end but per generation
	DisplayFinalOutput(); 
	
	delete [](globalSetup->penaltyWeights);
	delete [](globalSetup->variableTypes);
	for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
		delete [](globalSetup->variableRanges[ii]);
	}
	delete [](globalSetup->variableRanges);
	delete [](globalSetup->typeOfOptimizations);
	if ((globalSetup->selectionType == TournamentWR)
			||(globalSetup->selectionType == TournamentWOR)
			||(globalSetup->selectionType == Truncation)) {
		delete [](int*)(globalSetup->selectionParameters);
	}
	if ((globalSetup->xOverType==Uniform)||(globalSetup->xOverType == SBX))
		delete [](double*)(globalSetup->xOverParameters);
	if (globalSetup->mutationType == Polynomial)
		delete [](int*)(globalSetup->mutationParameters);
	else if (globalSetup->mutationType == Genewise)
		delete [](double*)(globalSetup->mutationParameters);
	if (globalSetup->nichingType == RTS)
		delete [](int*)(globalSetup->nichingParameters);
	if (globalSetup->nichingType == Sharing)
		delete [](double*)(globalSetup->nichingParameters);
	if (globalSetup->scalingMethod == SigmaScaling)
		delete [](double*)(globalSetup->scalingParameters);
	if (globalSetup->noOfStoppingCriterias > 0) {
		delete [] globalSetup->otherStoppingCriteria;
		delete [](double*)(globalSetup->stoppingParameter);
	}
	if (globalSetup->loadPopulation)
		delete [](globalSetup->populationFileName);
	if (globalSetup->savePopulation)
		delete [](globalSetup->saveEvalSolutions);

	return 1;
}


//////////////////////////////////

void AlterBoxCoordinates(double my_x[], int my_h)
{
	// i'm gonna take a chance and do this the easy elegant way
	double highest_previous=0;
	
	// 1st loop: find height of last set
	for (int xxj=0; xxj < numCells*9; xxj++)
	{
		// looking for largest value in a z-value slot
		if ( ((xxj+1)%3==0) and (my_x[xxj] > highest_previous) ) {	
			highest_previous = my_x[xxj];	
		}
	}

	// 2nd loop: replace heights and half-heights
	for (int xxj=0; xxj < numCells*9; xxj++)
	{
		if ((xxj+1)%3==0) {
			if (my_x[xxj] == highest_previous) {my_x[xxj] = my_h;}		// a z value at top height	
		}
	}
}

void AlterFunnelCoordinates(double my_x[], int my_h)
{
	// i'm gonna take a chance and do this the easy elegant way
	double highest_previous=0;
	
	// 1st loop: find height of last set
	for (int xxj=0; xxj < numCells*9; xxj++)
	{
		// looking for largest value in a z-value slot
		if ( ((xxj+1)%3==0) and (my_x[xxj] > highest_previous) ) {	
			highest_previous = my_x[xxj];	
		}
	}

	// 2nd loop: replace heights and half-heights
	for (int xxj=0; xxj < numCells*9; xxj++)
	{
		// first check if it is a z-value! otherwise it may alter x,y dimensions, undesirable!!!
		if ((xxj+1)%3==0) {
			if (my_x[xxj] == highest_previous) {my_x[xxj] = my_h;}		// a z value at top height
			if (my_x[xxj] == highest_previous/2) {my_x[xxj] = my_h/2;}	// a z value in that center point	
		}
	}
}

void AlterFunnelWaist(double my_x[], int my_y, int my_z)
{
	// i'm gonna take a chance and do this the easy elegant way
	//double x_previous=0;
	//double y_previous=0;
	
	// 1st loop: find 
	for (int xxj=1; xxj < numCells*9; xxj++)
	{
		// looking 
		if ((xxj+2)%3==0) { // first, is this a y-value?
			if ((my_x[xxj-1]) == xBound/2) { // next, is the vertex's x value in the center?
				my_x[xxj] = (double)my_y;
			}	
		}

		if ((xxj+1)%3==0) { // first, is this a z-value?
			if ((my_x[xxj-2]) == xBound/2) { // next, is the vertex's x value in the center of structure?
				my_x[xxj] = (double)my_z;
			}	
		}
	}
}

void AlterRoofStructure(double output_roof[], double unit_coords[], double roof_pos[], int num_insts, int copy_numcells)
{
	double tempadd =0;
	
	for (int ars=0; ars<num_insts; ars++) {
		for (int cel=0; cel<copy_numcells; cel++) {
			for(int var=0; var<9; var++) {
				
				if ( (var==0) or (var==3) || (var==6) ) {tempadd = roof_pos[ars*2+0];} // x
				if ( (var==1) or (var==4) || (var==7) ) {tempadd = roof_pos[ars*2+1];} // y
				if ( (var==2) or (var==5) || (var==8) ) {tempadd = 0;} //z

				// now set coordinates of current instance to the regular UNIT coords + translation (tempadd)
				output_roof[ars*copy_numcells*9 + cel*9 + var] = unit_coords[cel*9 +var] + tempadd;
			}
		}
	}
}


//int main() {
//	input_i = MIT;
//
//	globalCellsOneSided = false;
//	cellEfficiency = 0.06;
//	numCells = 12;
//
//	//dimension of the solar panel
//	double l = 0.624;
//	double w = 0.468;
//
//	raysPerCell = 100;
//	Get2DTriangleGrid(); // precalculate grid points
//
//	longitude = INPUT_examples[input_i][0];   // longitude of location
//	latitude    = INPUT_examples[input_i][1];   // latitude of location
//	avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//	avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//	elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//	timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//
//	xBound = 10;
//	yBound = 10;
//	zBound = 10;
//
//	startMonth = 8;
//	startDay = 30;
//	startYear = 2010;
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	double modelCoordinates[9*numCells];
//	modelPointer = modelCoordinates;
//	Cell_3 cellArray[numCells];
//
//	{
//		{
//			{
//				double m[] = {
//						0.000001,-0.000001,10.000001,10.000000,10.000000,10.000001,-0.000001,10.000000,10.000001,
//						0.000001,-0.000001,10.000001,10.000000,0.000000,10.000002,10.000000,10.000000,10.000001,
//						10.000000,0.000000,0.000000,-0.000001,10.000000,0.000000,10.000000,10.000000,-0.000000,
//						10.000000,0.000000,0.000000,0.000001,-0.000000,0.000001,-0.000001,10.000000,0.000000,
//						-0.000001,10.000000,10.000001,-0.000001,10.000000,0.000000,0.000001,-0.000001,10.000001,
//						-0.000001,10.000000,0.000000,0.000001,-0.000000,0.000001,0.000001,-0.000001,10.000001,
//						10.000000,10.000000,-0.000000,10.000000,10.000000,10.000001,10.000000,0.000000,0.000000,
//						10.000000,10.000000,10.000001,10.000000,0.000000,10.000002,10.000000,0.000000,0.000000,
//						0.000001,-0.000000,0.000001,10.000000,0.000000,0.000000,0.000001,-0.000001,10.000001,
//						10.000000,0.000000,0.000000,10.000000,0.000000,10.000002,0.000001,-0.000001,10.000001,
//						10.000000,10.000000,-0.000000,-0.000001,10.000000,0.000000,10.000000,10.000000,10.000001,
//						-0.000001,10.000000,0.000000,-0.000001,10.000000,10.000001,10.000000,10.000000,10.000001,
//				};
//
//				modelPointer = m;
//
//				int v_index_triflat[2][3] = {{0,1,2},{0,2,1}};  // stores integer indices for 3 verts of each triangle
//
//				Triangle *temptri;
//				double pn[3];
//				double newvertex1[3];
//				double newvertex2[3];
//				double newvertex3[3];
//
//				for(int i=0; i<numCells; i++) {
//					Cell_3 *cel = &cellArray[i];
//
//					cel->numVertices = 3;
//					cel->numTriangles = 2; // like TRIFLAT, it has a top and bottom triangle
//
//
//					for (int nt=0; nt < cel->numTriangles; nt++)
//					{
//						cel->triangleList[nt] = new Triangle;
//					}
//
//					cel->vertices[0] = new Position;
//					cel->vertices[0]->x = modelPointer[i*9+0]; // get coordinate from global pointer to array modelCoordinates
//					cel->vertices[0]->y = modelPointer[i*9+1];
//					cel->vertices[0]->z = modelPointer[i*9+2];
//
//					cel->vertices[1] = new Position;
//					cel->vertices[1]->x = modelPointer[i*9+3];
//					cel->vertices[1]->y = modelPointer[i*9+4];
//					cel->vertices[1]->z = modelPointer[i*9+5];
//
//					cel->vertices[2] = new Position;
//					cel->vertices[2]->x = modelPointer[i*9+6];
//					cel->vertices[2]->y = modelPointer[i*9+7];
//					cel->vertices[2]->z = modelPointer[i*9+8];
//
//
//					// for former complicated shapes there could be many more triangles. but a flat cell has two sides-->numTriangles=2
//					for (int m=0; m < (cel->numTriangles); ++m)
//					{
//						temptri = cel->triangleList[m];
//						temptri->index[0] = v_index_triflat[m][0];
//						temptri->index[1] = v_index_triflat[m][1];
//						temptri->index[2] = v_index_triflat[m][2];
//
//						//cal the normal to a particular face
//
//						pn[0]=0; pn[1]=0; pn[2]=0;
//						newvertex1[0] = (cel->vertices[temptri->index[0]])->x;
//						newvertex1[1] = (cel->vertices[temptri->index[0]])->y;
//						newvertex1[2] = (cel->vertices[temptri->index[0]])->z;
//						newvertex2[0] = (cel->vertices[temptri->index[1]])->x;
//						newvertex2[1] = (cel->vertices[temptri->index[1]])->y;
//						newvertex2[2] = (cel->vertices[temptri->index[1]])->z;
//						newvertex3[0] = (cel->vertices[temptri->index[2]])->x;
//						newvertex3[1] = (cel->vertices[temptri->index[2]])->y;
//						newvertex3[2] = (cel->vertices[temptri->index[2]])->z;
//						plane_exp_normal_3d(newvertex1, newvertex2, newvertex3, pn);
//						temptri->Normal[0]=pn[0]; temptri->Normal[1]=pn[1]; temptri->Normal[2]=pn[2];
//
//
//						temptri->vertList[0] = newvertex1[0]; temptri->vertList[1] = newvertex1[1]; temptri->vertList[2]=newvertex1[2];
//						temptri->vertList[3] = newvertex2[0]; temptri->vertList[4] = newvertex2[1]; temptri->vertList[5]=newvertex2[2];
//						temptri->vertList[6] = newvertex3[0]; temptri->vertList[7] = newvertex3[1]; temptri->vertList[8]=newvertex3[2];
//
//						temptri->Area = triangle_area_3d(temptri->vertList); // call to geometry.cpp function to get triangles area
//
//
//						if (m==0)
//						{
//							totalSurfaceArea+=temptri->Area;  // get the total area
//							temptri->is_positive=true;
//						} // finally, add area to the total structure area (m==0 => postive norm side only)
//						if (m==1) {
//							temptri->is_positive=false;
//						}
//					}
//
//				}
//
//				double day_energee=0;
//				day_energee = newCalculatePowerDayCycle(cellArray, modelPointer);   // calculate power at several instants in the day
//			}
//		}
//	}
//
//	return 0;
//}

//int main() {
//	input_i = MIT;
//
//	globalCellsOneSided = true;
//	cellEfficiency = 0.116;
//	numCells = 4;
//
//	//dimension of the solar panel
//	double l = 0.624;
//	double w = 0.468;
//
//	raysPerCell = 100;
//	Get2DTriangleGrid(); // precalculate grid points
//
//	longitude = INPUT_examples[input_i][0];   // longitude of location
//	latitude    = INPUT_examples[input_i][1];   // latitude of location
//	avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//	avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//	elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//	timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//
//	xBound = 1;
//	yBound = 1;
//	zBound = 1;
//
//	startMonth = 8;
//	startDay = 30;
//	startYear = 2010;
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	double modelCoordinates[9*numCells];
//	modelPointer = modelCoordinates;
//	Cell_3 cellArray[numCells];
//
////	std::ofstream myfile;
////	myfile.open("v_shape.txt");
//
//	double increment = PI/36.0;
//
////	for (double rho = 0; rho<2*PI; rho+=increment) //angle of rotation
//	double rho = 0*(PI/180);
//	{
////		for (double theta = 0; theta<=PI/2.0; theta+=increment) //angle of inclination
//		double theta = 0*(PI/180);
//		{
////			for (double phi = 0; phi<=PI; phi+=increment) //angle between two panels
//			double phi = 50*(PI/180);
//			{
//				double lt = w*sin(theta)*cos(phi/2.0); //find the offsets caused by the tilt
//				double wt = w*sin(phi/2.0);
//				double ht = w*cos(theta)*cos(phi/2.0);
//				//Create Tent
//				//use offsets and side lengths to calculate coordinates
//				double m[] = {0,lt,0, wt,0,ht, 0,lt+l*cos(theta),l*sin(theta),
//							  wt,0,ht, wt,l*cos(theta),ht+l*sin(theta), 0,lt+l*cos(theta),l*sin(theta),
//							  wt,0,ht, 2*wt,lt,0, 2*wt,lt+l*cos(theta),l*sin(theta),
//							  wt,0,ht, 2*wt,lt+l*cos(theta),l*sin(theta), wt,l*cos(theta),ht+l*sin(theta)};
//
//				//Create V
////						    double m[] = {0,0,ht, wt,lt,0, 0,l*cos(theta),ht+l*sin(theta),
////						    		      wt,lt,0, wt,lt+l*cos(theta),l*sin(theta), 0,l*cos(theta),ht+l*sin(theta),
////						       	          wt,lt,0, 2*wt,0,ht, 2*wt,l*cos(theta),ht+l*sin(theta),
////						    	          wt,lt,0, 2*wt,l*cos(theta),ht+l*sin(theta), wt,lt+l*cos(theta),l*sin(theta)};
//				//Rotate shape
//				for (int i = 0; i<3*numCells; i++){
//					double x = m[i*3]; //extract x,y,z values
//					double y = m[i*3+1];
//					double z = m[i*3+2];
//					x-=wt; //center shape at origin cubeFootprint_xmaxStructure_zmax
//				    y-=(lt+l*cos(theta))/2.0;
//				    z-=(ht+l*sin(theta))/2.0;
//					double rx = x*cos(rho)-y*sin(rho); //rotate shape counterclockwise w/respect to south
//					double ry = x*sin(rho)+y*cos(rho);
//					rx+=xBound/2.0; //center shape at 5,5,5
//				    ry+=yBound/2.0;
//				    z+=zBound/2.0;
//				    m[3*i] = rx;
//				    m[3*i+1] = ry;
//				    m[3*i+2] = z;
//				}
//
//				modelPointer = m;
//
////				//print coordinates of the structure
////				for(int i=0; i<numCells; i++) {
////					//print coordinates of the ith triangular cell
////					for(int j=0; j<9; j++) {
////						std::cout << m[9*i+j];
////						if(j%3 == 2) {
////							std::cout << " ";
////						} else {
////							std::cout << ",";
////						}
////					}
////					std::cout << std::endl;
////				}
////				std::cout << std::endl;
//
//				int v_index_triflat[2][3] = {{0,1,2},{0,2,1}};  // stores integer indices for 3 verts of each triangle
//
//				Triangle *temptri;
//				double pn[3];
//				double newvertex1[3];
//				double newvertex2[3];
//				double newvertex3[3];
//
//				for(int i=0; i<numCells; i++) {
//					Cell_3 *cel = &cellArray[i];
//
//					cel->numVertices = 3;
//					cel->numTriangles = 2; // like TRIFLAT, it has a top and bottom triangle
//
//
//					for (int nt=0; nt < cel->numTriangles; nt++)
//					{
//						cel->triangleList[nt] = new Triangle;
//					}
//
//					cel->vertices[0] = new Position;
//					cel->vertices[0]->x = modelPointer[i*9+0]; // get coordinate from global pointer to array modelCoordinates
//					cel->vertices[0]->y = modelPointer[i*9+1];
//					cel->vertices[0]->z = modelPointer[i*9+2];
//
//					cel->vertices[1] = new Position;
//					cel->vertices[1]->x = modelPointer[i*9+3];
//					cel->vertices[1]->y = modelPointer[i*9+4];
//					cel->vertices[1]->z = modelPointer[i*9+5];
//
//					cel->vertices[2] = new Position;
//					cel->vertices[2]->x = modelPointer[i*9+6];
//					cel->vertices[2]->y = modelPointer[i*9+7];
//					cel->vertices[2]->z = modelPointer[i*9+8];
//
//
//					// for former complicated shapes there could be many more triangles. but a flat cell has two sides-->numTriangles=2
//					for (int m=0; m < (cel->numTriangles); ++m)
//					{
//						temptri = cel->triangleList[m];
//						temptri->index[0] = v_index_triflat[m][0];
//						temptri->index[1] = v_index_triflat[m][1];
//						temptri->index[2] = v_index_triflat[m][2];
//
//						//cal the normal to a particular face
//
//						pn[0]=0; pn[1]=0; pn[2]=0;
//						newvertex1[0] = (cel->vertices[temptri->index[0]])->x;
//						newvertex1[1] = (cel->vertices[temptri->index[0]])->y;
//						newvertex1[2] = (cel->vertices[temptri->index[0]])->z;
//						newvertex2[0] = (cel->vertices[temptri->index[1]])->x;
//						newvertex2[1] = (cel->vertices[temptri->index[1]])->y;
//						newvertex2[2] = (cel->vertices[temptri->index[1]])->z;
//						newvertex3[0] = (cel->vertices[temptri->index[2]])->x;
//						newvertex3[1] = (cel->vertices[temptri->index[2]])->y;
//						newvertex3[2] = (cel->vertices[temptri->index[2]])->z;
//						plane_exp_normal_3d(newvertex1, newvertex2, newvertex3, pn);
//						temptri->Normal[0]=pn[0]; temptri->Normal[1]=pn[1]; temptri->Normal[2]=pn[2];
//
//
//						temptri->vertList[0] = newvertex1[0]; temptri->vertList[1] = newvertex1[1]; temptri->vertList[2]=newvertex1[2];
//						temptri->vertList[3] = newvertex2[0]; temptri->vertList[4] = newvertex2[1]; temptri->vertList[5]=newvertex2[2];
//						temptri->vertList[6] = newvertex3[0]; temptri->vertList[7] = newvertex3[1]; temptri->vertList[8]=newvertex3[2];
//
//						temptri->Area = triangle_area_3d(temptri->vertList); // call to geometry.cpp function to get triangles area
//
//
//						if (m==0)
//						{
//							totalSurfaceArea+=temptri->Area;  // get the total area
//							temptri->is_positive=true;
//						} // finally, add area to the total structure area (m==0 => postive norm side only)
//						if (m==1) {
//							temptri->is_positive=false;
//						}
//					}
//
//				}
//
//				double day_energee=0;
//				day_energee = newCalculatePowerDayCycle(cellArray, modelPointer);   // calculate power at several instants in the day
//
//				std::cout << setiosflags(std::ios::fixed) << std::setprecision(6);
//				std::cout << "rho: " << radtodeg(rho) << ", theta: " << radtodeg(theta) << ", phi: " << radtodeg(phi)
//						<< ", energy: " << day_energee << "\n";
//
////				myfile << setiosflags(std::ios::fixed) << std::setprecision(6);
////				myfile << radtodeg(rho) << " " << radtodeg(theta) << " " << radtodeg(phi) << " " << day_energee << "\n";
//			}
//		}
//	}
//
////	myfile.close();
//
//	return 0;
//}

//int main(int argc, char *argv[]) {
//	strcpy(Input_CMD, argv[1]);
//
//	getcwd(Main_Dir, 255); //initializes the main directory
//	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	strcat(Clone_Dir, argv[2]);
//	strcpy(Unique_Dir, Clone_Dir);
//	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
//
//	strcpy(Clone_Dir, Main_Dir);
//
//	ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
//
//	Get2DTriangleGrid(); // precalculate grid points
//
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//
//
//
//
//	numCoordinates=ReadStructureVertices_Count(0);
//
//	double modelCoordinates[numCoordinates];
//
//	modelPointer = modelCoordinates;
//
//	ReadStructureVertices_Retrieve(modelCoordinates, 0);
//
//	Cell_3 cellArray[numCells];
//
//	BuildStructure(cellArray, MODEL);
//
//	double day_energee=0;
//
//	day_energee = newCalculatePowerDayCycle(cellArray, modelPointer);   // calculate power at several instants in the day
//
//	return 0;
//}

//###############################################################################
//###############################################################################
//###############################################################################




//xBound, yBound and zBound must not have more than numSigFigs (default 6) digits after the decimal point
void gen_random_structure(double *modelCoordinates, int numCells, double xBound, double yBound, double zBound){
	int xIntBound = xBound*numSigFigs;
	int yIntBound = yBound*numSigFigs;
	int zIntBound = zBound*numSigFigs;

	double x, y, z;

	for(int i = 0; i < numCells; i++) {
		for(int j = 0; j < 3; j++) {
			x = (rand() % (xIntBound+1) ) / numSigFigs;
			y = (rand() % (yIntBound+1) ) / numSigFigs;
			z = (rand() % (zIntBound+1) ) / numSigFigs;

			modelCoordinates[9*i + 3*j] 	= x;
			modelCoordinates[9*i + 3*j + 1] = y;
			modelCoordinates[9*i + 3*j + 2] = z;
		}
	}
}






//// MAIN LOOP

int main(int argc, char *argv[])
{
	// depending on the first argument,
	// "sim" => regular simulation. followed by arguments "input_file" and "unique output directory"
	// or "view" => just view the outputs from an already completed simulation.
		// just follwed by argument "unique output directory" which is holding the stuff you want to look at

	string mainInputFile(argv[1]);
	string outDirectory = argv[2];
	int mode = atoi(argv[3]);

	strcpy(Input_CMD, argv[1]);
	getcwd(Main_Dir, 255); //initializes the main directory
	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
	strcpy(Clone_Dir, Main_Dir);
	strcat(Clone_Dir, argv[2]);
	strcpy(Unique_Dir, Clone_Dir);
	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
	strcpy(Clone_Dir, Main_Dir);

	ReadMainInputFile(mainInputFile, mode);

	if(mode == 0) { //GA
		if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points

		printf("Index of refraction: %f\n", index_ref_plastic);
		printf("latitude = %f\n",  latitude);

/*JIN: there is a better way to delete old data in C++*/
//		/////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
//
//		strcpy(Clone_Dir, Unique_Dir);
//		char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
//		FILE *fileBestF= fopen(output_catab, "w");
//		fclose(fileBestF);

		////////////////////////////////////////////////////////////////////////////////////////////////////

		// 8/04/09: set all this qwikStore stuff to 0...to start out
		for (int qi0=0; qi0<qwikIndMax; qi0++) {
			qwikStoreEnergies[qi0]=0;
			qwikStoreReflEnergies[qi0]=0;
			for (int qi1=0; qi1<qwikCoordMax; qi1++) {
				qwikStoreCoords[qi0][qi1]=0;
			}
		}

		if (GA_willload_population==1)
		{
			printf("Loading initial population from model files\n");
			numCoordinates=ReadStructureVertices_Count(0);
			printf("number of coord: %i\n", numCoordinates);
			double modelCoordinates[numCoordinates];
			modelPointer=modelCoordinates;
			strcpy(Clone_Dir, Unique_Dir);
			char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);

			FILE *inputPop = fopen(myGAinput, "w");
			fprintf(inputPop, "%i\n", GA_population);

			printf("%i distinct models\n",GA_numDistinctModels);
			for (int nm=0; nm < GA_numDistinctModels; nm++)
			{
				printf("first %i\n", nm);
				ReadStructureVertices_Retrieve(modelCoordinates, nm);
				printf("2nd %i\n", nm);

				for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
				{
					for (int nc=0; nc < numCells; nc++)
					{
						// format: tab \t between all dec vars. \n after objective value written, before next individal
						for (int np=0; np < 9; np++) {
							fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);
						}
						//printf("%i\n",nc);
					}
					if (GA_numberof_objectives==2)
						{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
					}
					if (GA_numberof_objectives==1) {
						fprintf(inputPop, "0\n");  // newline then next individual
					}
				}
				// close just to have one file open at a time (ReadStructure... opens a file)
			}
			fclose(inputPop);
		} else {
			printf("Generating random starting population within coordinate bounds: models in input file not used...\n");

			numCoordinates = 9*numCells;
			printf("number of coord: %i\n", numCoordinates);
			double modelCoordinates[numCoordinates]; //array of coordinates in the cell
			modelPointer=modelCoordinates;

			ofstream inputPop;
			inputPop.open("test_output/initialPopulation.txt");

			inputPop << GA_population << endl;
			printf("%i distinct models\n",GA_numDistinctModels);


			/* initialize random seed: */
			srand (time(NULL));

			for (int i=0; i < GA_numDistinctModels; i++)
			{
				//generate a new random structure for each distinct model
				gen_random_structure(modelCoordinates, numCells, xBound, yBound, zBound);

				for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
				{
					for (int j=0; j < numCells; j++)
					{
						// format: tab \t between all dec vars. \n after objective value written, before next individal
						for (int k=0; k < 9; k++) {
							inputPop << modelCoordinates[9*j + k] << "\t";
						}
					}

					if (GA_numberof_objectives==2) {
						inputPop << "0" << "\t";
						inputPop << "0" << endl;
					}

					if (GA_numberof_objectives==1) {
						inputPop << "0" << endl;
					}
				}

			}
			inputPop.close();
		}

		printf("Assuming a plastic solar cell efficiency of %f \n", cellEfficiency);
		printf("Starting simulation for %i generations.\n", GA_maxgenerations);
		WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
		WriteMainOutputFile();
		char mcommand[200];
		sprintf(mcommand, "mkdir ");
		strcpy(Clone_Dir, Unique_Dir);
		strcat(mcommand, Clone_Dir);
		strcat(mcommand, "models");
		system(mcommand);
		mainGA();

		Ggen_show=GA_gen_current; // OPENGL stuff follow this....
		// display (also print to file) desired output



	} else if(mode == 1) { //fixed shape
		Get2DTriangleGrid(); // precalculate grid points

		// 8/04/09: set all this qwikStore stuff to 0...to start out
		for (int qi0=0; qi0<qwikIndMax; qi0++) {
			qwikStoreEnergies[qi0]=0;
			qwikStoreReflEnergies[qi0]=0;
			for (int qi1=0; qi1<qwikCoordMax; qi1++) {
				qwikStoreCoords[qi0][qi1]=0;
			}
		}

		numCoordinates=ReadStructureVertices_Count(0);
		double modelCoordinates[numCoordinates];
		modelPointer = modelCoordinates;
		ReadStructureVertices_Retrieve(modelCoordinates, 0);
		Cell_3 cellArray[numCells];
		BuildStructure(cellArray, MODEL);

		double efficiencies[numCells], reflIndices[numCells];

		for(int i=0; i<numCells; i++) {
			efficiencies[i] = cellEfficiency;
			reflIndices[i] = index_ref_plastic;
		}

		double day_energee=0;
		day_energee = newCalculatePowerDayCycle(cellArray, modelPointer, efficiencies, reflIndices, true, true, outDirectory);   // calculate power at several instants in the day


	} else if(mode == 2) { //fixed shape with cell efficiency and reflectivity specifications

		Get2DTriangleGrid(); // precalculate grid points

		// 8/04/09: set all this qwikStore stuff to 0...to start out
		for (int qi0=0; qi0<qwikIndMax; qi0++) {
			qwikStoreEnergies[qi0]=0;
			qwikStoreReflEnergies[qi0]=0;
			for (int qi1=0; qi1<qwikCoordMax; qi1++) {
				qwikStoreCoords[qi0][qi1]=0;
			}
		}

		numCoordinates=ReadStructureVertices_Count(0);
		double modelCoordinates[numCoordinates];
		modelPointer = modelCoordinates;
		ReadStructureVertices_Retrieve(modelCoordinates, 0);
		Cell_3 cellArray[numCells];
		BuildStructure(cellArray, MODEL);

		double efficiencies[numCells], reflIndices[numCells];

		ifstream ifs(modelSpecDir.data());

		for(int i=0; i<numCells; i++) {
			ifs >> efficiencies[i];
			ifs >> reflIndices[i];
		}

		ifs.close();

		double day_energee=0;
		day_energee = newCalculatePowerDayCycle(cellArray, modelPointer, efficiencies, reflIndices, true, true, outDirectory);   // calculate power at several instants in the day

	} else if(mode == 3) { //Monte Carlo

		Get2DTriangleGrid(); // precalculate grid points

		// 8/04/09: set all this qwikStore stuff to 0...to start out
		for (int qi0=0; qi0<qwikIndMax; qi0++) {
			qwikStoreEnergies[qi0]=0;
			qwikStoreReflEnergies[qi0]=0;
			for (int qi1=0; qi1<qwikCoordMax; qi1++) {
				qwikStoreCoords[qi0][qi1]=0;
			}
		}

		numCoordinates=ReadStructureVertices_Count(0);
		double modelCoordinates[numCoordinates];
		modelPointer = modelCoordinates;
		ReadStructureVertices_Retrieve(modelCoordinates, 0); //0 because only a single fixed run
		Cell_3 cellArray[numCells];
		BuildStructure(cellArray, MODEL);

		double efficiencies[numCells], reflIndices[numCells];

		ifstream ifs(modelSpecDir.data());

		for(int i=0; i<numCells; i++) {
			ifs >> efficiencies[i];
			ifs >> reflIndices[i];
		}

		ifs.close();

		double day_energee=0;
		day_energee = newCalculatePowerDayCycle(cellArray, modelPointer, efficiencies, reflIndices, true, true, outDirectory);

	}


	return 0;


//	if (strcmp(argv[1],"sim")==0)
//	{
//	printf("damnit\n");
//	// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
//	strcpy(Input_CMD, argv[2]);
//	getcwd(Main_Dir, 255);
//	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	strcat(Clone_Dir, argv[3]);
//	strcpy(Unique_Dir, Clone_Dir);
//	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
//	/*cout << Main_Dir << endl;
//	cout << Unique_Dir << endl;
//	cout << Input_CMD << endl;*/
//	//printf("%s\n%s\n%s\n", Main_Dir, Unique_Dir, Input_CMD);
//	strcpy(Clone_Dir, Main_Dir);
//	ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
////	cout << "got it" << endl;
//	if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
//	printf("Index of refraction: %f\n", index_ref_plastic);
//
//	// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
//	doSingleReflections = atoi(argv[4]);
//
//	/*for (int aaa=0; aaa<=9; aaa++) {
//		for (int bbb=0; bbb<=9; bbb++) {
//			for (int ccc=0; ccc<=9; ccc++) {
//				for (int ddd=0; ddd<=9; ddd++) {
//					for (int eee=0; eee<=9; eee++) {
//						for (int fff=0; fff<=9; fff++) {
//					printf("%i%i%i%i%i%i\n", aaa,bbb,ccc,ddd,eee,fff);}}}}}}*/
//
//	if (input_i==0)
//	{
//		// if not an example case
//	}
//	else
//	{
//		longitude   = INPUT_examples[input_i][0];   // longitude of location
//		if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
//			latitude    = INPUT_examples[input_i][1];   // latitude of location
//
//			if (globalInputLatitude==-401) {globalCellsOneSided=true;}
//		}
//		else {
//			latitude = globalInputLatitude;
//		}
//		avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//		avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//		elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//		timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//	}
//	printf("latitude = %f\n",  latitude);
//	//printf("Day %i, Month %i, year %i\n", Input_startday, Input_startmonth, Input_startyear);
//	//printf("longitude %f latitude %f\n", longitude, latitude);
//
//	/////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
//
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
//	FILE *fileBestF= fopen(output_catab, "w");
//	fclose(fileBestF);
//
//	////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	if (GA_willload_population==1)
//	{
//		printf("Loading initial population from model files\n");
//		numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
//		printf("number of coord: %i\n", numCoordinates);
//		double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
//	//	double allModelDecisionVars[numCoordinates/*6/9*/GA_numDistinctModels]; // holds all decision vars for all models
//		modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
//	//	allDecVarPointer=allModelDecisionVars;                   // global pointer
//		/*FILE *inputPop = fopen(GA_initial_population_file,"w");
//		fclose (inputPop);*/
////		cout << Unique_Dir <<endl;
//		strcpy(Clone_Dir, Unique_Dir);
//		char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
//		//printf("%s\n", myGAinput);
////		cout << myGAinput << endl;
//		FILE *inputPop = fopen(myGAinput, "w"); // write the population file
//		fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
//		printf("%i distinct models\n",GA_numDistinctModels);
//		for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
//		{
//			printf("first %i\n", nm);
//			ReadStructureVertices_Retrieve(modelCoordinates, nm);
//			printf("2nd %i\n", nm);
//			for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
//			{
//				for (int nc=0; nc < numCells; nc++)
//				{
//					// format: tab \t between all dec vars. \n after objective value written, before next individal
//					for (int np=0; np < 9; np++) {
//						fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);
//					}
//					//printf("%i\n",nc);
//				}
//				if (GA_numberof_objectives==2)
//					{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
//				}
//				if (GA_numberof_objectives==1) {
//					fprintf(inputPop, "0\n");  // newline then next individual
//				}
//			}
//			// close just to have one file open at a time (ReadStructure... opens a file)
//		}
//		fclose(inputPop);
//	}
//	else
//	{
//		printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
//	}
//
//	printf("Assuming a plastic solar cell efficiency of 0.06\n");
//	printf("Starting simulation for %i generations.\n", GA_maxgenerations);
//	WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
//	WriteMainOutputFile();
//	char mcommand[200];
//	sprintf(mcommand, "mkdir ");
//	strcpy(Clone_Dir, Unique_Dir);
//	strcat(mcommand, Clone_Dir);
//	strcat(mcommand, "models");
//	system(mcommand);
//	mainGA();
//
//	Ggen_show=GA_gen_current; // OPENGL stuff follow this....
//	// display (also print to file) desired output
//	} // end SIMULATION option mainloop
//
//	if (strcmp(argv[1], "grid")==0) // expects two arguments: 1) output directory from previous run, 2) # grids (will be squared)
//	{
//		getcwd(Main_Dir, 255);
//		strcpy(Clone_Dir, Main_Dir);
//		strcpy(Main_Dir, ""); 		// make it blank just to avoid error later in model function
//		strcat(Clone_Dir, "/");
//		strcat(Clone_Dir, argv[2]);
//		printf("Searching in simulation output directory: %s\n", Clone_Dir);
//		strcpy(Unique_Dir, Clone_Dir);
//		strcat(Unique_Dir, "/");
//
//		// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
//		doSingleReflections = atoi(argv[4]);
//
//		globalInputLatitude=-400; //added 8/13/09...temporary fix, assume im always "gridding" at berkeley for now
//		if (input_i==0)
//		{
//			// if not an example case
//		}
//		else
//		{
//			longitude   = INPUT_examples[input_i][0];   // longitude of location
//			if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
//				latitude    = INPUT_examples[input_i][1];   // latitude of location
//
//				if (globalInputLatitude==-401) {globalCellsOneSided=true;}
//			}
//			else {
//				latitude = globalInputLatitude;
//			}
//			avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//			avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//			elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//			timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//		}
//		printf("latitude = %f\n",  latitude);
//
//		strcpy(Clone_Dir, Unique_Dir);
//		strcat(Clone_Dir, OutMain_name);
//		ReadMainOutputFile( Clone_Dir);			// 1) get some typical variables, e.g., like # generations, if needed
//		raysPerCell=atoi(argv[3]); 			// whatever this is, the number is squared in next function
//		if (use_trigrid==1) {Get2DTriangleGrid();} 	// precalculate grid points
//
//		///////////// create new version of these  files//////////////
//		/*strcpy(Clone_Dir, Unique_Dir);
//		char *output_catsc = strcat(Clone_Dir, OutPowerPerSubcell);
//		FILE *fileSubPOW= fopen(output_catsc, "w");
//		fclose(fileSubPOW);*/	// commented out 5/13/09 when added split by day iteration and "w" each one just once...no "a"
//
//		strcpy(Clone_Dir, Unique_Dir);
//		char *output_catgc = strcat(Clone_Dir, OutGridCoordinates);
//		FILE *fileGridCoord= fopen(output_catgc, "w");
//		fclose(fileGridCoord);
//		/////////////////////////////////////////////////////////
//		strcpy(Clone_Dir, Unique_Dir);
//		//char *GAmodel_prefix = "Model_gen";
//		char *GAmodel_suffix = "_cells.txt";
//		strcat(Clone_Dir, OutModel_prefix);
//		char  gnstring[200];
//		sprintf(gnstring,"%i", (GA_maxgenerations-1) ); 	// generation value was obtained by reading main output file
//		strcat(Clone_Dir, gnstring);
//		strcat(Clone_Dir, GAmodel_suffix);
//		strcpy(Input_modelALL[0], Clone_Dir);				// finally, full location of the model..
//		numCoordinates=ReadStructureVertices_Count(0);
//		double modelCoordinates[numCoordinates];
//		ReadStructureVertices_Retrieve(modelCoordinates, 0);
//
//		// I could just as well have obtained the coordinates from the AllBestGen file, but I think getting
//		// the coordinates from a model will be more versatile...or the same...
//		getSubCellData = true;
//		globalSetup = new GlobalSetup;			// just to make globalEvaluate happy (actually useless 'if' branch...)
//		globalHaveDonePreEval=1;			// fake this, since no GA is being started, just want one iteration
//		Write_tsampledata=true;		// to yes, write the coordinates of gridpoints first time
//		globalSetup->gaType = SGA;
//		globalSetup->finalNoOfObjectives=1;
//		globalSetup->noOfDecisionVariables=numCells*9;
//		double dummyObjective[5];
//		double dummyDouble[5];
//		int dummyInteger[5];
//		globalEvaluate(modelCoordinates, dummyObjective, dummyDouble, dummyDouble, dummyInteger);
//		DisplayPowerOutputGridVersion(); // output a power vs. frac-hour data set with #gridpoints appended to filename
//
//		// get & write some angle analysis data. distribution weighted according to sum of energy over the day
//		GetZenithAngleDistribution(modelCoordinates);
//		GetAzimuthAngleDistribution(modelCoordinates);
//		Get2dAngleDistribution(modelCoordinates);	// for 2dimensional-bins (x,y) distribution
//		WriteAngleDistrFile();
//	}
//
//	/////////////STRETCH WAIST STUFF!!!!!!!!!!!!!!1
//	if (strcmp(argv[1],"hf")==0)
//	{
//	// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
//	strcpy(Input_CMD, argv[2]);
//	getcwd(Main_Dir, 255);
//	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	strcat(Clone_Dir, argv[3]);
//	strcpy(Unique_Dir, Clone_Dir);
//	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
//	if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
//
//	// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
//	doSingleReflections = atoi(argv[4]);
//
//	if (input_i==0)
//	{}
//	else
//	{
//		longitude   = INPUT_examples[input_i][0];   // longitude of location
//		if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
//			latitude    = INPUT_examples[input_i][1];   // latitude of location
//
//			if (globalInputLatitude==-401) {globalCellsOneSided=true;}
//		}
//		else {
//			latitude = globalInputLatitude;
//		}
//		avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//		avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//		elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//		timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//	}
//	printf("latitude = %f\n",  latitude);
//
//	/////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
//	FILE *fileBestF= fopen(output_catab, "w");
//	fclose(fileBestF);
//	////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	double globalFunCoordinates[1000];
//	double globalBoxCoordinates[1000];
//
//	if (GA_willload_population==1)
//	{
//		printf("Loading initial population from model files\n");
//		numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
//		printf("number of coord: %i\n", numCoordinates);
//		double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
//		modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
//		strcpy(Clone_Dir, Unique_Dir);
//		char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
//		FILE *inputPop = fopen(myGAinput, "w"); // write the population file
//		fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
//		for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
//		{
//			ReadStructureVertices_Retrieve(modelCoordinates, nm);
//
//			for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
//			{
//				for (int nc=0; nc < numCells; nc++)
//				{
//					// format: tab \t between all dec vars. \n after objective value written, before next individal
//					for (int np=0; np < 9; np++) {
//						fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);
//
//						if (nm==0) {
//							globalFunCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
//						if (nm==1) {
//							globalBoxCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
//
//					}
//				}
//				if (GA_numberof_objectives==2)
//					{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
//				}
//				if (GA_numberof_objectives==1) {
//					fprintf(inputPop, "0\n");  // newline then next individual
//				}
//			}
//			// close just to have one file open at a time (ReadStructure... opens a file)
//		}
//		fclose(inputPop);
//	}
//	else
//	{
//		printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
//	}
//
//	printf("Assuming a plastic solar cell efficiency of 0.06\n");
//	//printf("Starting simulation for %i generations.\n", GA_maxgenerations);
//	WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
//	WriteMainOutputFile();
//	char mcommand[200];
//	sprintf(mcommand, "mkdir ");
//	strcpy(Clone_Dir, Unique_Dir);
//	strcat(mcommand, Clone_Dir);
//	strcat(mcommand, "models");
//	system(mcommand);
//
//	double heighttest_BoxPower[201]; // box power at each height
//	double heighttest_FunPower[201]; // funnel power at each height
//	for (int jij=0; jij<201; jij++) {
//		heighttest_BoxPower[jij]=0;
//		heighttest_FunPower[jij]=0;
//	}
//
//	int maxHeightGoto=100;
//	for (int ghi=2; ghi <= maxHeightGoto; ghi+=2) {
//		printf("\nheight=%i\n",ghi);
//		AlterBoxCoordinates(globalBoxCoordinates, ghi);
//		Cell_3::cell_counter = 0;
//		Cell_3 cellArray[numCells];
//		modelPointer = globalBoxCoordinates;
//		BuildStructure(cellArray, MODEL);
//		heighttest_BoxPower[ghi]=newCalculatePowerDayCycle(cellArray, globalBoxCoordinates);
//
//		AlterFunnelCoordinates(globalFunCoordinates, ghi); // takes coordinates and height as arguments
//		printf("%f\n",globalFunCoordinates[2]);
//		Cell_3::cell_counter = 0;
//		Cell_3 cellArray2[numCells];
//		modelPointer = globalFunCoordinates;
//		BuildStructure(cellArray2, MODEL);
//		heighttest_FunPower[ghi]=newCalculatePowerDayCycle(cellArray2, globalFunCoordinates);
//	}
//
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_concat1 = strcat(Clone_Dir, "HeightBoxFunnelCompare");
//	FILE *hOUTfile = fopen(output_concat1, "w");
//
//	fprintf(hOUTfile, "Height\t"); fprintf(hOUTfile, "Box (kWh)\t"); fprintf(hOUTfile, "Funnel (kWh)\t");
//	fprintf(hOUTfile, "percent diff\n");
//	for (int cff=2; cff<= maxHeightGoto; cff+=2)
//	{
//		double boxkWh = heighttest_BoxPower[cff]/(3600*1000);
//		double funkWh = heighttest_FunPower[cff]/(3600*1000);
//		double perC = ( 100*fabs(heighttest_FunPower[cff]-heighttest_BoxPower[cff])/heighttest_BoxPower[cff] );
//		fprintf(hOUTfile, "%i\t%f\t%f\t%f\n", cff, boxkWh, funkWh, perC);
//	}
//	fclose(hOUTfile);
//
//	Ggen_show=GA_gen_current; // OPENGL stuff follow this....
//	// display (also print to file) desired output
//	} // end SIMULATION option mainloop
//
//
//	////////?WAIST FUNNEL STUFF////////////////
//	if (strcmp(argv[1],"wf")==0)
//	{
//	// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
//	strcpy(Input_CMD, argv[2]);
//	getcwd(Main_Dir, 255);
//	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	strcat(Clone_Dir, argv[3]);
//	strcpy(Unique_Dir, Clone_Dir);
//	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
//	if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
//
//	// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
//	doSingleReflections = atoi(argv[4]);
//
//	if (input_i==0)
//	{}
//	else
//	{
//		longitude   = INPUT_examples[input_i][0];   // longitude of location
//		if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
//			latitude    = INPUT_examples[input_i][1];   // latitude of location
//
//			if (globalInputLatitude==-401) {globalCellsOneSided=true;}
//		}
//		else {
//			latitude = globalInputLatitude;
//		}
//		avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//		avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//		elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//		timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//	}
//	printf("latitude = %f\n",  latitude);
//
//	/////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
//	FILE *fileBestF= fopen(output_catab, "w");
//	fclose(fileBestF);
//	////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	double globalFunCoordinates[1000];
//	//double globalBoxCoordinates[1000];
//
//	if (GA_willload_population==1)
//	{
//		printf("Loading initial population from model files\n");
//		numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
//		printf("number of coord: %i\n", numCoordinates);
//		double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
//		modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
//		strcpy(Clone_Dir, Unique_Dir);
//		char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
//		FILE *inputPop = fopen(myGAinput, "w"); // write the population file
//		fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
//		for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
//		{
//			ReadStructureVertices_Retrieve(modelCoordinates, nm);
//
//			for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
//			{
//				for (int nc=0; nc < numCells; nc++)
//				{
//					// format: tab \t between all dec vars. \n after objective value written, before next individal
//					for (int np=0; np < 9; np++) {
//						fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);
//
//						if (nm==0) {
//							globalFunCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
//
//					}
//				}
//				if (GA_numberof_objectives==2)
//					{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
//				}
//				if (GA_numberof_objectives==1) {
//					fprintf(inputPop, "0\n");  // newline then next individual
//				}
//			}
//			// close just to have one file open at a time (ReadStructure... opens a file)
//		}
//		fclose(inputPop);
//	}
//	else
//	{
//		printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
//	}
//
//	printf("Assuming a plastic solar cell efficiency of 0.06\n");
//	//printf("Starting simulation for %i generations.\n", GA_maxgenerations);
//	WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
//	WriteMainOutputFile();
//	char mcommand[200];
//	sprintf(mcommand, "mkdir ");
//	strcpy(Clone_Dir, Unique_Dir);
//	strcat(mcommand, Clone_Dir);
//	strcat(mcommand, "models");
//	system(mcommand);
//
//	double heighttest_FunPower[50][50]; // funnel power at each height
//
//	int maxYGoto=9;
//	int maxZGoto=21;
//	for (int fiy=1; fiy <= maxYGoto; fiy+=1) {
//		for (int fiz=1; fiz <= maxZGoto; fiz+=1) {
//			printf("\nwaist_y=%i, waist_z=%i\n",fiy,fiz);
//
//			AlterFunnelWaist(globalFunCoordinates, fiy, fiz); // takes coordinates and height as arguments
//			//printf("%f\n",globalFunCoordinates[2]);
//			Cell_3::cell_counter = 0;
//			Cell_3 cellArray2[numCells];
//			modelPointer = globalFunCoordinates;
//			BuildStructure(cellArray2, MODEL);
//			heighttest_FunPower[fiy][fiz]=newCalculatePowerDayCycle(cellArray2, globalFunCoordinates);
//		}
//	}
//
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_concat1 = strcat(Clone_Dir, "WaistFunnelCompare");
//	FILE *hOUTfile = fopen(output_concat1, "w");
//
//	fprintf(hOUTfile, "y\t"); fprintf(hOUTfile, "z\t"); fprintf(hOUTfile, "Energy (kWh)\n");
//	for (int fiy=1; fiy <= maxYGoto; fiy+=1) {
//		for (int fiz=1; fiz <= maxZGoto; fiz+=1) {
//			double funkWh = heighttest_FunPower[fiy][fiz]/(3600*1000);
//			fprintf(hOUTfile, "%i\t%i\t%f\n", fiy, fiz, funkWh);
//		}
//	}
//	fclose(hOUTfile);
//
//	Ggen_show=GA_gen_current; // OPENGL stuff follow this....
//	// display (also print to file) desired output
//	} // end SIMULATION option mainloop
//
//	//////////////////////////////////////////////////////////////////
//	////////?ROOF ARRANGEMENT SUMULATION OPTIMIZATION ////////////////
//	//////////////////////////////////////////////////////////////////
//	if (strcmp(argv[1],"rf")==0)
//	{
//	// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
//	strcpy(Input_CMD, argv[2]);
//	getcwd(Main_Dir, 255);
//	strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	strcat(Clone_Dir, argv[3]);
//	strcpy(Unique_Dir, Clone_Dir);
//	strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
//	strcpy(Clone_Dir, Main_Dir);
//	ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
//	if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
//
//	// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
//	// argument 4 = roofx dimension
//	// argument 5 = roofy dimension
//	// argument 6 = # of instances of unit model to scatter around
//	// argument 7 = reflections turned on
//	//double roof_xmin = 0;		// bounds of the roof on x,y plane
//	//double roof_ymin = 0;
//	double roof_xmax = atof(argv[4]);
//	double roof_ymax = atof(argv[5]);
//	int InstCopies_numberof = atoi(argv[6]);
//
//	doSingleReflections = atoi(argv[7]);
//
//	if (input_i==0)
//	{}
//	else
//	{
//		longitude   = INPUT_examples[input_i][0];   // longitude of location
//		if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
//			latitude    = INPUT_examples[input_i][1];   // latitude of location
//
//			if (globalInputLatitude==-401) {globalCellsOneSided=true;}
//		}
//		else {
//			latitude = globalInputLatitude;
//		}
//		avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
//		avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
//		elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
//		timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
//	}
//	printf("latitude = %f\n",  latitude);
//
//	/////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
//	strcpy(Clone_Dir, Unique_Dir);
//	char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
//	FILE *fileBestF= fopen(output_catab, "w");
//	fclose(fileBestF);
//	////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	// 8/04/09: set all this qwikStore stuff to 0...to start out
//	for (int qi0=0; qi0<qwikIndMax; qi0++) {
//		qwikStoreEnergies[qi0]=0;
//		qwikStoreReflEnergies[qi0]=0;
//		for (int qi1=0; qi1<qwikCoordMax; qi1++) {
//			qwikStoreCoords[qi0][qi1]=0;
//		}
//	}
//
//	double globalBasicCoordinates[1000]; // basic coordinates of unit structure (at x,y = 0,0 corner)
//	//double globalBoxCoordinates[1000];
//
//	if (GA_willload_population==1)
//	{
//		printf("Loading initial population from model files\n");
//		numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
//		printf("number of coord: %i\n", numCoordinates);
//		double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
//		modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
//		strcpy(Clone_Dir, Unique_Dir);
//		char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
//		FILE *inputPop = fopen(myGAinput, "w"); // write the population file
//		fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
//		for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
//		{
//
//			ReadStructureVertices_Retrieve(modelCoordinates, nm);
//
//			for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
//			{
//				for (int nc=0; nc < numCells; nc++)
//				{
//					// format: tab \t between all dec vars. \n after objective value written, before next individal
//					for (int np=0; np < 9; np++) {
//						fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);
//
//						if (nm==0) {
//							globalBasicCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
//
//					}
//				}
//				if (GA_numberof_objectives==2)
//					{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
//				}
//				if (GA_numberof_objectives==1) {
//					fprintf(inputPop, "0\n");  // newline then next individual
//				}
//			}
//			// close just to have one file open at a time (ReadStructure... opens a file)
//		}
//		fclose(inputPop);
//	}
//	else
//	{
//		printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
//	}
//	printf("check1\n");
//
//	int simORmodel = 1; // 1 for just print out best model
//
//	double wholeRoofModel[InstCopies_numberof*numCells*9];
//	int copyCells_numberof = numCells;
//	int num_iterations = GA_maxgenerations; // may as well use same number
//
//	if (simORmodel==1)
//	{
//		double got_positions/*[num_iterations+100]*/[InstCopies_numberof*2];
//		int Best_config = atoi(argv[7]);// just for simplicity, use the singleDoreflections input for this as any integer
//		printf("check2\n");
//		ReadRoofPositions(got_positions, InstCopies_numberof*2, Best_config);
//		printf("check3\n");
//
//		AlterRoofStructure(wholeRoofModel, globalBasicCoordinates, got_positions, InstCopies_numberof, numCells);
//
//		strcpy(Clone_Dir, Unique_Dir);
//
//		char *output_catM = strcat(Clone_Dir, "BestFullModel.txt");
//		//strcat(output_cat, itoa(whichgen));
//		/*char istring[200];
//		sprintf(istring,"%i",whichgen);
//		strcat(output_cat, istring);
//		strcat(output_cat, "_cells.txt");*/
//		//strcat(output_cat, itoa(numCells);
//		FILE *OUTfileM = fopen(output_catM, "w");
//
//		for (int cff=0; cff<InstCopies_numberof; cff++)
//		{
//			for (int ccc=0; ccc<copyCells_numberof; ccc++)
//			{
//				for (int cdd=0; cdd<9; cdd++)
//				{
//					//fprintf(OUTfile, "%f", globalParamsAvgBest[whichgen*numCells*9+cff*9+cdd]);
//					fprintf(OUTfileM, "%f", wholeRoofModel[cff*copyCells_numberof*9 + ccc*9 + cdd]);
//					if (cdd<8) {fprintf(OUTfileM, ",");}
//					if (cdd==8) {fprintf(OUTfileM, "\n");}
//				}
//			}
//		}
//
//		fclose(OUTfileM);
//	}
//
//	if (simORmodel==0)
//	{
//
//		printf("Assuming a plastic solar cell efficiency of 0.06\n");
//		//printf("Starting simulation for %i generations.\n", GA_maxgenerations);
//		WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
//		WriteMainOutputFile();
//		char mcommand[200];
//		sprintf(mcommand, "mkdir ");
//		strcpy(Clone_Dir, Unique_Dir);
//		strcat(mcommand, Clone_Dir);
//		strcat(mcommand, "models");
//		system(mcommand);
//
//		double roof_TotalPower[num_iterations+100]; // total power for each configuration (iteration)
//		//double roofArea = 0;	// calculated once only at the beginning of simulation. how much rooftop area
//
//		double allPositions[num_iterations+100][InstCopies_numberof*2]; // [iteration index][two coordinates (x,y) for each instance's roof position)
//
//		globalDoRoof=true; // to let "calculatepower..." know to set the || rays reference point to center of roof, not 5,5,5
//		globalRoofx = roof_xmax;
//		globalRoofy = roof_ymax;
//
//		// need to initialize files BEFORE the iterations, and write during them, NOT AT END WHEN DATA MAY BE ALL LOST!
//		// energy output file
//		strcpy(Clone_Dir, Unique_Dir);
//		char *output_concat1 = strcat(Clone_Dir, "RoofTestEnergies");
//		FILE *hOUTfile = fopen(output_concat1, "w");
//		fprintf(hOUTfile, "Configuration\t"); fprintf(hOUTfile, "Energy (kWh)\n");
//		fclose(hOUTfile);
//
//		// (x,y) corner position output file
//		strcpy(Clone_Dir, Unique_Dir);
//		char *output_concat2 = strcat(Clone_Dir, "RoofTestPositions");
//		FILE *eOUTfile = fopen(output_concat2, "w");
//		fclose(eOUTfile);
//
//
//		int pXwidth=5; // an XxY size grid with gridlengths of the unit structure width
//		int pYwidth=5;
//		double unit_Xwidth=10;
//		double unit_Ywidth=10;
//
//		bool position_taken[pXwidth][pYwidth];
//		for (int pxi=0; pxi<pXwidth; pxi++) {
//			for (int pyi=0; pyi<pYwidth; pyi++) {
//				position_taken[pxi][pyi] = false;
//			}
//		}
//
//		for (int ri=0; ri<num_iterations; ri++) {
//		// choose the coordinates
//			printf("Trial %i of %i\n", ri, num_iterations-1);
//
//			// new iteration, so reset these boolean checks
//			for (int pxi=0; pxi<pXwidth; pxi++) {
//				for (int pyi=0; pyi<pYwidth; pyi++) {
//					position_taken[pxi][pyi] = false;
//				}
//			}
//
//			for (int guy=0; guy<InstCopies_numberof; guy++) {
//
//				// its (0,0) corner may sit anywhere within roof-area such that no part of instance is outside roof boundary.
//				bool position_chosen=false; // false until the current unit gets a unique grid position
//				int pXind=0; int pYind=0;   // final unique position (indices for "position_taken") of this unit
//				while (position_chosen==false)
//				{
//					int guess_X = myRandom.boundedIntegerRandom(0,pXwidth);
//					int guess_Y = myRandom.boundedIntegerRandom(0,pYwidth);
//					if (position_taken[guess_X][guess_Y] == false) {
//						if ((guess_X!=pXwidth) and (guess_Y!=pYwidth))  // don't want it to be at limit
//						{
//							// this is so far a unique untaken position, so set it
//							pXind = guess_X;
//							pYind = guess_Y;
//							position_chosen=true;
//							position_taken[guess_X][guess_Y] = true; // this position is now taken
//						}
//					}
//					else {
//						// grid position is taken already, avoid overlap by redoing loop
//					}
//				}
//				// in exchange for commented-out code below
//				allPositions[ri][guy*2+0] = ((double)pXind)*unit_Xwidth;
//				allPositions[ri][guy*2+1] = ((double)pYind)*unit_Ywidth;
//
//				// 9/30/09 commented out these two lines which are all that is needed for RANDOM POSITIONS
//				//allPositions[ri][guy*2+0] = myRandom.random01() * (roof_xmax-xBound);
//				//allPositions[ri][guy*2+1] = myRandom.random01() * (roof_ymax-yBound);
//			}
//			AlterRoofStructure(wholeRoofModel, globalBasicCoordinates, allPositions[ri], InstCopies_numberof, copyCells_numberof);
//			Cell_3::cell_counter = 0;
//
//			numCells = copyCells_numberof * InstCopies_numberof;// must change the number of cells, since used by power-calculate
//			Cell_3 cellArray2[numCells];
//
//			modelPointer = wholeRoofModel;		// the actual model is the whole thing
//
//			BuildStructure(cellArray2, MODEL);
//			roof_TotalPower[ri]=newCalculatePowerDayCycle(cellArray2, wholeRoofModel);
//
//			// write datas
//			strcpy(Clone_Dir, Unique_Dir);
//			char *output_concat1a = strcat(Clone_Dir, "RoofTestEnergies");
//			FILE *haOUTfile = fopen(output_concat1a, "a");
//			double roofkWh = roof_TotalPower[ri]/(3600*1000);
//			fprintf(haOUTfile, "%i\t%f\n", ri, roofkWh);
//			fclose(haOUTfile);
//
//			strcpy(Clone_Dir, Unique_Dir);
//			char *output_concat2a = strcat(Clone_Dir, "RoofTestPositions");
//			FILE *eaOUTfile = fopen(output_concat2a, "a");
//			for (int dude=0; dude < InstCopies_numberof; dude++) {
//				fprintf(eaOUTfile, "%f\t%f\t", allPositions[ri][dude*2+0], allPositions[ri][dude*2+1]);
//			}
//			fprintf(eaOUTfile, "\n");
//			fclose(eaOUTfile);
//		}
//
//		int bestest_iteration=0;
//		double bestest_power=0;
//		for (int rii=0; rii<num_iterations; rii++)
//		{
//			if (roof_TotalPower[rii] > bestest_power)
//			{
//				bestest_power = roof_TotalPower[rii];
//				bestest_iteration = rii;
//			}
//		}
//		printf("best iteration was %i\n", bestest_iteration);
//
//		// now, also OUTPUT a .txt version of the BEST rooftop configuration model. (eg 12 cells per unit, 15 units--> 180 triangle)
//
//		Ggen_show=GA_gen_current; // OPENGL stuff follow this....
//	}
//	// display (also print to file) desired output
//	} // end SIMULATION option mainloop


	/*if (strcmp(argv[1],"rfmodel")==0)
	{
		// argv[2]= unit model path location/filenames
		// argv[3]= working directory for the location of the "RoofTestPositions" file
		// argv[4]= index integer of the BEST configuration for this job
		// argv[5]=number of units used


		strcpy(Input_CMD, argv[2]);
		getcwd(Main_Dir, 255);
		strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
		strcpy(Clone_Dir, Main_Dir);
		strcat(Clone_Dir, argv[3]);
		strcpy(Unique_Dir, Clone_Dir);
		strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
		strcpy(Clone_Dir, Main_Dir);

		int best_config = atoi(argv[4]);
		int num_units = atoi(argv[5]);

		strcpy(Clone_Dir, Unique_Dir);

		char *output_cat = strcat(Clone_Dir, OutModel_prefix);
		//strcat(output_cat, itoa(whichgen));
		char istring[200];
		sprintf(istring,"%i",whichgen);
		strcat(output_cat, istring);
		strcat(output_cat, "_cells.txt");
		//strcat(output_cat, itoa(numCells);

		FILE *OUTfile = fopen(output_cat, "w");

		for (int cff=0; cff<numCells; cff++)
		{
			for (int cdd=0; cdd<9; cdd++)
			{
				fprintf(OUTfile, "%f", globalParamsAvgBest[whichgen*numCells*9+cff*9+cdd]);
				if (cdd<8) {fprintf(OUTfile, ",");}
				if (cdd==8) {fprintf(OUTfile, "\n");}
			}
		}

		fclose(OUTfile);
	}*/
}





//int main(int argc, char *argv[])
//{
	//// depending on the first argument, 
	//// "sim" => regular simulation. followed by arguments "input_file" and "unique output directory"	
	//// or "view" => just view the outputs from an already completed simulation. 
		//// just follwed by argument "unique output directory" which is holding the stuff you want to look at

	//if (strcmp(argv[1],"test")==0)
	//{
		//// for teseting random things
		//double xxxx=1.89;
		//printf("test integer: %i\n", (int)xxxx);
	//}

	//if (strcmp(argv[1],"sim")==0)
	//{
	//printf("damnit\n");
	//// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
	//strcpy(Input_CMD, argv[2]);
	//getcwd(Main_Dir, 255);
	//strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//strcat(Clone_Dir, argv[3]);
	//strcpy(Unique_Dir, Clone_Dir);
	//strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
	///*cout << Main_Dir << endl;
	//cout << Unique_Dir << endl;
	//cout << Input_CMD << endl;*/
	////printf("%s\n%s\n%s\n", Main_Dir, Unique_Dir, Input_CMD);
	//strcpy(Clone_Dir, Main_Dir);
	//ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
////	cout << "got it" << endl;
	//if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
	//printf("Index of refraction: %f\n", index_ref_plastic);
	
	//// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
	//doSingleReflections = atoi(argv[4]);	

	///*for (int aaa=0; aaa<=9; aaa++) {
		//for (int bbb=0; bbb<=9; bbb++) {
			//for (int ccc=0; ccc<=9; ccc++) {
				//for (int ddd=0; ddd<=9; ddd++) {
					//for (int eee=0; eee<=9; eee++) {
						//for (int fff=0; fff<=9; fff++) {
					//printf("%i%i%i%i%i%i\n", aaa,bbb,ccc,ddd,eee,fff);}}}}}}*/
	
	//if (input_i==0)
	//{
		//// if not an example case
	//}
	//else
	//{
		//longitude   = INPUT_examples[input_i][0];   // longitude of location
		//if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
			//latitude    = INPUT_examples[input_i][1];   // latitude of location
	
			//if (globalInputLatitude==-401) {globalCellsOneSided=true;} 
		//}
		//else {
			//latitude = globalInputLatitude;
		//}
		//avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
		//avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
		//elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
		//timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
	//}
	//printf("latitude = %f\n",  latitude);
	////printf("Day %i, Month %i, year %i\n", startDay, startMonth, startYear);
	////printf("longitude %f latitude %f\n", longitude, latitude);

	///////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////

	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
	//FILE *fileBestF= fopen(output_catab, "w");
	//fclose(fileBestF);

	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//// 8/04/09: set all this qwikStore stuff to 0...to start out
	//for (int qi0=0; qi0<qwikIndMax; qi0++) {
		//qwikStoreEnergies[qi0]=0;
		//qwikStoreReflEnergies[qi0]=0;
		//for (int qi1=0; qi1<qwikCoordMax; qi1++) {
			//qwikStoreCoords[qi0][qi1]=0;
		//}
	//}

	//if (GA_willload_population==1)
	//{
		//printf("Loading initial population from model files\n");
		//numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
		//printf("number of coord: %i\n", numCoordinates);
		//double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
	////	double allModelDecisionVars[numCoordinates/*6/9*/GA_numDistinctModels]; // holds all decision vars for all models
		////modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
	////	allDecVarPointer=allModelDecisionVars;                   // global pointer
		///*FILE *inputPop = fopen(GA_initial_population_file,"w");
		//fclose (inputPop);*/
////		cout << Unique_Dir <<endl;
		//strcpy(Clone_Dir, Unique_Dir);
		//char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
		////printf("%s\n", myGAinput);
////		cout << myGAinput << endl;
		//FILE *inputPop = fopen(myGAinput, "w"); // write the population file
		//fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
		//printf("%i distinct models\n",GA_numDistinctModels);
		//for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
		//{
			//printf("first %i\n", nm);
			//ReadStructureVertices_Retrieve(modelCoordinates, nm);
			//printf("2nd %i\n", nm);
			//for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
			//{
				//for (int nc=0; nc < numCells; nc++)
				//{
					//// format: tab \t between all dec vars. \n after objective value written, before next individal
					//for (int np=0; np < 9; np++) {
						//fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);	
					//}
					////printf("%i\n",nc);
				//}
				//if (GA_numberof_objectives==2) 
					//{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
				//}
				//if (GA_numberof_objectives==1) {
					//fprintf(inputPop, "0\n");  // newline then next individual
				//}
			//}
			//// close just to have one file open at a time (ReadStructure... opens a file)
		//}
		//fclose(inputPop);
	//}
	//else
	//{
		//printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
	//}
	
	//printf("Assuming a plastic solar cell efficiency of 0.06\n");
	//printf("Starting simulation for %i generations.\n", GA_maxgenerations);
	//WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
	//WriteMainOutputFile();
	//char mcommand[200];
	//sprintf(mcommand, "mkdir ");
	//strcpy(Clone_Dir, Unique_Dir);
	//strcat(mcommand, Clone_Dir);
	//strcat(mcommand, "models");
	//system(mcommand);
	//mainGA();

	//Ggen_show=GA_gen_current; // OPENGL stuff follow this....
	//// display (also print to file) desired output
	//} // end SIMULATION option mainloop

	//if (strcmp(argv[1], "grid")==0) // expects two arguments: 1) output directory from previous run, 2) # grids (will be squared)
	//{
		//getcwd(Main_Dir, 255);
		//strcpy(Clone_Dir, Main_Dir);
		//strcpy(Main_Dir, ""); 		// make it blank just to avoid error later in model function
		//strcat(Clone_Dir, "/");
		//strcat(Clone_Dir, argv[2]);
		//printf("Searching in simulation output directory: %s\n", Clone_Dir);
		//strcpy(Unique_Dir, Clone_Dir);
		//strcat(Unique_Dir, "/");
		
		//// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
		//doSingleReflections = atoi(argv[4]);	
		
		//globalInputLatitude=-400; //added 8/13/09...temporary fix, assume im always "gridding" at berkeley for now
		//if (input_i==0)
		//{
			//// if not an example case
		//}
		//else
		//{
			//longitude   = INPUT_examples[input_i][0];   // longitude of location
			//if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
				//latitude    = INPUT_examples[input_i][1];   // latitude of location
		
				//if (globalInputLatitude==-401) {globalCellsOneSided=true;} 
			//}
			//else {
				//latitude = globalInputLatitude;
			//}
			//avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
			//avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
			//elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
			//timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
		//}
		//printf("latitude = %f\n",  latitude);

		//strcpy(Clone_Dir, Unique_Dir);
		//strcat(Clone_Dir, OutMain_name);	
		//ReadMainOutputFile( Clone_Dir);			// 1) get some typical variables, e.g., like # generations, if needed
		//raysPerCell=atoi(argv[3]); 			// whatever this is, the number is squared in next function
		//if (use_trigrid==1) {Get2DTriangleGrid();} 	// precalculate grid points
		
		/////////////// create new version of these  files//////////////
		///*strcpy(Clone_Dir, Unique_Dir);
		//char *output_catsc = strcat(Clone_Dir, OutPowerPerSubcell);
		//FILE *fileSubPOW= fopen(output_catsc, "w");
		//fclose(fileSubPOW);*/	// commented out 5/13/09 when added split by day iteration and "w" each one just once...no "a"

		//strcpy(Clone_Dir, Unique_Dir);
		//char *output_catgc = strcat(Clone_Dir, OutGridCoordinates);
		//FILE *fileGridCoord= fopen(output_catgc, "w");
		//fclose(fileGridCoord);
		///////////////////////////////////////////////////////////
		//strcpy(Clone_Dir, Unique_Dir);
		////char *GAmodel_prefix = "Model_gen";
		//char *GAmodel_suffix = "_cells.txt";
		//strcat(Clone_Dir, OutModel_prefix);
		//char  gnstring[200];
		//sprintf(gnstring,"%i", (GA_maxgenerations-1) ); 	// generation value was obtained by reading main output file
		//strcat(Clone_Dir, gnstring);
		//strcat(Clone_Dir, GAmodel_suffix);
		//strcpy(Input_modelALL[0], Clone_Dir);				// finally, full location of the model..
		//numCoordinates=ReadStructureVertices_Count(0);
		//double modelCoordinates[numCoordinates];
		//ReadStructureVertices_Retrieve(modelCoordinates, 0);

		//// I could just as well have obtained the coordinates from the AllBestGen file, but I think getting
		//// the coordinates from a model will be more versatile...or the same...
		//getSubCellData = true;
		//globalSetup = new GlobalSetup;			// just to make globalEvaluate happy (actually useless 'if' branch...)
		//globalHaveDonePreEval=1;			// fake this, since no GA is being started, just want one iteration
		//Write_tsampledata=true;		// to yes, write the coordinates of gridpoints first time
		//globalSetup->gaType = SGA;
		//globalSetup->finalNoOfObjectives=1;
		//globalSetup->noOfDecisionVariables=numCells*9;
		//double dummyObjective[5];
		//double dummyDouble[5];
		//int dummyInteger[5];
		//globalEvaluate(modelCoordinates, dummyObjective, dummyDouble, dummyDouble, dummyInteger); 
		//DisplayPowerOutputGridVersion(); // output a power vs. frac-hour data set with #gridpoints appended to filename
		
		//// get & write some angle analysis data. distribution weighted according to sum of energy over the day
		//GetZenithAngleDistribution(modelCoordinates);
		//GetAzimuthAngleDistribution(modelCoordinates);
		//Get2dAngleDistribution(modelCoordinates);	// for 2dimensional-bins (x,y) distribution
		//WriteAngleDistrFile();
	//}

	///////////////STRETCH WAIST STUFF!!!!!!!!!!!!!!1
	//if (strcmp(argv[1],"hf")==0)
	//{
	//// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
	//strcpy(Input_CMD, argv[2]);
	//getcwd(Main_Dir, 255);
	//strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//strcat(Clone_Dir, argv[3]);
	//strcpy(Unique_Dir, Clone_Dir);
	//strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
	//if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
	
	//// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
	//doSingleReflections = atoi(argv[4]);	
	
	//if (input_i==0)
	//{}
	//else
	//{
		//longitude   = INPUT_examples[input_i][0];   // longitude of location
		//if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
			//latitude    = INPUT_examples[input_i][1];   // latitude of location
	
			//if (globalInputLatitude==-401) {globalCellsOneSided=true;} 
		//}
		//else {
			//latitude = globalInputLatitude;
		//}
		//avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
		//avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
		//elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
		//timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
	//}
	//printf("latitude = %f\n",  latitude);

	///////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
	//FILE *fileBestF= fopen(output_catab, "w");
	//fclose(fileBestF);
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//// 8/04/09: set all this qwikStore stuff to 0...to start out
	//for (int qi0=0; qi0<qwikIndMax; qi0++) {
		//qwikStoreEnergies[qi0]=0;
		//qwikStoreReflEnergies[qi0]=0;
		//for (int qi1=0; qi1<qwikCoordMax; qi1++) {
			//qwikStoreCoords[qi0][qi1]=0;
		//}
	//}

	//double globalFunCoordinates[1000];
	//double globalBoxCoordinates[1000];
	
	//if (GA_willload_population==1)
	//{
		//printf("Loading initial population from model files\n");
		//numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
		//printf("number of coord: %i\n", numCoordinates);
		//double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
		//modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
		//strcpy(Clone_Dir, Unique_Dir);
		//char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
		//FILE *inputPop = fopen(myGAinput, "w"); // write the population file
		//fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
		//for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
		//{
			//ReadStructureVertices_Retrieve(modelCoordinates, nm);
			
			//for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
			//{
				//for (int nc=0; nc < numCells; nc++)
				//{
					//// format: tab \t between all dec vars. \n after objective value written, before next individal
					//for (int np=0; np < 9; np++) {
						//fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);

						//if (nm==0) {
							//globalFunCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
						//if (nm==1) {
							//globalBoxCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
	
					//}
				//}
				//if (GA_numberof_objectives==2) 
					//{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
				//}
				//if (GA_numberof_objectives==1) {
					//fprintf(inputPop, "0\n");  // newline then next individual
				//}
			//}
			//// close just to have one file open at a time (ReadStructure... opens a file)
		//}
		//fclose(inputPop);
	//}
	//else
	//{
		//printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
	//}
	
	//printf("Assuming a plastic solar cell efficiency of 0.06\n");
	////printf("Starting simulation for %i generations.\n", GA_maxgenerations);
	//WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
	//WriteMainOutputFile();
	//char mcommand[200];
	//sprintf(mcommand, "mkdir ");
	//strcpy(Clone_Dir, Unique_Dir);
	//strcat(mcommand, Clone_Dir);
	//strcat(mcommand, "models");
	//system(mcommand);
	
	//double heighttest_BoxPower[201]; // box power at each height
	//double heighttest_FunPower[201]; // funnel power at each height
	//for (int jij=0; jij<201; jij++) {
		//heighttest_BoxPower[jij]=0;
		//heighttest_FunPower[jij]=0;
	//}
	
	//int maxHeightGoto=100;
	//for (int ghi=2; ghi <= maxHeightGoto; ghi+=2) {
		//printf("\nheight=%i\n",ghi);
		//AlterBoxCoordinates(globalBoxCoordinates, ghi);
		//Cell_3::cell_counter = 0;            
		//Cell_3 cellArray[numCells];
		//modelPointer = globalBoxCoordinates;
		//BuildStructure(cellArray, MODEL);   
		//heighttest_BoxPower[ghi]=newCalculatePowerDayCycle(cellArray, globalBoxCoordinates);

		//AlterFunnelCoordinates(globalFunCoordinates, ghi); // takes coordinates and height as arguments
		//printf("%f\n",globalFunCoordinates[2]);
		//Cell_3::cell_counter = 0;            
		//Cell_3 cellArray2[numCells];
		//modelPointer = globalFunCoordinates;
		//BuildStructure(cellArray2, MODEL);   
		//heighttest_FunPower[ghi]=newCalculatePowerDayCycle(cellArray2, globalFunCoordinates);
	//}

	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_concat1 = strcat(Clone_Dir, "HeightBoxFunnelCompare");
	//FILE *hOUTfile = fopen(output_concat1, "w");
	
	//fprintf(hOUTfile, "Height\t"); fprintf(hOUTfile, "Box (kWh)\t"); fprintf(hOUTfile, "Funnel (kWh)\t"); 
	//fprintf(hOUTfile, "percent diff\n"); 
	//for (int cff=2; cff<= maxHeightGoto; cff+=2)	
	//{	
		//double boxkWh = heighttest_BoxPower[cff]/(3600*1000);
		//double funkWh = heighttest_FunPower[cff]/(3600*1000);
		//double perC = ( 100*fabs(heighttest_FunPower[cff]-heighttest_BoxPower[cff])/heighttest_BoxPower[cff] );
		//fprintf(hOUTfile, "%i\t%f\t%f\t%f\n", cff, boxkWh, funkWh, perC); 
	//}
	//fclose(hOUTfile);

	//Ggen_show=GA_gen_current; // OPENGL stuff follow this....
	//// display (also print to file) desired output
	//} // end SIMULATION option mainloop

	
	//////////?WAIST FUNNEL STUFF////////////////
	//if (strcmp(argv[1],"wf")==0)
	//{
	//// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
	//strcpy(Input_CMD, argv[2]);
	//getcwd(Main_Dir, 255);
	//strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//strcat(Clone_Dir, argv[3]);
	//strcpy(Unique_Dir, Clone_Dir);
	//strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
	//if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
	
	//// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
	//doSingleReflections = atoi(argv[4]);	
	
	//if (input_i==0)
	//{}
	//else
	//{
		//longitude   = INPUT_examples[input_i][0];   // longitude of location
		//if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
			//latitude    = INPUT_examples[input_i][1];   // latitude of location
	
			//if (globalInputLatitude==-401) {globalCellsOneSided=true;} 
		//}
		//else {
			//latitude = globalInputLatitude;
		//}
		//avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
		//avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
		//elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
		//timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
	//}
	//printf("latitude = %f\n",  latitude);

	///////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
	//FILE *fileBestF= fopen(output_catab, "w");
	//fclose(fileBestF);
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//// 8/04/09: set all this qwikStore stuff to 0...to start out
	//for (int qi0=0; qi0<qwikIndMax; qi0++) {
		//qwikStoreEnergies[qi0]=0;
		//qwikStoreReflEnergies[qi0]=0;
		//for (int qi1=0; qi1<qwikCoordMax; qi1++) {
			//qwikStoreCoords[qi0][qi1]=0;
		//}
	//}

	//double globalFunCoordinates[1000];
	////double globalBoxCoordinates[1000];
	
	//if (GA_willload_population==1)
	//{
		//printf("Loading initial population from model files\n");
		//numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
		//printf("number of coord: %i\n", numCoordinates);
		//double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
		//modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
		//strcpy(Clone_Dir, Unique_Dir);
		//char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
		//FILE *inputPop = fopen(myGAinput, "w"); // write the population file
		//fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
		//for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
		//{
			//ReadStructureVertices_Retrieve(modelCoordinates, nm);
			
			//for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
			//{
				//for (int nc=0; nc < numCells; nc++)
				//{
					//// format: tab \t between all dec vars. \n after objective value written, before next individal
					//for (int np=0; np < 9; np++) {
						//fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);

						//if (nm==0) {
							//globalFunCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
	
					//}
				//}
				//if (GA_numberof_objectives==2) 
					//{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
				//}
				//if (GA_numberof_objectives==1) {
					//fprintf(inputPop, "0\n");  // newline then next individual
				//}
			//}
			//// close just to have one file open at a time (ReadStructure... opens a file)
		//}
		//fclose(inputPop);
	//}
	//else
	//{
		//printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
	//}
	
	//printf("Assuming a plastic solar cell efficiency of 0.06\n");
	////printf("Starting simulation for %i generations.\n", GA_maxgenerations);
	//WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
	//WriteMainOutputFile();
	//char mcommand[200];
	//sprintf(mcommand, "mkdir ");
	//strcpy(Clone_Dir, Unique_Dir);
	//strcat(mcommand, Clone_Dir);
	//strcat(mcommand, "models");
	//system(mcommand);
	
	//double heighttest_FunPower[50][50]; // funnel power at each height
	
	//int maxYGoto=9;
	//int maxZGoto=21;
	//for (int fiy=1; fiy <= maxYGoto; fiy+=1) {
		//for (int fiz=1; fiz <= maxZGoto; fiz+=1) {
			//printf("\nwaist_y=%i, waist_z=%i\n",fiy,fiz);

			//AlterFunnelWaist(globalFunCoordinates, fiy, fiz); // takes coordinates and height as arguments
			////printf("%f\n",globalFunCoordinates[2]);
			//Cell_3::cell_counter = 0;            
			//Cell_3 cellArray2[numCells];
			//modelPointer = globalFunCoordinates;
			//BuildStructure(cellArray2, MODEL);   
			//heighttest_FunPower[fiy][fiz]=newCalculatePowerDayCycle(cellArray2, globalFunCoordinates);
		//}
	//}

	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_concat1 = strcat(Clone_Dir, "WaistFunnelCompare");
	//FILE *hOUTfile = fopen(output_concat1, "w");
	
	//fprintf(hOUTfile, "y\t"); fprintf(hOUTfile, "z\t"); fprintf(hOUTfile, "Energy (kWh)\n");  
	//for (int fiy=1; fiy <= maxYGoto; fiy+=1) {
		//for (int fiz=1; fiz <= maxZGoto; fiz+=1) {	
			//double funkWh = heighttest_FunPower[fiy][fiz]/(3600*1000);
			//fprintf(hOUTfile, "%i\t%i\t%f\n", fiy, fiz, funkWh); 
		//}
	//}
	//fclose(hOUTfile);

	//Ggen_show=GA_gen_current; // OPENGL stuff follow this....
	//// display (also print to file) desired output
	//} // end SIMULATION option mainloop
	
	////////////////////////////////////////////////////////////////////
	//////////?ROOF ARRANGEMENT SUMULATION OPTIMIZATION ////////////////
	////////////////////////////////////////////////////////////////////
	//if (strcmp(argv[1],"rf")==0)
	//{
	//// if im doing above, i need to increase index of the below stuffs by 1, and put this in conditional block for "sim"
	//strcpy(Input_CMD, argv[2]);
	//getcwd(Main_Dir, 255);
	//strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//strcat(Clone_Dir, argv[3]);
	//strcpy(Unique_Dir, Clone_Dir);
	//strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
	//strcpy(Clone_Dir, Main_Dir);
	//ReadMainInputFile(strcat(Clone_Dir, Input_CMD));
	//if (use_trigrid==1) {Get2DTriangleGrid();} // precalculate grid points
	
	//// newly-added (6/4/09) command line parameter: doReflections = true or false (1 or 0)
	//// argument 4 = roofx dimension
	//// argument 5 = roofy dimension
	//// argument 6 = # of instances of unit model to scatter around
	//// argument 7 = reflections turned on
	////double roof_xmin = 0;		// bounds of the roof on x,y plane
	////double roof_ymin = 0;		
	//double roof_xmax = atof(argv[4]);
	//double roof_ymax = atof(argv[5]);
	//int InstCopies_numberof = atoi(argv[6]);
	
	//doSingleReflections = atoi(argv[7]);	
	
	//if (input_i==0)
	//{}
	//else
	//{
		//longitude   = INPUT_examples[input_i][0];   // longitude of location
		//if ( (globalInputLatitude==-400) || (globalInputLatitude==-401) ) {
			//latitude    = INPUT_examples[input_i][1];   // latitude of location
	
			//if (globalInputLatitude==-401) {globalCellsOneSided=true;} 
		//}
		//else {
			//latitude = globalInputLatitude;
		//}
		//avgTemp     = INPUT_examples[input_i][2];   // local annual temp in location (C)
		//avgPressure = INPUT_examples[input_i][3];   // local annual pressure in location (millibars)
		//elevation   = INPUT_examples[input_i][4];   // elevation of location (meters)
		//timeZone    = INPUT_examples[input_i][5];   // timeZone (hours west of GMT are negative)
	//}
	//printf("latitude = %f\n",  latitude);

	///////DELETE OLD DATA BEFORE APPENDING TO THIS FILENAME//////////////////////////////////////////////
	//strcpy(Clone_Dir, Unique_Dir);
	//char *output_catab = strcat(Clone_Dir, OutAppBest_Obj);
	//FILE *fileBestF= fopen(output_catab, "w");
	//fclose(fileBestF);
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//// 8/04/09: set all this qwikStore stuff to 0...to start out
	//for (int qi0=0; qi0<qwikIndMax; qi0++) {
		//qwikStoreEnergies[qi0]=0;
		//qwikStoreReflEnergies[qi0]=0;
		//for (int qi1=0; qi1<qwikCoordMax; qi1++) {
			//qwikStoreCoords[qi0][qi1]=0;
		//}
	//}

	//double globalBasicCoordinates[1000]; // basic coordinates of unit structure (at x,y = 0,0 corner)
	////double globalBoxCoordinates[1000];
	
	//if (GA_willload_population==1)
	//{
		//printf("Loading initial population from model files\n");
		//numCoordinates=ReadStructureVertices_Count(0);            // returns # of coordinates stored in model, (must be multiple of 9)
		//printf("number of coord: %i\n", numCoordinates);
		//double modelCoordinates[numCoordinates];                 // holds just coordinates of ONE model individual
		//modelPointer=modelCoordinates;                           // global pointer to first element in modelCoordinates[..]
		//strcpy(Clone_Dir, Unique_Dir);
		//char *myGAinput = strcat(Clone_Dir, GA_initial_population_file);
		//FILE *inputPop = fopen(myGAinput, "w"); // write the population file
		//fprintf(inputPop, "%i\n", GA_population); // first line is number of individuals
		//for (int nm=0; nm < GA_numDistinctModels; nm++)          // for each distinct model, load its coordinates
		//{
			
			//ReadStructureVertices_Retrieve(modelCoordinates, nm);
			
			//for (int nclone=0; nclone < GA_population/GA_numDistinctModels; nclone++)
			//{
				//for (int nc=0; nc < numCells; nc++)
				//{
					//// format: tab \t between all dec vars. \n after objective value written, before next individal
					//for (int np=0; np < 9; np++) {
						//fprintf(inputPop, "%f\t", modelCoordinates[nc*9 + np]);

						//if (nm==0) {
							//globalBasicCoordinates[nc*9+np] = modelCoordinates[nc*9+np];}
	
					//}
				//}
				//if (GA_numberof_objectives==2) 
					//{fprintf(inputPop, "0\t"); fprintf(inputPop, "0\n"); // already distiguishes NSA or SGA
				//}
				//if (GA_numberof_objectives==1) {
					//fprintf(inputPop, "0\n");  // newline then next individual
				//}
			//}
			//// close just to have one file open at a time (ReadStructure... opens a file)
		//}
		//fclose(inputPop);
	//}
	//else
	//{
		//printf("Generating random starting population within coordinate bounds: models in input file not used...\n");
	//}
	//printf("check1\n");
	
	//int simORmodel = 1; // 1 for just print out best model

	//double wholeRoofModel[InstCopies_numberof*numCells*9];
	//int copyCells_numberof = numCells;
	//int num_iterations = GA_maxgenerations; // may as well use same number

	//if (simORmodel==1)
	//{
		//double got_positions/*[num_iterations+100]*/[InstCopies_numberof*2];
		//int Best_config = atoi(argv[7]);// just for simplicity, use the singleDoreflections input for this as any integer
		//printf("check2\n");
		//ReadRoofPositions(got_positions, InstCopies_numberof*2, Best_config);
		//printf("check3\n");
		
		//AlterRoofStructure(wholeRoofModel, globalBasicCoordinates, got_positions, InstCopies_numberof, numCells);

		//strcpy(Clone_Dir, Unique_Dir);

		//char *output_catM = strcat(Clone_Dir, "BestFullModel.txt");
		////strcat(output_cat, itoa(whichgen));
		///*char istring[200];
		//sprintf(istring,"%i",whichgen);	
		//strcat(output_cat, istring);
		//strcat(output_cat, "_cells.txt");*/
		////strcat(output_cat, itoa(numCells);
		//FILE *OUTfileM = fopen(output_catM, "w");
	
		//for (int cff=0; cff<InstCopies_numberof; cff++)
		//{	
			//for (int ccc=0; ccc<copyCells_numberof; ccc++)
			//{
				//for (int cdd=0; cdd<9; cdd++)
				//{
					////fprintf(OUTfile, "%f", globalParamsAvgBest[whichgen*numCells*9+cff*9+cdd]);
					//fprintf(OUTfileM, "%f", wholeRoofModel[cff*copyCells_numberof*9 + ccc*9 + cdd]);
					//if (cdd<8) {fprintf(OUTfileM, ",");}
					//if (cdd==8) {fprintf(OUTfileM, "\n");}
				//}
			//}
		//}

		//fclose(OUTfileM);
	//}

	//if (simORmodel==0)
	//{

		//printf("Assuming a plastic solar cell efficiency of 0.06\n");
		////printf("Starting simulation for %i generations.\n", GA_maxgenerations);
		//WriteGAInputFile(numCells, 9);  //  times 6 for position/rotation picture. x9 for coordinate picture
		//WriteMainOutputFile();
		//char mcommand[200];
		//sprintf(mcommand, "mkdir ");
		//strcpy(Clone_Dir, Unique_Dir);
		//strcat(mcommand, Clone_Dir);
		//strcat(mcommand, "models");
		//system(mcommand);
	
		//double roof_TotalPower[num_iterations+100]; // total power for each configuration (iteration)
		////double roofArea = 0;	// calculated once only at the beginning of simulation. how much rooftop area 

		//double allPositions[num_iterations+100][InstCopies_numberof*2]; // [iteration index][two coordinates (x,y) for each instance's roof position)

		//globalDoRoof=true; // to let "calculatepower..." know to set the || rays reference point to center of roof, not 5,5,5
		//globalRoofx = roof_xmax;
		//globalRoofy = roof_ymax;

		//// need to initialize files BEFORE the iterations, and write during them, NOT AT END WHEN DATA MAY BE ALL LOST!
		//// energy output file
		//strcpy(Clone_Dir, Unique_Dir);
		//char *output_concat1 = strcat(Clone_Dir, "RoofTestEnergies");
		//FILE *hOUTfile = fopen(output_concat1, "w");
		//fprintf(hOUTfile, "Configuration\t"); fprintf(hOUTfile, "Energy (kWh)\n");
		//fclose(hOUTfile);

		//// (x,y) corner position output file
		//strcpy(Clone_Dir, Unique_Dir);
		//char *output_concat2 = strcat(Clone_Dir, "RoofTestPositions");
		//FILE *eOUTfile = fopen(output_concat2, "w");
		//fclose(eOUTfile);


		//int pXwidth=5; // an XxY size grid with gridlengths of the unit structure width
		//int pYwidth=5;
		//double unit_Xwidth=10;
		//double unit_Ywidth=10;

		//bool position_taken[pXwidth][pYwidth];
		//for (int pxi=0; pxi<pXwidth; pxi++) {
			//for (int pyi=0; pyi<pYwidth; pyi++) {
				//position_taken[pxi][pyi] = false;			
			//}
		//}

		//for (int ri=0; ri<num_iterations; ri++) {
		//// choose the coordinates
			//printf("Trial %i of %i\n", ri, num_iterations-1);
		
			//// new iteration, so reset these boolean checks
			//for (int pxi=0; pxi<pXwidth; pxi++) {
				//for (int pyi=0; pyi<pYwidth; pyi++) {
					//position_taken[pxi][pyi] = false;			
				//}
			//}

			//for (int guy=0; guy<InstCopies_numberof; guy++) {
			
				//// its (0,0) corner may sit anywhere within roof-area such that no part of instance is outside roof boundary.			
				//bool position_chosen=false; // false until the current unit gets a unique grid position
				//int pXind=0; int pYind=0;   // final unique position (indices for "position_taken") of this unit
				//while (position_chosen==false) 
				//{
					//int guess_X = myRandom.boundedIntegerRandom(0,pXwidth);
					//int guess_Y = myRandom.boundedIntegerRandom(0,pYwidth); 
					//if (position_taken[guess_X][guess_Y] == false) {
						//if ((guess_X!=pXwidth) and (guess_Y!=pYwidth))  // don't want it to be at limit
						//{
							//// this is so far a unique untaken position, so set it
							//pXind = guess_X;
							//pYind = guess_Y;
							//position_chosen=true;
							//position_taken[guess_X][guess_Y] = true; // this position is now taken
						//}
					//}
					//else {
						//// grid position is taken already, avoid overlap by redoing loop
					//}
				//}
				//// in exchange for commented-out code below				
				//allPositions[ri][guy*2+0] = ((double)pXind)*unit_Xwidth;
				//allPositions[ri][guy*2+1] = ((double)pYind)*unit_Ywidth;

				//// 9/30/09 commented out these two lines which are all that is needed for RANDOM POSITIONS
				////allPositions[ri][guy*2+0] = myRandom.random01() * (roof_xmax-xBound);
				////allPositions[ri][guy*2+1] = myRandom.random01() * (roof_ymax-yBound);
			//}
			//AlterRoofStructure(wholeRoofModel, globalBasicCoordinates, allPositions[ri], InstCopies_numberof, copyCells_numberof);
			//Cell_3::cell_counter = 0;            
		
			//numCells = copyCells_numberof * InstCopies_numberof;// must change the number of cells, since used by power-calculate
			//Cell_3 cellArray2[numCells];
		
			//modelPointer = wholeRoofModel;		// the actual model is the whole thing
		
			//BuildStructure(cellArray2, MODEL);   
			//roof_TotalPower[ri]=newCalculatePowerDayCycle(cellArray2, wholeRoofModel);

			//// write datas
			//strcpy(Clone_Dir, Unique_Dir);
			//char *output_concat1a = strcat(Clone_Dir, "RoofTestEnergies");
			//FILE *haOUTfile = fopen(output_concat1a, "a");
			//double roofkWh = roof_TotalPower[ri]/(3600*1000);
			//fprintf(haOUTfile, "%i\t%f\n", ri, roofkWh);
			//fclose(haOUTfile);

			//strcpy(Clone_Dir, Unique_Dir);
			//char *output_concat2a = strcat(Clone_Dir, "RoofTestPositions");
			//FILE *eaOUTfile = fopen(output_concat2a, "a");
			//for (int dude=0; dude < InstCopies_numberof; dude++) {
				//fprintf(eaOUTfile, "%f\t%f\t", allPositions[ri][dude*2+0], allPositions[ri][dude*2+1]);
			//}
			//fprintf(eaOUTfile, "\n");
			//fclose(eaOUTfile);
		//}	

		//int bestest_iteration=0;
		//double bestest_power=0;
		//for (int rii=0; rii<num_iterations; rii++)
		//{
			//if (roof_TotalPower[rii] > bestest_power)
			//{
				//bestest_power = roof_TotalPower[rii];
				//bestest_iteration = rii;
			//}
		//}
		//printf("best iteration was %i\n", bestest_iteration);

		//// now, also OUTPUT a .txt version of the BEST rooftop configuration model. (eg 12 cells per unit, 15 units--> 180 triangle)

		//Ggen_show=GA_gen_current; // OPENGL stuff follow this....
	//}
	//// display (also print to file) desired output
	//} // end SIMULATION option mainloop

	
	///*if (strcmp(argv[1],"rfmodel")==0) 
	//{
		//// argv[2]= unit model path location/filenames
		//// argv[3]= working directory for the location of the "RoofTestPositions" file
		//// argv[4]= index integer of the BEST configuration for this job
		//// argv[5]=number of units used


		//strcpy(Input_CMD, argv[2]);
		//getcwd(Main_Dir, 255);
		//strcat(Main_Dir,"/");               // final form of Main_Dir which will be copied many times
		//strcpy(Clone_Dir, Main_Dir);
		//strcat(Clone_Dir, argv[3]);
		//strcpy(Unique_Dir, Clone_Dir);
		//strcat(Unique_Dir,"/");             // final form of Unique_Dir which will be copied many times
		//strcpy(Clone_Dir, Main_Dir);

		//int best_config = atoi(argv[4]);
		//int num_units = atoi(argv[5]);

		//strcpy(Clone_Dir, Unique_Dir);

		//char *output_cat = strcat(Clone_Dir, OutModel_prefix);
		////strcat(output_cat, itoa(whichgen));
		//char istring[200];
		//sprintf(istring,"%i",whichgen);	
		//strcat(output_cat, istring);
		//strcat(output_cat, "_cells.txt");
		////strcat(output_cat, itoa(numCells);

		//FILE *OUTfile = fopen(output_cat, "w");
	
		//for (int cff=0; cff<numCells; cff++)
		//{	
			//for (int cdd=0; cdd<9; cdd++)
			//{
				//fprintf(OUTfile, "%f", globalParamsAvgBest[whichgen*numCells*9+cff*9+cdd]);
				//if (cdd<8) {fprintf(OUTfile, ",");}
				//if (cdd==8) {fprintf(OUTfile, "\n");}
			//}
		//}

		//fclose(OUTfile);
	//}*/
//}

