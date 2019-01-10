/*UCNTransport by A. T. Holley*/
//Version: 29.9

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>

//Running Parameters
#define N 1000  //The number of particles to propagate through the geometry.
#define BATCH "OFF"  //Flag to turn batch mode on and off.
#define CHECK "OFF" //Flag to control whether detailed information on intermediate solutions etc. is displayed.
#define CONSOLE "ON"  //Flag to control whether WARNINGS and ERRORS are displayed in the console ("ON") or only ERRORS ("OFF").
#define EVENTS "ON"  //Flag to control whether all events (ON) or just some events (currently only detector and trajectory events) are written to file.
#define TRAJFLAG "OFF"  //Flag to indicate whether or not to record detailed trajectory information.
#define TRAJTS .001  //How often to record a trajectory point (in seconds).
#define REGFILE "Regionfile"  //The name of the 'regions' file to be read.
#define CONFILE "Connexfile"  //The name of the 'connex' file to be read.
#define RSEED "RAND"  //Random seed for the random number generator, which should be an unsigned integer or "RAND" if a new seed is desired each simulation
#define MZERO 1e-15  //The zero boundary value used in 'mathzero' to control roundoff error.
#define TZERO 1e-9  //The zero boundary value used in 'timezero' to control roundoff error.
#define GZERO 1e-9  //The zero boundary value used in 'move' to check for correct intersections-- essentially the fuzziness of the geometry.
#define VCUTOFF 7.8  //Cut-off speed for a v^2 dv speed dsitribution.
#define MONOENERGY 5.0 //Speed for a monoenergetic energy distribution
#define BEAMTIME 0 //Number of seconds that neutrons are being produced.
#define SHUTTERTIME 0 //Time when the shutter opens.
#define DETCUTOFF 0 //maximum speed at which a neutron can be detected for specified detector below. No cutoff if set to zero
#define CUTOFFREG 0 //region for which speed cutoff is implemented

//Constants
#define nMASS 1.67493e-27 //neutron mass in [kg]
#define PI 3.1415927
#define G 9.8 //acceleration due to gravity in [m/s^2]
#define nGAMMA 1.8325e8 //neutron gyromagnetic ratio in [1/(s*T)]
#define hBAR 1.05457e-34 //Plank's constant over 2pi in [J*s]
#define neV2J 1.602177e-28 //Conversion factor from neV to Joules


time_t time(time_t *t); //declare time() so that it can be used to generate a new random seed each simulation if that option is specified

struct complex {  //Structure to hold a complex number for use in 'solve3' and 'ccuberoot'.
    double r;
    double i;
};

double regions[100][7]; //Holds the geometry information for each region. The component [region#][#] is assigned as follows:
                      //0 -> surface type ; 1 -> Dim x or Diameter ; 2 -> Dim y or Length ; 3 -> Dim z
                     //4 -> theta orientation of the region in RADIANS ; 5 -> phi orientation in RADIANS ; 6 -> psi orientation in RADIANS.
double basepoints[100][3]; //Holds each region's basepoint. The component [region#][#] is assigned as follows:
                         //0 -> bpx ; 1 -> bpy ; 2 -> bpz.
double cplanes[100][8]; //Holds the orientation information for each cut-plane. A cut-plane with theta=0 and phi=0 is defined to be perpendicular to the z-axis.
                      //The component [region# associated with cut-plane][#] is assigned as follows:
                     //0 -> theta orientation of cut-plane [RADIANS] ; 1 -> phi orientation [RADIANS] ; 2 -> special handling code ; 3 -> Special-handling code 2,7
                    //value ; 4 -> internal special handling code ; 5 -> T Cut-plane offset ; 6 -> Region-shift Flag ; 7 -> psi orientation [RADIANS]
                   //T-junc. intersection codes: 0 -> Region has not been changed
                  //                            -1 -> Region has been shrunk.
                 //                              1 -> Region has been expanded.
                //Internal Special Handling Codes: 0 => No internal special handling
               //                                  1 => Cylinder region - Box region connection
              //                                   2 => The region is T'ed into (i.e. a region connects into its side).
             //                                    3 => The region T's into another region, BEGINNING in the side of a contiguous section.
            //                                     4 => The region T's into another region, ENDING in the side of a contiguous section.
double rparams[100][16]; //Holds the physics parameters for each region. The component [region#][#] is assigned as follows:
                       //0 -> dB/dx ; 1 -> dB/dy ; 2 -> dB/dz ; 3 -> Spec ; 4 -> Losspb ; 5 -> Depol. ; 6 -> Wall Potential ; 7 -> Bulk Potential ; 8 -> Surf/Scatt
                      //9 -> Absorption ; 10 -> Det ; 11 -> Prop Mode ; 12 -> Scatt Model ; 13 -> Loss Model ; 14 -> Depol Model ; 15 -> open slot for future
int connex[100][6]; //Holds the region connection information. For each region (first slot), the array holds up to six regions that are connected to it.
                   //A '-1' is placed in the first empty slot.

int simexec(char [],char [],int *); //Function to execute the main simulation logic.
                                                                          //Argument #1: Passes name of 'events' file.
                                                                         // Argument #2: Passes name of 'detectors' file.
                                                                        //  Argument #3: Allows the function to report back the integrated counts for five detectors.
                                                                       //                Note that what will be passed is a pointer to an array.
                                                                      //The function returns -1 on error and 0 otherwise.
void geomread(char [],char []); //Function to fill the above arrays from 'regions' and 'connex'.
                               //The arguments pass the names of the files to read.
void geomprint(char []); //Function to write the above array values into a text file for checking purposes.
                                                              //Arugment #1: Name of the file to write to 
void rota(double,double,double,double,double,double,double *,double *,double *); //Function to rotate a vector in the following order:
                                                                                // 1. Rotate by psi around global y-axis.
                                                                               //  2. Rotate by theta around x-axis that would result by rotating global x-axis via 1
                                                                              //   3. Rotate by phi around z-axis that would result by rotating global z-axis via 1
                                                                             //Arguments: (psi,theta,phi,xi,yi,zi,xf,yf,zf)
                                                                            //Note that the last three are pointers so that values in calling function can be set.
                                                                           //NB: This carries out an ACTIVE rotation. This is also the INVERSE of the associated
                                                                          //passive rotation. [Angles in RADIANS.]
void rotp(double,double,double,double,double,double,double *,double *,double *); //Function to rotate the global coordinate system in the follwing order:
                                                                                // 1. Rotate by psi around global y-axis.
                                                                               //  2. Rotate by theta around the x-axis resulting from 1
                                                                              //   3. Rotate by phi around z-axis which resulted after step 1 (and before step 2)
                                                                             //Arguments: (psi,theta,phi,xi,yi,zi,xf,yf,zf)
                                                                            //Note that the last three are pointers so that values in calling function can be set.
                                                                           //NB: This carries out a PASSIVE rotation. This is also the INVERSE of the associated
                                                                          //active rotation. [Angles in RADIANS.]
double solve2(double,double,double); //Function to find the smallest non-zero real root of: at^2 + bt + c = 0.
                                    //Arguments are: (a,b,c) | Returns: Smallest non-zero real root.
double solve4(double,double,double,double,double); //Function to find the smallest non-zero real root of: at^4 + bt^3 + ct^2 + dt + e = 0.
                                                  //Arguments are: (a,b,c,d,e) | Returns: Smallest non-zero real root or zero.
int newton(double,double,double,double,double,double*,double*); //Function to apply Newton's Method to a root of a polynomial of the form:
                                                               //                     t^4 + a3 t^3 + a2 t^2 + a1 t + a0 = 0.
                                                              //The first four arguments are, respectively, the coeficients a3, a2, a2, a0. The fifth argument 
                                                             //passes the starting value for the root. The first pointer argument allows the function to return the
                                                            //value of the root after at most 100 iterations of Newton's method. The second pointer argument allows
                                                           //the function to to return the polynomial at the improved value of the root (which of course should be
                                                          //zero). The funtion returns: 0 => No errors | -1 => error.
void solve3(double,double,double,double,double*,double*,double*); //Function to find the real roots of a cubic polynomial of the form ax^3 + bx^2 + cx + d = 0.
                                                                 //Arguments are: (a,b,c,d,r1,r2,r3)
struct complex ccuberoot(struct complex); //Function to calculate the complex cube root of a complex number.
                                         //Argument: Input Value | Returns: Complex number which is the cube root of input value.
double grn(void); //Function which returns a random number in [0,1].
double timezero(double); //Function which checks to see if a double-precision number is consistent with zero (i.e. within a certain range of zero)
                        //and if so sets the value exactly equal to zero. This function is used when a zero result is expected but roundoff error might
                       //have left a residual small value. (For example, any intersection time that is too small could be the result of roundoff error
                      //allowing the solvers to find the current intersection as a very small non-zero time solution.) Note that the zero boundary set
                     //in 'timezero' corresponds to distances of less than an angstrom for particle's with UCN-type speeds. The argument is the number to check,
                    //and the function returns the number if it is not consistent with zero and returns zero if the number is consistent with zero.
double mathzero(double); //This function is exactly like 'timezero' except that the zero boundary is defined differently. This function is used for testing
                        //math parameters for zero consistency (e.g. for sign determination in 'solve4') rather than to eliminate solution with times that are
                       //too small.
void poof(int,double,int,int); //Function to create a particle. The first two arguments are the following pieces of information:
                       //         Argument #1 => Region number in which to start the particle. NB: No error trapping for an undefined region.
                      //          Argument #2 => Particle will be created in a plane parallel to the region's cut-plane a distance into the region from the
                     //                          cut-plane given by this argument. Distance is in [m]. NB: No error trapping for offstes longer than the region.
                    //The third argument allows one to select from different speed distributions:
                   //          0 => A v^2 dv distribution of speeds with cutoff defined in header.
                  //           1 => A distribution (defined in the function) to be calculated via the Monte Carlo method.
                //             2 => A monoenergetic distribution with speed defined in header
                 //The fourth argument allows one to select from different angular distributions:
                //             0 => A cos weighted isotropic distribution (appropriate for isotropic illumination of a plane) directed into 2pi into the region.
               //              1 => A truly isotropic distribution directed into 4pi: sin(theta) d(theta) d(phi), i.e. no cos weighting.
              //Note that for 4pi  starting distributions one must always specify some non-zero initial offset.
int propagate(void); //Function to find the smallest-time intersection of the particle's trajectory with the geometry, update the particle's dynamical variables
                    //(but NOT region number) and write trajectory events if tracking is turned ON. The function returns an event code (see below).
int bounce(int,double,int,double,double,double); //Function to perform a bounce based on the models and values specified for the region number given by the first agument.
                           //The function returns -1 on error, 1 for wall penetration, 2 for physical wall loss, and 0 otherwise.
                          //The third argument which is passed is a flag indicating whether calculated normals should be reversed. (This needs to be done when the
                         //intersection with a cut-plane is from outside the cut-plane's region.)
                        //             1 => Calculated normal is an INWARD NORMAL | -1 => Calculated normal is an OUTWARD NORMAL.
                       //The second arument passed is the value of the material potential used to check for wall penetration.
                      //last three arguments, normally zero, are used if you want to pass a predetermined normal
int nsbounce(int,int,int,double,double,double); //Function to calculate the results of a non-specular bounce based on the selected scattering model.
                                               //Argument #1: The region number from which to obtain values for spec and scattering model.
                                              // Argument #2: The direction flag passed to 'bounce'.
                                             //  Argument #2: Flag indicating whether the particle depolarized during the bounce.
                                            //   Argument #4,5,6: vxb,vyb,vzb i.e. the particle's velocity in the bounce system.
                                           //The function sets the post-bounce particle velocity and return 0 or returns -1 for an error.
int bouncetype(int,double); //Function to determine if a bounce is specular or non-specular. It returns 0 for a specular bounce, 1 for a non-specular bounce, and
                           //returns -1 for an error. The first argument passes the region number specifying the model/values to use and the second argument passes
                          //the angle that the velocity vector makes with the local normal.
int loss(int,double,double,double,double); //Function to determine if a particle is lost during a wall collision. The first argument passes the region number specifying the model/
                     //parameters to utilize in the determination and the second argument passes the material potential to be used to check for wall penetration.
                     //RETURNS: -1 on error, 1 for wall penetration, 2 for physical wall loss, and 0 otherwise.
                    //last three arguments, normally zero, are used if you want to pass a predetermined normal
int depol(int); //Function to determine if a particle is depolarized during a collision. The argument passes the region number specifying the model/values.
               //RETURNS: -1 (and alters the particle's spin appropriately) if depolarization occurs, 0 for no depolarization, and -5 on error.
int propinmed(double,double [3][3]); //Function to handle propagation of a particle through a medium. The first argument which is passed is the time to the next
                                  //intersection point and the second is the particle's trajectory matrix. The function returns an event code (or -2 if there is no
                                 //interaction).
double intersection(int*,double [3][3],int *); //Function to take the particle's current dynamical information and find the earliest-time intersection with any
                                            //accessible region or cut-plane. It returns the elapsed time between the current location and the intersection.
                                           //The first pointer argument allows an eventcode to be returned. The second pointer argument returns -1 if the
                                          //intersection is with the current region. If the intersection is with a cut-plane then it returns the number of
                                         //the region associated with the cut-plane. The non-pointer argument passes a trajectory coefficient matrix
                                        //(see below) to the function.
double planesect(int,double [3][3]); //Function to find the intersection of a trajectory with a cut-plane.
                                  //First argument: Region # associated with the cut-plane to be checked for intersection.
                                 //Second argument: Trajectory coefficient matrix in the form
                                //                        cx2  cx1  cx0
                               //                         cy2  cy1  cy0
                              //                          cz2  cz1  cz0
                             //for a trajectory component equation c2 t^2 + c1 t + c0
                            //Function returns: Elapsed time from starting point to intersection.
double boxsect(int,double [3][3]); //Function to find the intersection of a trajectory with a box region.
                                //First argument: Region #
                               //Second argument: Trajectory coefficient matrix (see above).
                              //Function returns: Elapsed time from starting point to intersection.
double cylindersect(int,double [3][3]); //Function to find the intersection of a trajectory with a cylinder region.
                                     //First argument: Region #
                                    //Second argument: Trajectory coefficient matrix (see above).
                                   //Function returns: Elapsed time from starting point to intersection.
void normal(double*,double*,double*,double,double,double); //Function to calculate the (inward) normal vector to a point on a surface.
                                     //For a cut-plane this normal will always be given into the region associated with the cut-plane.
                                    //Pointer (returned) arguments: (nx,ny,nz)
                                   //Note that the function returns (0,0,0) in case of error.
                                  //last three arguments, normally zero, are used if you want to pass a predetermined normal - this is generally done when simulating surface roughness
void gsys2bsys(double,double,double,double*,double*,double*,double,double,double); //Function to transform the components of a vector into a system with the +x-dir along
                                                             //the inward normal at the current point of intersection, the +z-dir in the longitudinal
                                                            //direction, and the +y-dir consistent with a right-handed coordinte system.
                                                           //Arguments: (xi,yi,zi) | Returned (pointer) arguments: (xf,yf,zf).
                                                          //last three arguments, normally zero, are used if you want to pass a predetermined normal
void bsys2gsys(double,double,double,double*,double*,double*,double,double,double); //Function that is the inverse of 'gsys2bsys'.
                                                             //Arguments: (xi,yi,zi) | Returned (pointer) arguments: (xf,yf,zf).
                                                            //NB: The definition of the axes perpendicular to the normal may change from intersection point
                                                           //to intersection point, so 'gsys2bsys' and 'bsys2gsys' can only be used as inverses when applied
                                                          //at the same point in the geometry.
                                                        //last three arguments, normally zero, are used if you want to pass a predetermined normal
int move(double,double [3][3],int); //Function to update the particle's position, velocity, and time. The first argument is the total time during which it
                             //follows the trajectory from its original point. The second argument is the trajectory coefficient matrix (see above).
                            //The third argument specifies a checking code so that 'move' knows which geometry checks to perform.
                           //Geometry Checking Codes:  0 => Ensure that particle has not left a region.
                          //                           1 => Ensure that particle has not left a region AND check for a valid cut-plane intersection.
                         //The function returns 0 if everything is OK, and -1 if the particle is iredeemably outside of the geometry.
void trajmatrix(double [3][3]); //Function to calculate the trajectory matrix (see above) for a particle. The matrix is returned through the array argument.
int tjunc(int); //Function to determine whether the intersection with a region that is T'ed into happens at the location of the T-junction. The argument gives
               //the region number of the region that is T-ing into the current region.
              //Function returns: -1 -> Error ; 0 -> Intersection is NOT at T-junction ; 1 -> Intersection IS at T-junction.
void cplaneshift(int,int); //Function to shift the basepoint and/or length of a T-ing region depending on which side of the cut-plane the particle is currently
                          //moving. The first argument should be +1 to expand the T-ing region into the connecting region and -1 to shrink the T-ing region so that
                         //it no longer penetrates the connecting region. The second argumet passes the region number of the region to be changed.
int cplanehandling(void); //Function to take care of handling cut-plane intersections. It passes, bounces, etc. the particle in accordance with the handling code.
                         //The funtion returns -1 on error, 1 for a physical loss, and 0 otherwise.
void recdet(int [5000][10],int*); //Function to write the detector histograms to a data file. The function passes the array which holds the detector histograms (first
                           //argument) and a pointer to an array that holds the integrated counts for each detector (second argument).
int paramadjust(char [],char [],int,int,double); //Function to adjust a physics parameter to a new value. The arguments are as follows:
                                                //(action, parameter name, start region, end region, new parameter value) where 'start region' and 'end region'
                                               //define the range of regions for which the parameter should change, and the 'action' string is "set" or "reset"
                                              //where "set" indicates to load the passed value and "reset" indicates that the original value should be loaded.
                                             //The function returns 0 for OK and -1 on error.
void event(int,double,double,double,double); //Function to write an event into the 'events.sim' file. The first argument is the event code:
                                            //Event Codes:  -1 => Error
                                           //                0 => Record a trajectory point
                                          //                 1 => Intersection with a region
                                         //                  2 => Intersection with a cut-plane
                                        //                   3 => Leaving a specular bounce
                                       //                    4 => Leaving a non-specular bounce
                                      //                     5 => Leaving a specular bounce where particle was depolarized
                                     //                      6 => Leaving a non-specular bounce where particle was depolarized
                                    //                       7 => Loss of particle upon wall interaction other than due to penetration.
                                   //                        8 => Absorption of particle in a region
                                  //                         9 => Absorption of particle in a detector region.
                                 //                         10 => Spin-flip
                                //                          11 => Beta-decay of a neutron
                               //                           12 => Entering a bounce or checking for reflection off an increase in Fermi potential.
                              //                            13 => Particle lost due to penetration into wall or penetrates a bulk medium.
                             //                             14 => A created particle's dynamical information.
                            //                              15 => Particle gets an energy shift due to an interaction with a bulk medium.
                           //                               16 => Intersection with a T-junction.
                          //                                17 => Particle scattered inside a bulk medium.
						  //								18 => Detection of a particle by a cut-plane
                          //                                19 => Particle bounced off a cutplane with Fermi Gradient
                        //                                  20 => Particle penetrated a cutplane with Fermi Gradient
                         //
                        //The subsequent arguments allow data to be passed to 'event'. The type of data will depend on the event code.

int perturbednormal(double *nx,double *ny,double *nz); //finds a normal vector and then perturbs it so as to simulate surface roughness, created by Erik Lutz
                                                      //returns and integer indicating the model used

struct particle {  //This structure holds all the dynamical information about a particular particle.
  double t; //time
  double x; //global x coordinate
  double y; //global y coordinate
  double z; //global z coordinate
  double vx; //x-component of velocity (relative to global coordinate system)
  double vy; //y-component of velocity (relative to global coordinate system)
  double vz; //z-component of velocity (relative to global coordinate system)
  double spin[3]; //three components of a UNIT VECTOR pointing in the direction of the particle's spin (relative to global coordinate system)
  int region; //region number in which particle currently resides
  int xcode; //This variable is set to -1 if the current intersection of the particle is with the current region. It is set to -2 if the particle
            //has no current intersection (i.e. it is being propagated through a complex magnetic field). If the current intersection is with a
           //cut-plane then this is set to the region number associated with the cut-plane.
  int num; //particle number
};
struct particle neutron;
double bparamhold[100]; //Global variable to hold the original values of a parameter changed by a batch file.
double previntersectime; //Global variable to hold the previous result from 'solve4' so that iteration from that root via Newton's Method may be performed.

FILE *eventsfp; //File pointer for the file to which we will write out event information.
FILE *detectorsfp; //File pointer for the file to which we will write out detector information.

int main(int argc, char *argv[]) {
  int h,i,j,k,nbatch,nsim,startreg,endreg,simskip;
  int counts[10]; //Variables to receive integrated counts from detectors.
  double pvalue;
  char skip[500]; //Character array for parsing a 'batch' file.
  char regionsfn[100],connexfn[100],batchcountsfn[100],eventsfn[100],detectorsfn[100],pname[100];
  
  for(k=0 ; k < 10 ; k++) counts[k] = 0; //Zero the variables for recording integrated counts in the detectors during a simulation run in batch mode.
  if(argc == 1){ //no path specified as argument
    printf("specify the path for the region and connex files OR batch file as the argument. './' is okay.");
    return 0;
  }
  
  if(argc == 2 && BATCH == "OFF"){ //Path specified AND not using batch mode
    char* pathRegFile; //full path of the region file
    pathRegFile = malloc(strlen(argv[1])+1+10); /* make space for the new string (should check the return value ...) */
    strcpy(pathRegFile, argv[1]); /* copy name into the new var */
    strcat(pathRegFile, "RegionFile"); /* add the extension */

    char* pathConFile; //full path of the connex file
    pathConFile = malloc(strlen(argv[1])+1+10); /* make space for the new string (should check the return value ...) */
    strcpy(pathConFile, argv[1]); /* copy name into the new var */
    strcat(pathConFile, "Connexfile"); /* add the extension */
    
    char* geomoutSimFile; //full path of the geomout.sim file to output. Put it in the same folder as the Region and Connex files
    geomoutSimFile = malloc(strlen(argv[1])+1+11); /* make space for the new string (should check the return value ...) */
    strcpy(geomoutSimFile, argv[1]); /* copy name into the new var */
    strcat(geomoutSimFile, "geomout.sim"); /* add the extension */
    
    char* pathEventsSimFile; //full path of the events.sim file to output. Put it in the same folder as the Region and Connex files
    pathEventsSimFile = malloc(strlen(argv[1])+1+10); /* make space for the new string (should check the return value ...) */
    strcpy(pathEventsSimFile, argv[1]); /* copy name into the new var */
    strcat(pathEventsSimFile, "events.sim"); /* add the extension */
    
    char* pathDetectorsSimFile; //full path of the detectors.sim file to output. Put it in the same folder as the Region and Connex files
    pathDetectorsSimFile = malloc(strlen(argv[1])+1+10); /* make space for the new string (should check the return value ...) */
    strcpy(pathDetectorsSimFile, argv[1]); /* copy name into the new var */
    strcat(pathDetectorsSimFile, "detectors.sim"); /* add the extension */
    
    geomread(pathRegFile,pathConFile);
    printf("+Geometry data successfully read from '%s'.\n",argv[1]);
    geomprint(geomoutSimFile);
    simexec(pathEventsSimFile,pathDetectorsSimFile,counts);
    return 0;
  }
    
  if(argc == 2 && BATCH == "ON") {  //We must use the 'batch' file to run multiple simulations.  
    char* pathBatchFile; //full path of the batch file
    pathBatchFile = malloc(strlen(argv[1])+1+5); /* make space for the new string (should check the return value ...) */
    strcpy(pathBatchFile, argv[1]); /* copy name into the new var */
    strcat(pathBatchFile, "batch"); /* add the extension */
    
    char* pathEventsSimFile; //full path of the events.sim file to output. Put it in the same folder as the Region and Connex files
    pathEventsSimFile = malloc(strlen(argv[1])+1+10); /* make space for the new string (should check the return value ...) */
    strcpy(pathEventsSimFile, argv[1]); /* copy name into the new var */
    strcat(pathEventsSimFile, "events.sim"); /* add the extension */
    
    FILE *batchfp; //File pointer for the batch file.
    FILE *batchcountsfp; //File pointer for the integrated counts output file.
    batchfp = fopen(pathBatchFile,"r");

    fgets(skip,500,batchfp); //Skip "Number of BATCHES to run:".
    fscanf(batchfp,"%d\n",&nbatch); //Get the number of batches to run.
    fgets(skip,500,batchfp); //Skip "-".
    for(i=0 ; i < nbatch ; i++) {  //Loop to run through the required number of batches.
      fgets(skip,500,batchfp); //Skip "Number of simulations this batch:".
      fscanf(batchfp,"%d\n",&nsim); //Get the number of simulations required for this batch.
      fgets(skip,500,batchfp); //Skip "Geometry files this batch (regions connex):".
      fscanf(batchfp,"%s %s\n",regionsfn,connexfn); //Get the geometry file names for this batch.
      printf("___________________________________________________________\n");
      geomread(regionsfn,connexfn);  //Process the geometry files for this batch.
      geomprint(pathEventsSimFile);
      fgets(skip,500,batchfp); //Skip "Integrated counts filename for this batch:".
      fscanf(batchfp,"%s\n",batchcountsfn); //Get the filename for the integrated counts output file for this batch.
      batchcountsfp = fopen(batchcountsfn,"w+");
      fgets(skip,500,batchfp); //Skip batch information header.
      fprintf(batchcountsfp,"#Startreg | Endreg | Pname |  Pvalue  | D1 Counts | D2 Counts | D3 Counts | D4 Counts | D5 Counts\n"); //Write a header for 'counts'.
        for(j=0 ; j < nsim ; j++) {  //Loop to run through the required number of simulations for this batch.
          simskip = 0; //Set the batch file error flag.
          fscanf(batchfp,"%d %d %s %s %s %lf\n",&startreg,&endreg,eventsfn,detectorsfn,pname,&pvalue);  //Get batch information.
          for(k=0 ; k < 10 ; k++) counts[k] = 0; //Zero the variables for recording integrated counts in the detectors for a particular simulation run.
          printf("-------------------------------------------------------\n");
          printf("\nBatch# = %d  |  Simulation# = %d  |  %s = %f\n\n",i,j,pname,pvalue); //Print batch information to console.
          if(paramadjust("set",pname,startreg,endreg,pvalue) != 0) {  // Adjust the relevant parameters, terminate this simulation on error.
            printf("+ERROR: Qued Simulation Terminated... continuing to next simulation.\n");
            simskip = 1;
          }
          if(simskip == 0) simexec(eventsfn,detectorsfn,counts); //No error in reading the batch file so execute the simulation.
          if(paramadjust("reset",pname,startreg,endreg,pvalue) != 0) {  //Return the adjusted parameter to its original value or deal with error conditions.
            if(simskip != 1) {
              printf("+ERROR: Error in batch file not caught in initial parameter adjustment!\n");
              break;
            }
          }
          fprintf(batchcountsfp,"%d           %d        %s    %f",startreg,endreg,pname,pvalue);
          fprintf(batchcountsfp,"   %d           %d           %d           %d           %d\n",counts[0],counts[1],counts[2],counts[3],counts[4]);
        }
      fclose(batchcountsfp);
      fgets(skip,500,batchfp); //Skip the "-" that demarks batch runs.
      if(skip[0] != 45) {  //45 is the ASCII code for the character "-".
        printf("+ERROR: Batch File Processing Error! Terminating all subsequent batches...\n");
        return 0;
      }
    }
    fclose(batchfp);
    return 0;
  }
}

int simexec(char *eventsfn,char *detectorsfn,int *counts) {
  int reinteract;
  time_t randomseed;
  double timeelapsed;
  int timeestimate;
  time_t *calendartime; //time_t is variable type appropriate for reading the calendar time from the system.
  struct tm *loctime; //tm is a pre-defined structure for holding the local time.
  //loctime = localtime(calendartime); //Converts the calendar time into the local time tm-structure.
  char *timestamp;
  char skip[500];
  //timestamp = asctime(loctime); //Takes info. from tm-structure and makes a character string with date/time.
  
  if (RSEED == "RAND") randomseed = time(NULL); //choose a new seed for the random number generator each time
  
  srand(randomseed); //Set the random seed to 'RSEED' specified in header or to a random number each time if that is specified
  detectorsfp = fopen(detectorsfn,"w+"); //Open the 'detectors' file.
  eventsfp = fopen(eventsfn,"w+"); //Open the 'events' file.
  fprintf(eventsfp,"#Neutrons Propagated = %d | Random Seed = %ld\n",N,randomseed); //Put header into events file.
  //fprintf(eventsfp,"#Lost Particle Tally:      of the particles propagated through the geometry were lost unphysically.                           \n\n");
  fprintf(eventsfp,"#Neutron# | ecode | Reg# | xcode |     t     |     x     |     y     |     z     |     vx    |     vy    |     vz    |");
  fprintf(eventsfp,"   Spinz   |   Data1   |   Data2   |   Data3   |   Data4   |\n");
    
  //printf("\nBegin Simulation... %s\n\n",timestamp);
  printf("+There will be %d neutrons propagated through the geometry.\n",N);
  printf("+Random seed is %ld.\n\n\n",randomseed);
  
  if (CONSOLE == "ON") printf("\nBEGIN CONSOLE-----------------------------------------------------------------------------------\n");
  
//-------------------------------------------------BEGIN SIMULATION LOGIC--------------------------------------------------//
  int its,n,ecode,i,j,lost,bcode,cpcode;
  int detectors[5000][10]; //Array to hold detector hits. It currently allows .1s resolution over a 5min. run and provides slots for up to ten detectors.
  
  for(i=0 ; i < 5000 ; i++) {  //Zero the 'detectors' array.
    for(j=0 ; j < 10 ; j++) {
      detectors[i][j] = 0;
    }
  }
  lost = 0; //Zero the lost particle counter.
  for(n=0 ; n < N ; n++) {
    //New particle so zero dynamical variables and set codes etc. appropriately.
    neutron.t = grn()*BEAMTIME; //For timing info.
    // neutron.t = 0;
    neutron.x = 0;
    neutron.y = 0;
    neutron.z = 0;
    neutron.vx = 0;
    neutron.vy = 0;
    neutron.vz = 0;
    neutron.spin[0] = 0;
    neutron.spin[1] = 0;
    neutron.spin[2] = 0;
    neutron.region = -1;
    neutron.xcode = -5;
    neutron.num = n;
    if (neutron.num % 1000 == 0 && neutron.num > 0) {
        timeelapsed = difftime(time(NULL),randomseed);
        timeestimate = (int)((timeelapsed/n)*(N-n)/60);
        printf("Creating neutron %d. Approximately %d minutes remaining.\n",neutron.num,timeestimate);
    }
    poof(0,grn()*0.68,0,1); //Create the particle.
    //neutron.region = 0;  neutron.vz = 1.5;  neutron.vx = 0.;  neutron.vy = 0.;  neutron.x = 0.;  neutron.y = 0.;  neutron.z = 0.05;

  
    for(its = 0 ; its < 100000 ; its++) {  //Run a particle for at most 100,000 iterations.
	
	     if(neutron.t<=SHUTTERTIME){  //Close Shutter at Region 9.
           paramadjust("set","bpot",9,9,180);
           paramadjust("set","abs",9,9,10000);
           paramadjust("set","PM",9,9,1);
           }
 
         if(neutron.t>SHUTTERTIME){  //Open Shutter at Region 9.
           paramadjust("reset","bpot",9,9,180);
           paramadjust("reset","abs",9,9,10000);
           paramadjust("reset","PM",9,9,1);
           }

      ecode = propagate(); //Propagate the particle to its next intersection.
      if(ecode == -1) {  //This particle encountered an error during propagation so abandon it and move on to the next particle.
        lost++; //Increment the lost particle counter.
        break;
      }
      if(ecode == 8 || ecode == 11) break; //Particle was lost through a physical mechanism... move on to the next particle.
      if(ecode == 9) {  //Particle was absorbed in a detector region.
        if(DETCUTOFF > 0 && neutron.region == CUTOFFREG && sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)) > DETCUTOFF) break;
        counts[(int)rparams[neutron.region][10]-1]++; //Increment the integrated counts for the appropriate detector.
        for(i = 0 ; i < 5000 ; i++) {  //Find the correct bin to increment.
          if(neutron.t >= i/10. && neutron.t < (i+1.)/10.) detectors[i][(int)rparams[neutron.region][10]-1]++; //Increment the bin for the region's detector number.
        }
        break; //Particle has been detected so move on to the next particle.
      }
      if(ecode == 2 || ecode == 16) {  //Particle intersected a cut-plane or a T-junction.
        cpcode = cplanehandling(); //Handle the cut-plane intersection.
        reinteract=0;
          while (cpcode == 2) {
              cpcode = cplanehandling(); //particle was bounced but continues in same direction b\c of roughness, call cphandling again to simulate it rehitting surface
              reinteract++;
          }
        //if (CONSOLE == "ON" && reinteract != 0) printf("Particle re-interacted with rough surface %d times.\n",reinteract);
        if(cpcode == 1) break; //The particle was lost due to physical interactions so move on to the next particle.
        if(cpcode == 9) { //Particle was absorbed in a detector region.
            if(DETCUTOFF > 0 && neutron.region == CUTOFFREG && sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)) > DETCUTOFF) break;
            counts[(int)rparams[neutron.xcode][10]-1]++; //Increment the integrated counts for the appropriate detector.
            for(i = 0 ; i < 5000 ; i++) {  //Find the correct bin to increment.
              if(neutron.t >= i/10. && neutron.t < (i+1.)/10.) detectors[i][(int)rparams[neutron.xcode][10]-1]++; //Increment the bin for the region's detector number.
            }
            break;
        }
        if(cpcode == -1) {  //The particle was lost unphysically.
          lost++; //Increment the lost particle counter.
          break; //Move on to the next particle.
        }
      }
      if(ecode == 17) {}  //The particle scattered in a bulk medium so no bounce is required... the next geometry intersection point may be sought.
      if(ecode == 1) {  //Particle intersected with the current region.
        bcode = bounce(neutron.region,rparams[neutron.region][6],1,0,0,0); //Bounce the particle using models/values and wall potential of current region.
        if(bcode == -1) {  //Particle was lost unphysically.
          lost++; //Increment the lost particle counter.
          break; //Abandon the particle.
        }
        if(bcode == 1 || bcode == 2) break; //Particle was lost physically, so move on to the next particle.
      }
    }
  }
  recdet(detectors,counts); //Record the detectors' histograms in a data file.
//-------------------------------------------------END SIMULATION LOGIC----------------------------------------------------//
  
  if (CONSOLE == "ON") printf("END CONSOLE-----------------------------------------------------------------------------------");
  printf("\n\n");
  printf("+Final Lost Particle Tally: %f%% (%d) of the particles propagated through the geometry were lost unphysically during this simulation.\n",(float)lost/N*100.,lost);
  printf("\n\n");
  
  printf("Total Counts:\n");
  for(i=0 ; i < 10 ; i++) printf("Detector %d -> %d\n",(i+1),counts[i]);
  printf("\n\n");
  
  printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
  printf("##################################################################################################\n");
  
  rewind(eventsfp); //Back to the beginning of 'events' so we can record lost particle percentage.
  fgets(skip,500,eventsfp); //Skip header line already in place.
  //fprintf(eventsfp,"#Lost Particle Tally: %f%% of the particles propagated through the geometry were lost unphysically.",(float)lost/N*100.);
  fclose(eventsfp);
  fclose(detectorsfp);
  return 0;
}

void geomprint(char *geomoutfn) {
  int i,j;  
  FILE *geomoutfp;
  geomoutfp = fopen(geomoutfn,"w");
  
  fprintf(geomoutfp,"regions:\n\n");
  for(i=0 ; i < 100 ; i++) {
    if(i!=0 && (int)regions[i][0] == 0) break; //Stop writing geometry data when all regions have been recorded.
    fprintf(geomoutfp,"region# = %d | surface type = %f | ",i,regions[i][0]);
    fprintf(geomoutfp,"Dimx or Diameter = %f | Dimy or Length = %f | Dimz = %f | ",regions[i][1],regions[i][2],regions[i][3]);
    fprintf(geomoutfp,"psi = %f | theta = %f | phi = %f\n",regions[i][6],regions[i][4],regions[i][5]);
    fprintf(geomoutfp,"\n");
  }
  fprintf(geomoutfp,"\n");
  
  fprintf(geomoutfp,"basepoints:\n\n");
  for(i=0 ; i < 100 ; i++) {
    if(i!=0 && (int)regions[i][0] == 0) break; //Stop writing geometry data when all regions have been recorded.
    fprintf(geomoutfp,"region# %d --> bpx = %f | bpy = %f | bpz = %f",i,basepoints[i][0],basepoints[i][1],basepoints[i][2]);
    fprintf(geomoutfp,"\n\n");
  }
  fprintf(geomoutfp,"\n");
  
  fprintf(geomoutfp,"cplanes:\n\n");
  for(i=0 ; i < 100 ; i++) {
    if(i!=0 && (int)regions[i][0] == 0) break; //Stop writing geometry data when all regions have been recorded.
    fprintf(geomoutfp,"region # = %d | psi = %f | theta = %f | phi = %f | s.h. code = %d | ",i,cplanes[i][7],cplanes[i][0],cplanes[i][1],(int)cplanes[i][2]);
    fprintf(geomoutfp,"int. s.h. code = %d | T cut-plane offset = %f | s.h. 2 value = %d\n",(int)cplanes[i][4],cplanes[i][5],(int)cplanes[i][3]);
    fprintf(geomoutfp,"\n");
  }
  fprintf(geomoutfp,"\n");
  
  fprintf(geomoutfp,"rparams:\n\n");
  for(i=0 ; i < 100 ; i++) {
    if(i!=0 && (int)regions[i][0] == 0) break; //Stop writing geometry data when all regions have been recorded.
    fprintf(geomoutfp,"%f|%f|%f|%f|%e|%e|%f|%f|",rparams[i][0],rparams[i][1],rparams[i][2],rparams[i][3],rparams[i][4],rparams[i][5],rparams[i][6],rparams[i][7]);
    fprintf(geomoutfp,"%e|%e|%f|%f|%f|%f|%f",rparams[i][8],rparams[i][9],rparams[i][10],rparams[i][11],rparams[i][12],rparams[i][13],rparams[i][14]);
    fprintf(geomoutfp,"\n\n");
  }
  fprintf(geomoutfp,"\n");
  
  fprintf(geomoutfp,"connex:\n\n");
  for(i=0 ; i < 100 ; i++) {
    if(i!=0 && (int)regions[i][0] == 0) break; //Stop writing geometry data when all regions have been recorded.
    fprintf(geomoutfp,"region %d --> %d and %d and %d and %d and %d and %d",i,connex[i][0],connex[i][1],connex[i][2],connex[i][3],connex[i][4],connex[i][5]);
    fprintf(geomoutfp,"\n\n");
  }
  fprintf(geomoutfp,"\n");
  
  fclose(geomoutfp);

  //printf("+Geometry data written to 'geomout.sim'.\n");
  return;
}

void geomread(char regionsfile[20],char connexfile[20]) {
  int i,j,k,l,reg,prevreg,surftype,prevsurftype,chflag,shcode,conreg,regmain;
  int regnum; //The number of regions defined in 'regions'.
  int conchk[100],bpflag,regflag;
  char skip[500],parse[2];
  long int fpos;
  double bpx,bpy,bpz,psir,thetar,phir,thetarcp,psircp,d1,d2,alpha,beta;
  double dbdxg,dbdyg,dbdzg;
  double bppx,bppy,bppz;
  
  fpos_t *posp; //fpos_t is the type for a file pointer and is used below to reposition pointer in 'regions'.
  
  FILE *regionsfp,*connexfp;
  regionsfp = fopen(regionsfile,"r");
  connexfp = fopen(connexfile,"r");
  
//Clear all geometry arrays  
  for(i=0 ; i < 100 ; i++) {
    for(j=0 ; j < 6 ; j++) {
      regions[i][j]=0;
    }
  }
  for(i=0 ; i < 100 ; i++) {
    for(j=0 ; j < 7 ; j++) {
      cplanes[i][j]=0;
    }
  }  
  for(i=0 ; i < 100 ; i++) {
      for(j=0 ; j < 3 ; j++) {
        basepoints[i][j]=0;
      }
    }
  for(i=0 ; i < 100 ; i++) {
    for(j=0 ; j < 6 ; j++) {
      connex[i][j]=0;
    }
  }
  for(i=0 ; i < 100 ; i++) {
    for(j=0 ; j < 16 ; j++) {
      rparams[i][j]=0;
    }
  }
  
  regnum = 0;
  
  for(i=0 ; i < 100 ; i++) conchk[i]=0; //This array keeps track below of which region's cut-planes have already been defined, and whether the cut-plane's
                                      //orientation needs to be calculated (if basepoint is inferred) or given the orientation of the region (if
                                     //basepoint is explicitly specified, in which case a T is assumed and the basepoint appropriately altered).
                                    //[-1 => cut-plane has been defined | -2 => Basepint explicitly defined for a region other than 0]
  
//Determine how many regions are defined.
  regflag = 0;
  for(i=0 ; i < 100 ; i++) {
    fgets(skip,500,regionsfp);
    if(skip[0] == '/') break;
    if(feof(regionsfp) != 0) {
      regflag = -1;
      printf("\n+WARNING: No end-of-file flag detected in 'regions' file. Check number of defined regions!\n");
      break;
    }
  }
  if(regflag == -1) regnum = i;
  else regnum = i-1;
  printf("\n\n##################################################################################################\n");
  printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
  printf("\n\nUCN TRANSPORT SIMULATION\n\n+There are %d regions defined.\n",regnum);
  rewind(regionsfp);

//Fill the 'rparams' array and place the region's dimensions in 'regions'.
//Key to 'rparams':
   //0 -> dB/dx ; 1 -> dB/dy ; 2 -> dB/dz ; 3 -> Spec ; 4 -> Losspb ; 5 -> Depol. ; 6 -> Wall Pot. ; 7 -> Bulk Pot.
   //8 -> Surf/Scatt ; 9 -> Absorp. ; 10 -> Det ; 11 -> Prop Mode ; 12 -> Scatt Model ; 13 -> Loss Model ; 14 -> Depol Model
  fgets(skip,500,regionsfp); //Skip the header line.
  for(i=0 ; i < regnum ; i++) {
    fscanf(regionsfp,"%d %d %s",&reg,&surftype,skip); //Get region number and surface type and skip over to the region's dimesnions.
    if((surftype != 1) && (surftype != 2) && (surftype != 3)) printf("+ERROR: Undefined Region Type!\n");
    regions[reg][0] = (double)surftype; //Place this region's surface type into 'regions'.
    if(surftype == 1) {
      fscanf(regionsfp,"%lf,%lf,%lf %s",&regions[reg][1],&regions[reg][2],&regions[reg][3],skip);
      fscanf(regionsfp,"%lf,%lf,%lf %lf",&rparams[reg][0],&rparams[reg][1],&rparams[reg][2],&rparams[reg][3]);
      fscanf(regionsfp,"%lf %lf %lf %lf",&rparams[reg][4],&rparams[reg][5],&rparams[reg][6],&rparams[reg][7]);
      fscanf(regionsfp,"%lf %lf %lf %lf %lf",&rparams[reg][8],&rparams[reg][9],&rparams[reg][10],&rparams[reg][11],&rparams[reg][12]);
      fscanf(regionsfp,"%lf %lf",&rparams[reg][13],&rparams[reg][14]);
    }
    if(surftype == 2) {
      regions[reg][3] = 0; //There is no Dimz entry for a cylinder-region.
      fscanf(regionsfp,"%lf,%lf %s",&regions[reg][1],&regions[reg][2],skip);
      fscanf(regionsfp,"%lf,%lf,%lf %lf",&rparams[reg][0],&rparams[reg][1],&rparams[reg][2],&rparams[reg][3]);
      fscanf(regionsfp,"%lf %lf %lf %lf",&rparams[reg][4],&rparams[reg][5],&rparams[reg][6],&rparams[reg][7]);
      fscanf(regionsfp,"%lf %lf %lf %lf %lf",&rparams[reg][8],&rparams[reg][9],&rparams[reg][10],&rparams[reg][11],&rparams[reg][12]);
      fscanf(regionsfp,"%lf %lf",&rparams[reg][13],&rparams[reg][14]);
    }
  }
  rewind(regionsfp);

//Place the region orientation information in 'regions', convert DEGREES to RADIANS, and transform input magnetic field gradients into the global system.
  fgets(skip,500,regionsfp); //Skip header line.
  for(i=0 ; i < regnum ; i++) {
    fscanf(regionsfp,"%d %d %s %s",&reg,&surftype,skip,skip); //Get region number and skip over to orientation information.
    fscanf(regionsfp,"%lf,%lf,%lf",&regions[reg][6],&regions[reg][4],&regions[reg][5]);
    regions[reg][4] = PI/180 * regions[reg][4]; //Convert theta from DEGREES to RADIANS.
    regions[reg][5] = PI/180 * regions[reg][5]; //Convert phi from DEGREES to RADIANS.
    regions[reg][6] = PI/180 * regions[reg][6]; //Convert psi from DEGREES to RADIANS.
    fgets(skip,500,regionsfp); //Skip to the next line.
    rota(regions[reg][6],regions[reg][4],regions[reg][5],rparams[reg][0],rparams[reg][1],rparams[reg][2],&dbdxg,&dbdyg,&dbdzg); //Transform the given (solving
                                                                                                                               //system) magnetic field gradients
                                                                                                                              //into the global system. (The
                                                                                                                             //trajectory matrix is formed
                                                                                                                            //assuming everything is relative to the
                                                                                                                           //global system.)
    rparams[reg][0] = dbdxg; //Store the global x field gradient component.
    rparams[reg][1] = dbdyg; //Store the global y field gradient component.
    rparams[reg][2] = dbdzg; //Store the global z field gradient component.
  }
  rewind(regionsfp);

//Fill the 'basepoints' array. (basepoints[][0] => x / basepoints[][1] => y / basepoints[][2] => z)
  fgets(skip,500,regionsfp); //Skip header line.
  reg = 0;
  surftype = 0;
  for(i=0 ; i < regnum ; i++) {
    bpflag = 0; //Reset the flag for determining how to interpret an explicitly defined basepoint.
    prevreg = reg; //Store the previous region number in case a basepoint needs to be inferred.
    prevsurftype = surftype; //Store the previous surface type in case a basepoint needs to be inferred.
    fscanf(regionsfp,"%d %d",&reg,&surftype); //Get new region # and surface type.
    fpos=ftell(regionsfp); //Read current position in 'regions' into fpos.
    fscanf(regionsfp,"%s",skip);
    if(*skip == '*' || *skip == '<') {  //'skip', being a character array, points to the memory location of the first element of the array.
                                       //So, '*skip' indicates the contents of that memory location, i.e. the first character of the character
                                      //string contained in 'skip'.
      bpflag = 1; //Set the flag indicating that this basepoint has been determined.
      if(*skip == '<') conchk[reg] = -4; //This region is the end of a contiguous path and T's into the side of another set of regions so its basepoint may be
                                        //deduced, but we must indicate that it is a T-ing region so that the appropriate codes will be set below.
      thetar = regions[prevreg][4]; //The new region's basepoint is at the end of the previous region's so that it must be displaced along the
      phir = regions[prevreg][5];  //orientation of the previous region.
      psir = regions[prevreg][6];
      if(prevsurftype==1) rota(psir,thetar,phir,0.,0.,regions[prevreg][3],&bpx,&bpy,&bpz);  //Rotate a displacement in the z-direction with the previous region's
      if(prevsurftype==2) rota(psir,thetar,phir,0.,0.,regions[prevreg][2],&bpx,&bpy,&bpz); //length into the orientation of the previous region.
      basepoints[reg][0] = basepoints[prevreg][0] + bpx;  //The new region's basepoint is obtained be displacing from the previous region's basepoint through
      basepoints[reg][1] = basepoints[prevreg][1] + bpy; //the rotated displacement just calculated.
      basepoints[reg][2] = basepoints[prevreg][2] + bpz;
    }
    if(*skip == '>') {  //The explicit basepoint defined occurs at a T-junction where the main region comes BEFORE (increasing z) the T-ing region's basepoint, i.e
                       //this region is the start of a new set of contiguous regions, and it begins in the side of a separate set of contiguous regions.
      bpflag = 1; //Set the flag indicating that this basepoint has been determined.
      fseek(regionsfp,fpos,SEEK_SET); //Backup to point in file before the basepoint was read.
      for(l=1 ; l<100 ; l++) {  //We must move forward in the file until we have gotten past the > symbol.
        fgets(parse,2,regionsfp);
        if(*parse == '>') break;
        if(l > 95) printf("+ERROR: Parsing for '>' symbol during determination of basepoints failed!\n");
      }
      fscanf(regionsfp,"%lf,%lf,%lf",&basepoints[reg][0],&basepoints[reg][1],&basepoints[reg][2]); //Read in the basepoint.
      conchk[reg] = -3; //Set the flag indicating a main region into T-junction.
    }
    if(bpflag != 1) {  //The explicitly defined basepoint does NOT indicate a T-juntion.
      conchk[reg] = -2; //Set the flag indicating an explicitly defined basepoint NOT associated with a T-juntion.
      fseek(regionsfp,fpos,SEEK_SET); //Backup to point in file before basepoint was read.
      fscanf(regionsfp,"%lf,%lf,%lf",&basepoints[reg][0],&basepoints[reg][1],&basepoints[reg][2]); //Read in the basepoint.
    }
    fgets(skip,500,regionsfp); //skip to the next line in 'regions' file.
  }
  rewind(regionsfp);
  
//Fill the 'connex' array.
  fgets(skip,500,connexfp); //Skip header line #1.
  fgets(skip,500,connexfp); //Skip header line #2.
  fgets(skip,500,connexfp); //Skip header line #3.
  for(i=0 ; i < regnum ; i++) {
    //Get the region number and parse for a special-handling code.
    fpos = ftell(connexfp); //Read the current position in 'connex' into 'fpos'.
    chflag = 0; //Clear the flag used to parse the inclusion of a special-handling identifier.
    fscanf(connexfp,"%s",skip);
    for(k=0 ; k<10 ; k++) {  //Scan through the first ten characters in 'skip'.
      if(skip[k] == ',') {  //There is a special-handling code.
        fseek(connexfp,fpos,SEEK_SET); //Backup to the point in the file before the string was read.
        fscanf(connexfp,"%d,%d",&reg,&shcode); //Get the current region number whose connection list we are parsing and special-handling code.
        cplanes[reg][2] = (double)shcode; //Load the special-handling code into 'cplanes'.
        if(cplanes[reg][2] == 2. || cplanes[reg][2] == 7. || cplanes[reg][2] == 9.) fscanf(connexfp,"(%lf)",&cplanes[reg][3]); //If the cut-plane has special-handling code 2,7,9 then read in the associated value.
        chflag = -1; //Set flag to indicate that a special-handling flag has been read.
        break;
      }
    }
    if(chflag != -1) {  //No special-handling flag was read so we only need to read in the region number.
      fseek(connexfp,fpos,SEEK_SET); //Backup to the point in the file before the string was read.
      fscanf(connexfp,"%d",&reg); //Get the current region number whose connection list we are parsing.
    }
    for(j=0 ; j < 6 ; j++) {  //A region can connect to at most six other regions.
      fpos=ftell(connexfp); //Read current position in 'connex' into fpos.
      fscanf(connexfp,"%s",skip); //Read in the next set of characters.
      if(skip[0] == '/') {  //We have reached the end of the list of connections for the current region.
        connex[reg][j] = -1; //This is the flag in 'connex' used to indicate the end of the connections list.
        break;
      }
      fseek(connexfp,fpos,SEEK_SET); //Backup to point in file before information was read since we aren't at the end of the line and must read in a new connection.
      fscanf(connexfp,"%d",&connex[reg][j]); //Load the connection into 'connex'.
      if(chflag != -1) cplanes[reg][2] = -1; //Indicate no special-handling only if no special-handling code has been set.
    }
    fgets(skip,500,connexfp); //Skip to the next line.
  }
  rewind(connexfp);
  
//Calculate cut-plane orientations and put them in 'cplanes'. Also find basepoints offsets for T connections and set appropriate internal special handling codes.
  fgets(skip,500,connexfp); //Skip header line #1.
  fgets(skip,500,connexfp); //Skip header line #2.
  fgets(skip,500,connexfp); //Skip header line #3.
  for(i = 0 ; i < regnum ; i++) {
    fscanf(connexfp,"%d",&reg); //Get the next region number.
    fgets(skip,500,connexfp); //Skip to the next line.
    if(reg == 0) {  //Region 0 is the start region so that its cut-plane does not depend on an adjoining region.
      cplanes[reg][0] = regions[reg][4]; //Make the cut-plane's theta-orientation the same as the region's.
      cplanes[reg][1] = regions[reg][5]; //Make the cut-plane's phi-orientation the same as the region's.
      cplanes[reg][7] = regions[reg][6]; //Make the cut-plane's psi-orientation the same as the region's.
      conchk[reg] = -1; //Indicate that region 0 has had its cut-plane defined.
    }
    for(j=0 ; j < 6 ; j++) {
      conreg = connex[reg][j]; //Region for which we wish to define a cut-plane.
      if(conreg == -1) break; //We have reached the end of this region's connection list.
      surftype = (int)regions[conreg][0]; //Get the surface type for the connecting region.
      if((int)regions[reg][0] == 2 && surftype == 1 && cplanes[conreg][4] == 0) cplanes[conreg][4] = 1;   //Set the internal special-handling flag for a cylinder
                                                                                                         //region connecting to a box region.
      if((int)regions[reg][0] == 1 && surftype == 2 && cplanes[conreg][4] == 0) cplanes[conreg][4] = 1; //or box-region connecting to a cylinder region.
      
      //Adjust the special-handling codes and store proper cut-plane offsets for T-juntions.
      if(conchk[conreg] == -3 || conchk[conreg] == -4) {  //We must calculate offsets and set internal special-handling flags for a T connection.
        cplanes[conreg][0] = regions[conreg][4]; //Make the cut-plane's theta-orientation the same as its region's.
        cplanes[conreg][1] = regions[conreg][5]; //Make the cut-plane's phi-orientation the same as its region's.
        cplanes[conreg][7] = regions[conreg][6]; //Make the cut-plane's psi-orientation the same as its region's.
        if(conchk[conreg] == -3) {   //Main Region -> T-ing Region so the region whose connections we are scanning is the main region.
          regmain = reg;
          cplanes[conreg][4] = 3; //Set the appropriate flag for the T-ing region.
        }
        if(conchk[conreg] == -4) {  //T-ing Region -> Main Region so the region that this T-ing region connects to NOT THROUGH ITS OWN CUT-PLANE
                                   //(so the second 'connex' entry) must be the main region.
          regmain = connex[conreg][1]; 
          cplanes[conreg][4] = 4; //Set the appropriate flag for the T-ing region.
        }
        if(connex[conreg][1] < 0 || connex[conreg][1] > 100) printf("+ERROR: T-ing region has no second entry in 'connex'!\n");
        cplanes[regmain][4] = 2; //Set the appropriate flag for the main region.
        if((int)regions[conreg][0] == 1 && (int)regions[regmain][0] == 2) {  //A box-region T's into a cylinder-region.
          if(pow(regions[regmain][1]/2.,2) - pow(regions[conreg][2]/2.,2) < 0) {  //Attempt to connect a smaller region into a larger region.
                                                                                 //NB: Here box's y-dim (height) is used to calculate the basepoint adjustment.
            printf("+ERROR: T-intersection defined with a box T'ed into a smaller cylinder! (Region #%d)\n",conreg);
          }
          else {
            cplanes[conreg][5] = regions[regmain][1]/2. - sqrt(pow(regions[regmain][1]/2.,2) - pow(regions[conreg][2]/2.,2)); //Store the T-junction offset.
            cplaneshift(1,conreg); //T-ing regions need to be expanded into the main region.
          }
        }
        if((int)regions[conreg][0] == 2 && (int)regions[regmain][0] == 1) {  //A cylinder-region T's into a box region.
          cplanes[conreg][5] = 1e-5; //Set a small offset to avoid simultaneous intersections with cut-plane and geometry.
          cplaneshift(1,conreg); //T-ing regions need to be expanded into the main region.
        }
        if((int)regions[conreg][0] == 1 && (int)regions[regmain][0] == 1) {  //A box-region T's into a box-region.
          cplanes[conreg][5] = 1e-5; //Set a small offset to avoid simultaneous intersections with cut-plane and geometry.
          cplaneshift(1,conreg); //T-ing regions need to be expanded into the main region.
        }
        if((int)regions[conreg][0] == 2 && (int)regions[regmain][0] == 2) {  //A cylinder-region T's ino a cylinder region.
          if(pow(regions[regmain][1]/2.,2) - pow(regions[conreg][1]/2.,2) < 0) {  //Attempt to connect a smaller region into a larger region.
            printf("+ERROR: T-intersection defined with a larger cylinder T'ed into a smaller cylinder! (Region #%d)\n",conreg);
          }
          else {
            cplanes[conreg][5] = regions[regmain][1]/2. - sqrt(pow(regions[regmain][1]/2.,2) - pow(regions[conreg][1]/2.,2)); //Store cut-plane offset.
            cplaneshift(1,conreg); //T-ing regions need to be expanded into the main region.
          }
        }
        conchk[conreg] = -1; //Set the connection check flag to indicate that we are finished with this region's cut-plane.
      }
      
      if(conchk[conreg] == -2) {  //The region has an explicitly defined basepoint not associated with a T-junction so that its cut-plane has the
                                 //same orientation as the region.
        cplanes[conreg][0] = regions[conreg][4]; //Make the cut-plane's theta-orientation the same as the region's.
        cplanes[conreg][1] = regions[conreg][5]; //Make the cut-plane's phi-orientation the same as the region's.
        cplanes[conreg][7] = regions[conreg][6]; //Make the cut-plane's psi-orientation the same as the region's.
        conchk[conreg] = -1; //Indicate that this region has had its cut-plane defined.
      }
      
      if(conchk[conreg] == 0 && abs(conreg-reg) == 1) {  //This region connects to a previous region and so we must calculate the appropriate cut-plane orientation.
                                                        //Also, we only want to calculate cases where the regions really do join together end-to-end!
        if(surftype == 1 || cplanes[conreg][4] == 1) { //Here we either have a box-region connected to a box-region or a cylinder-region connected to a box-
                                                      //region. In either case the region should have the same orientation as the previous region and there should be
                                                     //no orientation differences between the two regions.
          cplanes[conreg][0] = regions[conreg][4]; //The cut-plane for a box region has the same orientation as the region.
          cplanes[conreg][1] = regions[conreg][5];
          cplanes[conreg][7] = regions[conreg][6];
          conchk[conreg] = -1; //Mark this region as having had its cut-plane defined.
          //Check for problems.
          thetar = regions[conreg][4] - regions[reg][4]; //delta theta between the connecting regions.
          psir = regions[conreg][6] - regions[reg][6]; //delta psi between the connecting regions.
          if(thetar != 0. || psir != 0.) printf("+ERROR: Regions %d and %d have different orientations but at least one is a box-region!\n",reg,conreg);
        }
        if(surftype == 2 && cplanes[conreg][4] != 1) {  //Two cylinder-regions are connected.
          conchk[conreg] = -1; //Mark this region as having its cut-plane defined.
          thetar = regions[conreg][4] - regions[reg][4]; //delta theta between the connecting regions.
          psir = regions[conreg][6] - regions[reg][6]; //delta psi between the connecting regions.
          d1 = regions[reg][1]; //Diameter of the 'ending' region.
          d2 = regions[conreg][1]; //Diameter of the 'beginning' region.
          if(thetar == 0. && psir == 0. && d1 != d2) {  //In this case special handling will be required, but cut-plane should still have orientation of its region.
            cplanes[conreg][0] = regions[conreg][4];
            cplanes[conreg][1] = regions[conreg][5];
            cplanes[conreg][7] = regions[conreg][6];
          }
          if(thetar == 0. && d1 == d2) { //Formula is undefined in this case, but the cut-plane should have the theta-orientation of its region.
            cplanes[conreg][0] = regions[conreg][4];
            cplanes[conreg][1] = regions[conreg][5];
          }
          if(psir == 0. && d1 == d2) { //Formula is undefined in this case, but the cut-plane should have the psi-orientation of its region.
            cplanes[conreg][7] = regions[conreg][6];
            cplanes[conreg][1] = regions[conreg][5];
          }
          if(thetar != 0.) {  //Calculate the theta-orientation of the cut-plane. The formula gives the angle from the +z-axis (positive->LHR). If the 'beginning'
                             //region has a non-zero orientation, then the angle from the formula (alpha) must be subtracted from it to give the cut-plane's
                            //orientation relative to the +z-axis (positive->RHR). However, cut-planes with theta=0 are defined perpendicular to the z-axis, not
                           //parallel to it, so the final orientation of the cut-plane must be shifted by PI/2.
            alpha = atan(d1*sin(thetar)/(d2 - d1*cos(thetar)));
            if(alpha < 0.) thetarcp = regions[reg][4] - alpha - PI/2.; //In this case the cut-plane (when observed from an orientation where the first region
                                                                      //is horizontal) should be rotated counterclockwise so that its orientation (before
                                                                     //correcting for the initial region's orientation) should be made negative.
            else thetarcp = regions[reg][4] - alpha + PI/2.;
            cplanes[conreg][0] = thetarcp; //Save the cut-plane's theta-orientation in 'cplanes'.
            cplanes[conreg][1] = regions[conreg][5]; //Store the cut-plane's phi-orientation in 'cplanes' (same as region's.)
          }
          if(psir != 0.) {  //Calculate the psi-orientation of the cut-plane. The formula gives the angle from the +z-axis (positive->LHR). If the 'beginning'
                           //region has a non-zero orientation, then the angle from the formula (alpha) must be subtracted from it to give the cut-plane's
                          //orientation relative to the +z-axis (positive->RHR). However, cut-planes with psi=0 are defined perpendicular to the z-axis, not
                         //parallel to it, so the final orientation of the cut-plane must be shifted by PI/2.
            beta = atan(d1*sin(psir)/(d2 - d1*cos(psir)));
            if(beta < 0.) psircp = regions[reg][6] - beta - PI/2.; //In this case the cut-plane (when observed from an orientation where the first region
                                                                  //is horizontal) should be rotated counterclockwise so that its orientation (before
                                                                 //correcting for the initial region's orientation) should be made negative.
            else psircp = regions[reg][6] - beta + PI/2.;
            cplanes[conreg][7] = psircp; //Save the cut-plane's psi-orientation in 'cplanes'.
            cplanes[conreg][1] = regions[conreg][5]; //Store the cut-plane's phi-orientation in 'cplanes' (same as region's.)
          }
        }
        if(conchk[conreg] != -1) printf("+ERROR: Cut-plane for region %d not defined!\n",conreg);
      }
    }
    if(cplanes[reg][2] == 9.) cplanes[reg][0] = cplanes[reg][0] + (PI/180 * cplanes[reg][3]);
  }
  rewind(connexfp);
  
  fclose(regionsfp);
  fclose(connexfp);
  //printf("+Geometry data successfully read from '%s' and '%s'.\n",regionsfile,connexfile);
  
  return;
}

void rota(double psir,double thetar,double phir,double xi,double yi,double zi,double *xf,double *yf,double *zf) {
  *xf = xi*cos(phir)*cos(psir)+yi*(-sin(phir)*cos(thetar)*cos(psir)+sin(thetar)*sin(psir))+zi*(sin(phir)*sin(thetar)*cos(psir)+cos(thetar)*sin(psir));
  *yf = xi*sin(phir)+yi*cos(phir)*cos(thetar)-zi*cos(phir)*sin(thetar);
  *zf = -xi*cos(phir)*sin(psir)+yi*(sin(thetar)*cos(psir)+cos(thetar)*sin(phir)*sin(psir))+zi*(cos(thetar)*cos(psir)-sin(thetar)*sin(phir)*sin(psir));
  return;
}

void rotp(double psir,double thetar,double phir,double xi,double yi,double zi,double *xf,double *yf,double *zf) {
  *xf = xi*cos(phir)*cos(psir)+yi*sin(phir)-zi*cos(phir)*sin(psir);
  *yf = xi*(-cos(psir)*sin(phir)*cos(thetar)+sin(thetar)*sin(psir))+yi*cos(phir)*cos(thetar)+zi*(sin(thetar)*cos(psir)+cos(thetar)*sin(phir)*sin(psir));
  *zf = xi*(sin(thetar)*sin(phir)*cos(psir)+cos(thetar)*sin(psir))-yi*sin(thetar)*cos(phir)+zi*(cos(thetar)*cos(psir)-sin(thetar)*sin(phir)*sin(psir));
  return;
}

double solve2(double c2, double c1, double c0) {
  double D,t1,t2;

  //If c2=0 then we just need to solve a linear equation.
  if(c2 == 0) {
    if(c1 != 0) {  //Particle is indeed moving so that solution will be defined.
      if(CHECK == "ON") printf(">>solve2 linear equation solution: %f\n",-c0/c1);
      if(timezero(-c0/c1) > 0) return -c0/c1; //Solution is positive and gives a physically reasonable time.
      else return 0; //Solution is negative or zero so return zero.
    }
    else {
      if(CONSOLE == "ON") printf("+WARNING: Static Particle detected in 'solve2'!\n");
      return 0;
    }
  }

  //Calculate the discriminant of the quadratic intersection equation.
  D = pow(c1,2) - 4.*c2*c0;
  if(D < 0) {  //No real solutions so return a 0.
    if(CHECK == "ON") printf(">>No real solutions from solve2.\n");
    return 0;
  }

  //Calculate the two intersectiion times.
  t1 = (-c1 + sqrt(D))/(2.*c2);
  t2 = (-c1 - sqrt(D))/(2.*c2);
  if(CHECK == "ON") printf(">>Real solve2 solutions: t1 = %e | t2 = %e\n",t1,t2);

  //Choose the smallest non-zero time to return, where note that 'non-zero' means outside of the zero-boundary set by 'zero'.
  if(timezero(t1) > 0 && timezero(t2) <= 0) return t1;
  if(timezero(t2) > 0 && timezero(t1) <= 0) return t2;
  if(t1 <= t2 && timezero(t1) > 0) return t1;
  if(t2 <= t1 && timezero(t2) > 0) return t2;
  return 0; //Neither root is positive so return 0.
}

double solve4(double c4, double c3, double c2, double c1, double c0) {
  int i,flag1,flag2;
  double a3,a2,a1,a0;
  double u2,u1,u0;
  double D1,D2,D1same,D2same,D1opp,D2opp,r1,r2,r3,r;
  double bn,cn,bp,cp,bc,cc;
  double sdp,sdn;
  double t1,t2,t3,t4,thold,t,error;
  double t1same,t2same,t3same,t4same,t1opp,t2opp,t3opp,t4opp,told;

  //Two of the t1, t2, t3, t4 variables might not be set and yet enter comparsions below so we should zero them.
  t1=0; t1same=0; t1opp=0;
  t2=0; t2same=0; t2opp=0;
  t3=0; t3same=0; t3opp=0;
  t4=0; t4same=0; t4opp=0;
    
  if(CHECK == "ON") printf(">>Quartic Equation: %e t^4 + %e t^3 + %e t^2 + %e t + %e = 0\n",c4,c3,c2,c1,c0);
  
//See Abramowitz and Stegun p.17-18 for the quartic formula.
  if(c4 == 0) {  //'solve4' was called with no quartic term!
    printf("+ERROR: 'solve4' was called for a quartic with no fourth-order term!\n");
    return 0;
  }
  
  //'Normalize' quartic equation.
  a3 = c3/c4;
  a2 = c2/c4;
  a1 = c1/c4;
  a0 = c0/c4;
  
  //Calculate coefficients of the associated cubic equation.
  u2 = -a2;
  u1 = a1*a3 - 4.*a0;
  u0 = -(pow(a1,2) + a0*pow(a3,2) - 4.*a0*a2);

  solve3(1.,u2,u1,u0,&r1,&r2,&r3);
  if(CHECK == "ON") printf(">>Associated Cubic Equation: t^3 + %e t^2 + %e t + %e = 0\n",u2,u1,u0);
  if(CHECK == "ON") printf(">>Associated Cubic Equation Roots: r1 = %e | r2 = %e | r3 = %e\n",r1,r2,r3);

  if(r2 == -999999 && r3 == -999999) r = r1; //Only one real root from the cubic equation.
  else {  //All three roots of the cubic equation are real.
    //Choose the largest root -- in case round-off makes a zero into a small non-zero value -- that makes the coefficients of the quadratic equation real.
    bc = pow(a3,2)/4. + r1 - a2;
    cc = pow(r1/2,2) - a0;
    if(bc >= 0 && cc >= 0) r = r1;
    bc = pow(a3,2)/4. + r2 - a2;
    cc = pow(r2/2,2) - a0;
    if(bc >= 0 && cc >= 0 && fabs(r2) > fabs(r)) r = r2;
    bc = pow(a3,2)/4. + r3 - a2;
    cc = pow(r3/2,2) - a0; 
    if(bc >= 0 && cc >= 0 && fabs(r3) > fabs(r)) r = r3;
  }
  if(CHECK == "ON") printf(">>Root Chosen from Associated Cubic: r = %e\n",r);

  //Set the coefficients to the real values that were determined.
  bc = pow(a3,2)/4. + r - a2;
  cc = pow(r/2.,2) - a0;
  
  if(bc < 0) {  //It is possible that 'solve3' returns one real root which (due to roundoff) produces a negative value for 'bc'.
    bc = mathzero(bc);
    if(bc != 0) {
      printf("+ERROR: Imaginary 'bc' coefficient for quadratic in 'solve4'! [t = %f]\n",neutron.t);
      return 0;
    }
  }
  if(cc < 0) {  //It is possible that 'solve3' returns one real root which (due to roundoff) produces a negative value for 'cc'.
    cc=mathzero(cc);
    if(cc != 0) {
      printf("+ERROR: Imaginary 'cc' coefficient for quadratic in 'solve4'! [t = %f]\n",neutron.t);
      return 0;
    }
  }
  
  //Calculate quadratic coefficients for the different possible sign choices.
  bn = a3/2. - sqrt(bc);
  cn = r/2. - sqrt(cc);
  bp = a3/2. + sqrt(bc);
  cp = r/2. + sqrt(cc);
  
   //Determine which combination of signs to use. (Three of the sign-determination equations -- the 'ancillary' equations -- simply force what A&S call p1 and p2
  //to have opposite signs and also force q1 and q2 to have opposite signs.)
  if(mathzero(a3-(bp+bn)) == 0 && mathzero(a2-(bp*bn+cp+cn)) == 0 && mathzero(a0-cp*cn) == 0); //Check the ancillary sign-determination equalities.
  else if(CHECK == "ON") printf(">>Ancillary sign-determination equalities in 'solve4' not satisfied!\n");
  
  if(-4.*a2 + pow(a3,2) + 4.*r >= 0 && -4*a0 + pow(r,2) >= 0) {  //Make certain that radicals in sign-determination equation are defined.
    sdp = -2.*a1 + r*a3 + sqrt(-4.*a2 + pow(a3,2) + 4.*r)*sqrt(-4.*a0 + pow(r,2)); //The p1*q2 + p2*q1 = a1 equation reduces to this for p+/q- or p-/q+.
    sdn = -2.*a1 + r*a3 - sqrt(-4.*a2 + pow(a3,2) + 4.*r)*sqrt(-4.*a0 + pow(r,2)); //The p1*q2 + p2*q1 = a1 equation reduces to this for p-/q- or p+/q+.
    if(CHECK == "ON") printf(">>Sign-Determination: p,q opposite = %e | p,q same = %e\n",sdp,sdn);
  }
  else if(CONSOLE == "ON" || CHECK == "ON") {  //If sign-determination equations are undefined then we'll want to check both possibilities below.
    printf("+WARNING: Sign determination equations undefined in 'solve4'! [%d]\n",neutron.num);
    sdp = -999999;
    sdn = -999999;
  }
  
  //Calculte the discriminants for both possibilities.
  D1same = pow(bp,2) - 4.*cp;
  D2same = pow(bn,2) - 4.*cn;
  D1opp = pow(bp,2) - 4.*cn;
  D2opp = pow(bn,2) - 4.*cp;
  
  if(mathzero(sdn) == 0) {  //Check to see if the value is consistent with zero.
    D1 = D1same;
    D2 = D2same;
    flag1 = -1;
  }
  if(mathzero(sdp) == 0) {  //Check to see if the value is consistent with zero and the discriminants are positive.
    D1 = D1opp;
    D2 = D2opp;
    flag1 = -1;
  }
  
  //Calculate the real roots.
  if(flag1 == -1) {  //A sign-determination was made so we may use that set of discriminants.
    if(D1 >= 0) {
      t1 = (-bp + sqrt(D1))/2.;
      t2 = (-bp - sqrt(D1))/2.;
    }
    if(D2 >= 0) {
      t3 = (-bn + sqrt(D2))/2.;
      t4 = (-bn - sqrt(D2))/2.;
    }
    if(CHECK == "ON") printf(">>Real Roots from solve4: t1 = %e | t2 = %e | t3 = %e | t4 = %e\n",t1,t2,t3,t4);
    
    //Find the smallest nonzero root.
    t=0; //It is physically possible for there to be no positive-time intersection so we don't want to write an error message, but we do want to return exactly zero,
        //so by setting t=0 here if the tests below fail to assign a root to t then zero will be returned.
    if(timezero(t1) > 0) t = t1;
    if(timezero(t2) > 0) t = t2;
    if(timezero(t3) > 0) t = t3;
    if(timezero(t4) > 0) t = t4;
    if(t1 < t && timezero(t1) > 0) t = t1;
    if(t2 < t && timezero(t2) > 0) t = t2;
    if(t3 < t && timezero(t3) > 0) t = t3;
    if(t4 < t && timezero(t4) > 0) t = t4;
    
    flag2 = newton(a3,a2,a1,a0,t,&t,&error); //Apply Newton's method to improve the root.
    if(flag2 == 0 && fabs(error) < 1e-3 && t >= 0) {
      previntersectime = t; //Save this intersection time in case iteration from this root is required in the next call.
      return t; //Return the value of the found and improved root;
    }
    if(flag2 == 0 && fabs(error) < 1e-3 && t < 0) return 0; //We return 0 if there are no physical (nonzero) solutions.
    if(flag2 == -1) flag1 = 0; //If Newton's method produced an error so we should expand our search for a root.
  }
  
  if(flag1 != -1) {  //No definite determination of which signs to use or root finding failed even with a valid sign determination.
    if(CHECK == "ON") printf(">>No clear sign determination in 'solve4'... checking both possibilities:\n");
    //Calculate sets of roots for each possibility and iterates them through Newton's Method. Also iterates the previous root through Newton's Method.
    if(D1same >= 0) {
      if(CHECK == "ON") printf(">>t1same:\n");
      newton(a3,a2,a1,a0,(-bp + sqrt(D1same))/2.,&t1same,&error);
      if(CHECK == "ON") printf(">>t2same:\n");
      newton(a3,a2,a1,a0,(-bp - sqrt(D1same))/2.,&t2same,&error);
    }
    if(D2same >= 0) {
      if(CHECK == "ON") printf(">>t3same:\n");
      newton(a3,a2,a1,a0,(-bn + sqrt(D2same))/2.,&t3same,&error);
      if(CHECK == "ON") printf(">>t4same:\n");
      newton(a3,a2,a1,a0,(-bn - sqrt(D2same))/2.,&t4same,&error);
    }
    if(D1opp >= 0) {
      if(CHECK == "ON") printf(">>t1opp:\n");
      newton(a3,a2,a1,a0,(-bp + sqrt(D1opp))/2.,&t1opp,&error);
      if(CHECK == "ON") printf(">>t1same:\n");
      newton(a3,a2,a1,a0,(-bp - sqrt(D1opp))/2.,&t2opp,&error);
    }
    if(D2opp >= 0) {
      if(CHECK == "ON") printf(">>t1same:\n");
      newton(a3,a2,a1,a0,(-bn + sqrt(D2opp))/2.,&t3opp,&error);
      if(CHECK == "ON") printf(">>t1same:\n");
      newton(a3,a2,a1,a0,(-bn - sqrt(D2opp))/2.,&t4opp,&error);
    }
	if(CHECK == "ON") printf(">>Previous intersection time:\n");
    newton(a3,a2,a1,a0,previntersectime,&told,&error);
    
    //Find the smallest nonzero root.
    t = 99999;
    if(t1same < t && timezero(t1same) > 0) t = t1same;
    if(t2same < t && timezero(t2same) > 0) t = t2same;
    if(t3same < t && timezero(t3same) > 0) t = t3same;
    if(t4same < t && timezero(t4same) > 0) t = t4same;
    if(t1opp < t && timezero(t1opp) > 0) t = t1opp;
    if(t2opp < t && timezero(t2opp) > 0) t = t2opp;
    if(t3opp < t && timezero(t3opp) > 0) t = t3opp;
    if(t4opp < t && timezero(t4opp) > 0) t = t4opp;
    if(told < t && timezero(told) > 0) t = told;
    
    if(t == 99999) return 0;
    error = pow(t,4)+a3*pow(t,3)+a2*pow(t,2)+a1*t+a0; //The value of the quartic at the final value of the approximate root.
    if(fabs(error) < 1e-3 && t >= 0) {
      previntersectime = t;
      return t;
    }
  }

  if(CONSOLE == "ON" || CHECK == "ON") printf("+WARNING: No root found by 'solve4'. [%d]\n",neutron.num);
  return 0;
}

int newton(double a3, double a2, double a1, double a0, double ri, double *rf, double *error) {
  int i,flag;
  double t,thold;
  
  flag = 0;
  t = ri;
  *rf = ri;
  *error = pow(t,4)+a3*pow(t,3)+a2*pow(t,2)+a1*t+a0; //The value of the quartic at the value of the new approximate root.
  if(CHECK == "ON") printf(">>Iterating to improve root...\n  Starting: t = %e | error = %e\n",t,*error);
  
  if(mathzero(*error) == 0) {  //Iteration is not required.
    if(CHECK == "ON") printf("  Ending (after no iterations): t = %e | error = %e.\n",*rf,*error);
    return 0;
  }
  
  for(i=0 ; i < 100 ; i++) {  //Iterate Newton's Rule at most 100 times in order to improve accuracy of the calculated root.
    if((4.*pow(t,3)+3.*a3*pow(t,2)+2.*a2*t+a1) == 0) {  //Derivative in Newton's Rule is zero leading to an undefined calculation.
        if(CHECK == "ON") printf(">>Zero derivative in Newton's Rule iteration after %d steps... iteration terminated.\n",i);
        flag = -1;
        break;
    }
    thold = t - (pow(t,4)+a3*pow(t,3)+a2*pow(t,2)+a1*t+a0)/(4.*pow(t,3)+3.*a3*pow(t,2)+2.*a2*t+a1); //Apply Newton's Rule to improve value of root.
    *error = pow(thold,4)+a3*pow(thold,3)+a2*pow(thold,2)+a1*thold+a0; //The value of the quartic at the value of the new approximate root.
    if(t == thold) break;  //Stop iterating if solution has stopped changing.
    t = thold; //Solution has not stopped changing so update t and keep iterating.
  }
  
  if(flag == -1) return -1;
  if(CHECK == "ON") printf("  Ending (after %d iterations): t = %e | error = %e.\n",i+1,thold,*error);
  *rf = thold;
  return 0;
}

void solve3(double c3, double c2, double c1, double c0,double *pr1,double *pr2,double *pr3) {
  
  struct complex r1,r2,r3,s1,s2,as1,as2,iv;
  double theta,a,b,c,q,r,D;
  
//See Abramowitz and Stegun p.17 for the cubic formula.
  if(c3 == 0) {  //'solve3' called for an equation without a third-order term!
    printf("+ERROR: 'solve3' called for an equation without a third-order term!\n");
    *pr1=0;
    *pr2=0;
    *pr3=0;
    return;
  }
  a = c2/c3;
  b = c1/c3;
  c = c0/c3;
  
  q = b/3. - pow(a,2)/9.;
  r = (b*a - 3.*c)/6. - pow(a,3)/27.;
  D = pow(q,3) + pow(r,2);

  if(D >= 0) {
    iv.r = sqrt(D);
    iv.i = 0;
  }
  if(D < 0) {
    iv.r = 0;
    iv.i = sqrt(-D);
  }

  as1.r = r + iv.r;
  as1.i = iv.i;
  as2.r = r - iv.r;
  as2.i = -iv.i;

  if(D > 0) { //There will be one real root which will come from the z1 equation. 'as1' and 'as2' will be real so the usual cube root should be used.
    if(mathzero(as1.i) != 0 || mathzero(as2.i) != 0) {
      if(CONSOLE == "ON") printf("+WARNING: Nonzero imaginary part for D>0 case in solve3!\n");
    }
    if(as1.r >= 0) s1.r = pow(as1.r,1./3.);
    if(as1.r < 0) s1.r = -pow(-as1.r,1./3.);
    if(as2.r >= 0) s2.r = pow(as2.r,1./3.);
    if(as2.r < 0) s2.r = -pow(-as2.r,1./3.);
    
    r1.r = s1.r + s2.r - a/3.;
    
    *pr1 = r1.r;
    *pr2 = -999999;
    *pr3 = -999999;
  }
  
  if(D == 0) {  //All roots will be real. 'as1' and 'as2' will be real so the usual cube root should be used.
    if(mathzero(as1.i) != 0 || mathzero(as2.i) != 0) {
      if(CONSOLE == "ON") printf("+WARNING: Nonzero imaginary part for D=0 case in solve3!\n");
    }
    if(as1.r >= 0) s1.r = pow(as1.r,1./3.);
    if(as1.r < 0) s1.r = -pow(-as1.r,1./3.);
    if(as2.r >= 0) s2.r = pow(as2.r,1./3.);
    if(as2.r < 0) s2.r = -pow(-as2.r,1./3.);
    
    r1.r = s1.r + s2.r - a/3;
    r1.i = s1.i + s2.i;
    r2.r = (-1./2.)*(s1.r+s2.r) - a/3. - (sqrt(3.)/2.)*(s1.i - s2.i);
    r2.i = (-1./2.)*(s1.i+s2.i) + (sqrt(3.)/2.)*(s1.r - s2.r);
    r3.r = (-1./2.)*(s1.r+s2.r) - a/3. + (sqrt(3.)/2.)*(s1.i - s2.i);
    r3.i = (-1./2.)*(s1.i+s2.i) - (sqrt(3.)/2.)*(s1.r - s2.r);
    if(mathzero(r1.i) != 0 || mathzero(r2.i) != 0 || mathzero(r3.i) != 0) {
      //if(CONSOLE == "ON") printf("+WARNING: A root has a nonzero imaginary part for D=0 case in solve3!\n");
    }
    
   *pr1 = r1.r;
   *pr2 = r2.r;
   *pr3 = r3.r;
 }
 
  if(D < 0) {  //All roots will be real. 'as1' and 'as2' will be imaginary so that the complex cube root should be used.
    s1 = ccuberoot(as1);
    s2 = ccuberoot(as2);
    r1.r = s1.r + s2.r - a/3.;
    r1.i = s1.i + s2.i;
    r2.r = (-1./2.)*(s1.r+s2.r) - a/3. - (sqrt(3.)/2.)*(s1.i - s2.i);
    r2.i = (-1./2.)*(s1.i+s2.i) + (sqrt(3.)/2.)*(s1.r - s2.r);
    r3.r = (-1./2.)*(s1.r+s2.r) - a/3. + (sqrt(3.)/2.)*(s1.i - s2.i);
    r3.i = (-1./2.)*(s1.i+s2.i) - (sqrt(3.)/2.)*(s1.r - s2.r);
    if(mathzero(r1.i) != 0 || mathzero(r2.i) != 0 || mathzero(r3.i) != 0) {
      if(CONSOLE == "ON") printf("+WARNING: A root has a nonzero imaginary part for D<0 case in solve3!\n");
    }
    
    *pr1 = r1.r;
    *pr2 = r2.r;
    *pr3 = r3.r;
  }
  return;
}

struct complex ccuberoot(struct complex z) {
  struct complex cr;
  double theta;
  
  if(z.r == 0) theta = PI/2.; //Prevents an undefined condition in 'atan' function.
  else theta = fabs(atan(z.i/z.r)); //The phase of the complex number as the angle to closest part of the x-axis. Note that 'fabs' is the absolute value
                                   //function that must be used with doubles.
  if(z.r < 0 && z.i >= 0) theta = PI - theta;  //Phase adjusted to be in (-pi,pi] measured from positive real-axis due to branch cut
                                              //in the cube-root function along negative imaginary axis.
  if(z.r <= 0 && z.i < 0) theta = -(PI - theta);
  if(z.r > 0 && z.i < 0) theta = -theta;
  cr.r = pow((pow(z.r,2)+pow(z.i,2)),1./6.)*cos(theta/3.); //Cube root taken in polar form and converted back to
  cr.i = pow((pow(z.r,2)+pow(z.i,2)),1./6.)*sin(theta/3.);//rectangular form.
  return cr; //Return structure holding the result.
}

void poof(int startregion, double poofplaneD, int vdistcode, int adistcode) {
  int i;
  int vdistflag,adistflag; //Set to -1 when a code is recognized.
  int regindex,mcflag;
  double vcutoff; //Sets the v^2 dv speed distribution cut-off.
  double cosrangle1,sinrangle1,mradius,rradius,rangle1,rangle2;
  double dx,dy,dz,dxr,dyr,dzr,thetar,phir,psir;
  double v,vx,vy,vz,vxr,vyr,vzr,test,PofV;
    
  
  vdistflag = 0;
  adistflag = 0;
  
  vcutoff = VCUTOFF; //Specify the cut-off speed for a v^2 dv distribution.
  
//Set Particle's Starting Position
    neutron.region = startregion; //Set starting region.
    thetar = regions[neutron.region][4]; //Get orientation information for the start region.
    phir = regions[neutron.region][5];
    psir = regions[neutron.region][6];
    
    if(regions[neutron.region][0] == 1) {  //The starting region is a box region.
      dx = (grn() - 0.5) * regions[neutron.region][1] - .000001;  //Generate a random displacement in the region's starting plane from basepoint
      dy = (grn() - 0.5) * regions[neutron.region][2] - .000001; //which remains just inside the plane.
      dz = poofplaneD;
      rota(psir,thetar,phir,dx,dy,dz,&dxr,&dyr,&dzr); //Rotate random displacement from basepoint by region's orientation.
    }
    if(regions[neutron.region][0] == 2) {  //The starting region is a cylinder region.
      rangle1 = 2.*PI*grn(); //Generate a random angle [0,2 PI].
      mradius = regions[neutron.region][1]/2. - .000001; //Get a maximum radius which is just inside the region.
      rradius = mradius * sqrt(grn()); //Generate a random radius with the correct r-weighting via the Inverse Transform method.
      dx = rradius * cos(rangle1); //Calculate random displacement from basepoint defined by random values just calculated.
      dy = rradius * sin(rangle1);
      dz = poofplaneD;
      rota(psir,thetar,phir,dx,dy,dz,&dxr,&dyr,&dzr); //Rotate displacement vector into the region's orientaion.
    }
  
  neutron.x = basepoints[neutron.region][0] + dxr; //Set the particle's initial position.
  neutron.y = basepoints[neutron.region][1] + dyr;
  neutron.z = basepoints[neutron.region][2] + dzr;

//Set the particle's starting velocity.
  //Set the starting direction. NB: We are assuming that the region orientation information loaded in getting starting point is still correct.
  if(adistcode == 0) {  //Here we only want velocities heading into the starting region, i.e. into the half-space.
    adistflag = -1; //Code recognized.
    cosrangle1 = sqrt(1. - grn()); //Generate cos(theta) for theta in [0,PI/2] consistent with an isotropic distribution illuminating a plane, i.e. proportional
                                  //to cos(theta) d(Omega) via the Inverse Transform method.
    sinrangle1 = sqrt(1 - pow(cosrangle1,2)); //Calculate the corresponding sin(theta) from the generated cos(theta).
    rangle2 = 2.*PI*grn(); //Random angle in [0,2 PI] (phi).
  }
  if(adistcode == 1) {  //Here we want a completely isotropic distribution, i.e. into the full 4pi.
    adistflag = -1; //Code recognized.
    cosrangle1 = 1. - 2. * grn(); //Generate cos(theta) for theta in [0,PI] consistent with an isotropic distribution, i.e. proportional to d(Omega) only
                                 //via the Inverse Transform method.
    sinrangle1 = sqrt(1 - pow(cosrangle1,2)); //Calculate the corresponding sin(theta) from the generated cos(theta).
    rangle2 = 2.*PI*grn(); //Random angle in [0,2 PI] (phi).
  }
    
  //Set the starting speed.
  if(vdistcode == 0) {  //We want a v^2 dv distribution.
    vdistflag = -1; //Code recognized.
    v = vcutoff*pow(grn(),1./3.);
  }
  if(vdistcode == 1) {  //We want to generate a speed distribution via the Monte Carlo method.
    vdistflag = -1; //Code recognized.
    mcflag = 0;
    while (mcflag == 0) {
      v = vcutoff*grn();
      test = grn();
      PofV = 3./pow(vcutoff,3) * pow(v,2);
      if (test < PofV) mcflag = 1;
    }
  }
  if (vdistcode == 2) {  //We want a monoenergetic energy distribution
    vdistflag = -1; //code recognized
    v = MONOENERGY;
  }
    
  //Find the pre-rotated components of the velocity.
  vx = v*sinrangle1*cos(rangle2);
  vy = v*sinrangle1*sin(rangle2);
  vz = v*cosrangle1;
  rota(psir,thetar,phir,vx,vy,vz,&vxr,&vyr,&vzr); //Rotate the velocity vector into the region's orientation.
  
  neutron.vx = vxr; //Set the particle's initial velocity.
  neutron.vy = vyr;
  neutron.vz = vzr;

  //if(EVENTS == "ON") event(14,cosrangle1,rangle2,v,0); //Write an event for the creation of a particle.
  event(14,cosrangle1,rangle2,.5*nMASS*(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2))/neV2J,.5*nMASS*pow(vz,2)/neV2J); //Write an event for the creation of a particle.
  
  
  if(vdistflag != -1 || adistflag != -1) printf("+ERROR: Code for particle velocity distribution not recognized in 'poof'!\n");
  return;
}

double grn(void) {
  return rand()/(double)RAND_MAX; //Generates a psuedo-random number in [0,1] using the seed defined in header. (The 'rand' function returns a
                                 //pseudo-random number between 0 and RAND_MAX, which is presumably defined in the 'stdlib.h' file which defines 'rand').
}

double timezero(double val) {
  if(val < TZERO && val > -TZERO) return 0;
  return val;
}

double mathzero(double val) {
  if(val < MZERO && val > -MZERO) return 0;
  return val;
}

int propagate(void) {
  int i,conreg,tsct,propmode;
  int eventcode, intrsct, blkmed;
  double t, dt;
  double pdymhold[7]; //Array to hold particle's initial t,x,y,z,vx,vy,vz in case trajectory tracking in ON.
  double traj[3][3];

//Particle must be transported through a non-trivial magnetic field.  
  if((int)rparams[neutron.region][11] == -1) {}
  
//Integration through a magnetic field not required.
  else {
    trajmatrix(traj); //Assemble trajectory matrix (relative to global coordinate system) for a particle moving with accelerations given in 'regions'.
    dt = intersection(&eventcode,traj,&intrsct); //Calculate the time to the nearest intersection with the gemetry.
    if(eventcode == -1) {  //An error occurred in calculating an intersection.
      event(-1,0,0,0,0);
      return -1;
    }
    
    neutron.xcode = intrsct; //Set neutron 'xcode' value to indicate any cut-plane intersections (or T-juntion intersections).
                            //NB: For T-juntion intersections 'xcode' holds the region number the particle is entering.
    
    //Check for absorption and/or scattering.
    propmode = (int)rparams[neutron.region][11];
    if(propmode == 1 || propmode == 2 || propmode == 3 || propmode == 4) {  //Region possesses a medium
      blkmed = propinmed(dt,traj); //Check for interactions with the bulk medium.
      if(blkmed == 8) {  //Particle absorbed in a non-detector region.
        if(EVENTS == "ON") event(8,0,0,0,0);
        return 8;
      }
      if(blkmed == 9) {  //Particle was absorbed in a detector region.
		if(EVENTS == "ON") event(9,0,0,0,0);
        return 9;
      }
      if(blkmed == 17) {  //Particle was scattered in the bulk medium.
        if(EVENTS == "ON") event(17,0,0,0,0);
        return 17;
      }
      if(blkmed == -1) {  //An error occurred in the medium interaction calculations.
        event(-1,0,0,0,0);
        return -1;
      }
    }
            
    //Check for Beta-Decay.
    if(grn() < 1./888. * dt) {
      if(EVENTS == "ON") event(11,0,0,0,0);
      return 11;
    }
    
    //Create trajectory events if trajectory tracking is turned ON.
    if(TRAJFLAG == "ON" || propmode == 5) {
      //Save initial dynmical information so that after trajectory tracking we can be sure to move particle properly to intersection.
      pdymhold[0] = neutron.t;
      pdymhold[1] = neutron.x;
      pdymhold[2] = neutron.y;
      pdymhold[3] = neutron.z;
      pdymhold[4] = neutron.vx;
      pdymhold[5] = neutron.vy;
      pdymhold[6] = neutron.vz;
      for(t = TRAJTS ; t <= dt ; t = t+TRAJTS) {  //Generate trajectory points over dt at intervals defined by TRAJTS in the header.
        if(move(t,traj,0) == -1) break; //Move particle along its trajectory by time TRAJTS. Note that while 'move' does change the particle's dynamical information,
                                       //it does not change the 'traj' array so that the particle's initial position in the trajectory equations is unchanged during
                                      //this loop. Thus we must 'move' through a t which increases with TRAJTS, not just by TRAJTS. If the particle is outside of
                                     //the geometry trajectory recording will cease.
        neutron.t = neutron.t - t + TRAJTS; //'move' increments the neutron time by the first argument, but in this case the first argument isn't a
                                           //time step from last point, so we must correct the neutron time.
        if(EVENTS == "ON") event(0,0,0,0,0); //Record a trajectory event.
      }
      //Reset particle its original state.
      neutron.t = pdymhold[0];
      neutron.x = pdymhold[1];
      neutron.y = pdymhold[2];
      neutron.z = pdymhold[3];
      neutron.vx = pdymhold[4];
      neutron.vy = pdymhold[5];
      neutron.vz = pdymhold[6];
    }

    if(eventcode != 1 && eventcode != 2 && eventcode != 16) {  //We should only be moving the particle to a cut-plane, T-junction, or geometry intersection.
      printf("+ERROR: Function 'propagate' is trying to move particle %d to a location associated with a non-intersection event!\n",neutron.num);
      return -1;
    }
    
    if(eventcode == 1 || eventcode == 16) if(move(dt,traj,0) == -1) return -1; //The particle has intersected a T-junction or an element of the geometry. Move the
                                                                              //particle to the new intersection, update time, position, and velocity.
                                                                              //Drop the particle if it is irredeemably outside of the geometry.
    if(eventcode == 2) if(move(dt,traj,1) == -1) return -1; //The particle has intersected a cut-plane. Move the particle to the new intersection, update time, 
                                                           //position, and velocity. Drop the particle if it is irredeemably outside of the geometry or off the
                                                          //cut-plane.

    //Now that the particle has been moved to its intersection, we must check to see if its current position puts it at a T-juntion.
    if((int)cplanes[neutron.region][4] == 2) {  //The particle intersected a region that is T'ed into so we must check to see if
                                               //the intersection is at a junction.
      for(i=0 ; i < 6 ; i++) {  //Scan through connex file for current region.
        tsct = 0; //Reset the T-junction intersection flag.
        conreg = connex[neutron.region][i]; //Connected region.
        if(conreg == -1) break; //We have reached the end of the list of connected regions.
        if((int)cplanes[conreg][4] == 3 || (int)cplanes[conreg][4] == 4) {  //The connected region T's into the current region.
          tsct = tjunc(conreg); //Determine whether the intersection is at this T-junction.
          if(tsct == -1) {  //An error occurred.
            event(-1,0,0,0,0);
            return -1;
          }
          if(tsct == 1) {  //We have intersected a T-junction.
            neutron.xcode = conreg; //Set the code indicating a transition between two regions.
            if(CHECK == "ON") {
              printf(">>Intersection with current (main) region #%d is at a T-junction: (%f,%f,%f)\n",neutron.region,neutron.x,neutron.y,neutron.z);
              printf("T-ing Region %d's cut-plane is at: %f,%f,%f",conreg,basepoints[conreg][0],basepoints[conreg][1],basepoints[conreg][2]);
              printf(" (t-junc. flag = %d) | Length of T-ing region: %f\n",(int)cplanes[conreg][6],regions[conreg][2]);
            }
            eventcode = 16; //We also must set the eventcode for a T-junction intersection so that we end up in 'cplanehandling'.
            break;
          }
        }
      }
    }
    if(EVENTS == "ON") event(eventcode,0,0,0,0); //Record the intersection event.
  }
  return eventcode;
}

int bounce(int parareg, double mpot, int dirflag, double passnx, double passny, double passnz) {
  int depolflag,btype,losschk;
  double vxb,vyb,vzb,vmag;
  double pang,azang;

  if(dirflag != 1 && dirflag != -1) {
    printf("+ERROR: 'bounce' called with undefined direction flag! [Particle#: %d]\n",neutron.num);
    return -1;
  }
  
  gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vxb,&vyb,&vzb,passnx,passny,passnz); //Transform the particle's velocity components into the bounce system.


  //Calculate the angle that the velocity vector makes with the appropriate normal to the surface (which in the bounce system is defined to be the x-axis).
  if(fabs(vxb) == 0) pang = PI/2.;
  else pang = atan(sqrt(pow(vyb,2)+pow(vzb,2))/fabs(vxb));
  azang = atan2(vyb,vzb); //Calculate the azimuthal angle of the velocity vector relative to the longitudinal axis (z-axis in the bounce system).
  if(EVENTS == "ON") event(12,pang,azang,vxb,vzb); //Write an event for entering a bounce.
  
  losschk = loss(parareg,mpot,passnx,passny,passnz); //Check for wall loss.
  if(losschk == -1) return -1; //Error in determining wall loss.
  if(losschk == 1) return 1; //Particle penetrates wall.
  if(losschk == 2) return 2; //Particle lost to wall through physical mechanism.
  btype = bouncetype(parareg,pang); //Determine whether the bounce is specular or non-specular.
  if(btype == -1) return -1; //Error in determining the type of bounce.
  depolflag = depol(parareg); //Check to see if the particle depolarized.
  if(depolflag == -5) return -1; //Error in determining the polarization state.

  if(btype == 0) {  //The bounce is specular.
    if((vxb >= 0 && dirflag == 1) || (vxb <= 0 && dirflag == -1)) {  //Particle has no velocity component into the wall.
      //printf("Vxb in bounce: %f\n",vxb);
      if(CONSOLE == "ON") printf("+WARNING: 'bounce' called for particle with no velocity component into wall! [Part#: %d | Reg#: %d | Pos: (%f,%f,%f) | Vx: %f | dirflag: %d]\n",neutron.num,neutron.region,neutron.x,neutron.y,neutron.z,vxb,dirflag);
	  event(-1,0,dirflag,vxb,0);
      return -1;
    }
    
    
    
    //printf("Vxb in bounce: %f, dirflag: %d, Reg: %d, xcode: %d\n",vxb,dirflag,neutron.region, neutron.xcode);
    vxb = -vxb;                                                  //Reverse the velocity vector's normal component (the x-component in the bounce system) for a specular bounce.
    bsys2gsys(vxb,vyb,vzb,&neutron.vx,&neutron.vy,&neutron.vz,passnx,passny,passnz); //Transform the "bounced" velocity vector back to the global system
                                                               //and update the particle's velocity vector.
      
                                                               //Calculate the exit polar and azimuthal angles and write the appropriate exit event.
    if(fabs(vxb) == 0) pang = PI/2.;                          //Calculate the angle that the velocity vector makes with the appropriate normal to the surface,
    else pang = atan(sqrt(pow(vyb,2)+pow(vzb,2))/fabs(vxb)); //which in the bounce system is defined to be the x-axis.
    azang = atan2(vyb,vzb);                                 //Calculate the azimuthal angle of the velocity vector relative to the longitudinal axis (z-axis in the bounce system).
                                                           //Write appropriate event depending on polarization state.
    if(depolflag == -1) {  //Particle was depolarized.
      if(EVENTS == "ON") event(5,pang,azang,vxb,vzb);
      return 0;
    }
    if(depolflag == 0) {  //Particle was not depolarized.
      if(EVENTS == "ON") event(3,pang,azang,vxb,vzb);
      return 0;
    }
  }
    
  if(btype == 1) {  //The bounce in non-specular.
    if(nsbounce(parareg,dirflag,depolflag,vxb,vyb,vzb) == 0) return 0; //Generate the appropriate particle velocity for a non-specular bounce.
    else {
      printf("+ERROR: Undefined scattering model in 'bounce'!\n");
      event(-1,0,0,0,0);
      return -1; //Error event code returned.
    }
  }
    
  printf("+ERROR: 'bounce' called with no result. Particle number %d abandoned!\n",neutron.num);
  return -1;
}

int nsbounce(int parareg, int dirflag, int depolflag, double vxb, double vyb, double vzb) {
  int check;
  double cosexpang,sinexpang,exazang,pang,azang,vmag;
  
  check = 0;
  if(rparams[parareg][12] == 0) {  //Diffuse scattering on a non-specular bounce.
    if(dirflag == 1) cosexpang = sqrt(1. - grn()); //Generate cos(polar exit angle) for exit angle in [0,Pi/2] consistent with an isotropic exit distribution via the
                                                  //Inverse Transform method. Note that this will be the cos of the angle from the +x-axis in the bounce system.
    if(dirflag == -1) cosexpang = -1.*sqrt(grn()); //If the normal vector is outward (as in the case of a cut-plane intersection where the particle comes from
                                                  //a region not associated with the cut-plane) then we want our polar angle to be in [Pi/2,Pi].
    sinexpang = sqrt(1-pow(cosexpang,2)); //Calculate the sin of the exit angle from the generated value of cos.
    exazang = 2.*PI*grn(); //Generate a random azimuthal angle in [0,2 Pi]. (This will be relative to the +z-axis in the bounce system.)
    check = 1; //Set the flag indicating that exit information has been calculated.
  }
  
  if(rparams[parareg][12] == 1) {  //cos-squared weighted exit angle on a non-specular bounce.
    if(dirflag == 1) cosexpang = pow((1. - grn()),1./3.); //Generate cos(polar exit angle) for exit angle in [0,Pi/2] consistent with a sin(theta)*cos(theta)
                                                         //exit distribution, i.e. isotropic cos(theta)-weighted distribution via the Inverse Transform method.
                                                        //Note that this will be the cos of the angle from the +x-axis in the bounce system.
    if(dirflag == -1) cosexpang = -1.*pow(grn(),1./3.); //If the normal vector is outward (as in the case of a cut-plane intersection where the particle comes from
                                                       //a region not associated with the cut-plane) then we want our polar angle to be in [Pi/2,Pi].
    sinexpang = sqrt(1-pow(cosexpang,2)); //Calculate the sin of the exit angle from the generated value of cos.
    exazang = 2.*PI*grn(); //Generate a random azimuthal angle in [0,2 Pi]. (This will be relative to the +z-axis in the bounce system.)
    check = 1; //Set the flag indicating that exit information has been calculated.
  }
  
  if(check == 1) {  //Exit angles were found... set exit velocity appropriately.
    vmag = sqrt(pow(vxb,2)+pow(vyb,2)+pow(vzb,2)); //Calculate the magnitude of the particle's velocity vector.
    vxb = vmag*cosexpang;               //Calculate the components of the exit velocity. Note that since the x-axis is our polar axis and the azimuthal angle is
    vyb = vmag*sin(exazang)*sinexpang; //measure relative to the +z-axis, the usual formulae have been modified.
    vzb = vmag*cos(exazang)*sinexpang;
      
    bsys2gsys(vxb,vyb,vzb,&neutron.vx,&neutron.vy,&neutron.vz,0,0,0); //Transform the "bounced" velocity vector back to the global system
                                                               //and update the particle's velocity vector.

    if(fabs(vxb) == 0) pang = PI/2.;                          //Calculate the angle that the velocity vector makes with the inward normal to the surface,
    else pang = atan(sqrt(pow(vyb,2)+pow(vzb,2))/fabs(vxb)); //which in the bounce system is defined to be the x-axis.
    azang = atan2(vyb,vzb); //Calculate the azimuthal angle of the velocity vector relative to the longitudinal axis (z-axis in the bounce system).
    if(depolflag == -1) {  //Particle was depolarized.
      if(EVENTS == "ON") event(6,pang,azang,vxb,vzb); //Write non-specular, depolarized bounce event.
      return 0;
    }
    if(depolflag == 0) {  //Particle was not depolarized.
      if(EVENTS == "ON") event(4,pang,azang,vxb,vzb); //Write non-specular, not depolarized bounce event.
      return 0;
    }
  }
  
  if(CONSOLE == "ON") printf("+WARNING: No non-specular scattering result for particle %d at time %f in 'nsbounce'!\n",neutron.num,neutron.t);
  return -1; //If we make it to the end without a result then an error occurred.
}

int bouncetype(int parareg, double nang) {
  
  if(rparams[parareg][12] == 0) {  //Fixed probability of a specular bounce.
    if(grn() < rparams[parareg][3]) return 0; //The bounce is specular.
    else return 1; //The bounce is non-specular.
  }
  
  if(rparams[parareg][12] == 1) {  //Probability of a non-specular bounce is (1-spec) with cos-weighting.
    if(grn() < (1. - rparams[parareg][3])*cos(nang)) return 1; //The bounce is non-specular.
    else return 0; //The bounce is specular.
  }
  
  if(CONSOLE == "ON") printf("+WARNING: No specular/non-specular determination for particle %d at time %f in 'bouncetype'!\n",neutron.num,neutron.t);
  return -1;
}
  
int loss(int parareg, double mpot, double passnx, double passny, double passnz) {
  double nx,ny,nz;
  double vperp,Eperp,TransCoeff,ReflCoeff;
  double rpotp2,ipotp2,pperp;
  double mth1,mth2,mth3,numtr,denmtr,lprob;
  int superflag;
  
  //Check to see if the particle's energy "perpendicular to the wall" is above the wall's Fermi potential.
  normal(&nx,&ny,&nz,passnx,passny,passnz); //Get normal vector to current intersection point.]
  //printf("loss pass = (%f,%f,%f), loss n = (%f,%f,%f).\n", passnx,passny,passnz,nx,ny,nz);
  vperp = neutron.vx*nx + neutron.vy*ny + neutron.vz*nz; //Calculate the component of the particle velocity perpendicular to the wall.
  //printf("vperp = %f, pass = (%f,%f,%f).\n",vperp,passnx,passny,passnz);
  Eperp = 1./2.*nMASS*pow(vperp,2); //Calculate the 'perpendicular' energy.
  superflag=0;

  if(mpot < 0) {
  	ReflCoeff = pow((sqrt(Eperp)-sqrt(Eperp-mpot*neV2J)),2)/pow((sqrt(Eperp)+sqrt(Eperp-mpot*neV2J)),2); //Reflection Coefficient from Golub's calculation.
	if(grn() > ReflCoeff) {
        if(EVENTS == "ON") event(13,Eperp,mpot,0,0);
        return 1;
	}
	else return 0;
	}
  if(Eperp > mpot*neV2J) {
	TransCoeff = 4.*sqrt(Eperp)*sqrt(Eperp-mpot*neV2J)/pow(sqrt(Eperp)+sqrt(Eperp-mpot*neV2J),2);
	if(grn() < TransCoeff){   //Particle penetrates the material and is lost.
    	if(EVENTS == "ON") event(13,Eperp,mpot,0,0);
        return 1;
	} //Particle penetrates wall.
	else superflag=1;
  }


  //Check to see if the particle is lost by other mechanisms.

  if((int)rparams[parareg][13] == 0) {  //Determination of loss using a fixed probability of loss from mechanisms other than wall penetration.
//    if(grn() < rparams[parareg][4]) {  //Random value is less than the probability of loss on a single bounce.
	  if(superflag != 1){
	  if(grn() < (rparams[parareg][4]*2*fabs(vperp)/(sqrt(2*mpot*neV2J/nMASS - pow(vperp,2))))) { //Ignatovich Wall Loss Formula (G. Palmquist)
      if(EVENTS == "ON") event(7,0,0,0,0);
      return 2; //Particle is lost to the wall.
      }
    }
    return 0; //Particle is not lost to the wall.
  }


  if((int)rparams[parareg][13] == 1) {  //Energy and angle dependent wall loss via an imaginary part to the wall potential using Morozov formula.
    pperp = -nMASS*vperp; //The component of particle momentum perpendicular to the wall.
    if(pperp < 0) {
      printf("+ERROR: Particle with normal velocity component not heading into the wall in 'loss' for neutron# %d!\n",neutron.num);
      return -1;
    }
    rpotp2 = 2.*nMASS*mpot*neV2J; //Real part of the wall potential times 2m. (Ur = pr^2/2m)
    ipotp2 = -2.*nMASS*rparams[parareg][4]*neV2J; //Imaginary part of the wall potential times 2m. (Note that the imaginary part of the wall potential
                                                        //will be negative so that Ui = -pi^2/2m.)
    mth1 = pow(pperp,2) - rpotp2;
    mth2 = pow(mth1,2) + pow(ipotp2,2);
    if (mth2 < 0) {
      printf("+ERROR: Imaginary Reflection Probability in 'loss' for neutron# %d!\n",neutron.num);
      return -1;
    }
    mth3 = sqrt(mth2);
    if(mth1+mth3 < 0) {
      printf("+ERROR: Imaginary Reflection Probability in 'loss' for neutron# %d!\n",neutron.num);
      return -1;
    }
    numtr = 2.*sqrt(2)*pperp*sqrt(mth1+mth3);
    denmtr = pow(pperp,2) + mth3 + numtr/2.;
    lprob = numtr/denmtr; //Probability that the particle is lost (i.e. NOT reflected) due to its interaction with the wall.
    if(grn() < lprob) {  //Random value is less than the calculated loss probability so the particle is lost.
      if(EVENTS == "ON") event(7,0,0,0,0);
      return 2;
    }
    return 0; //Particle was not lost to the wall.
  }
  
  printf("+ERROR: Undefined Loss Model!\n");
  return -1;
}

int depol(parareg) {

  if((int)rparams[parareg][14] == 0) {  //Determination of depolarization via a specified probability.
    if(grn() < rparams[parareg][5]) {  //Random number is less than probability of depolarization on a single bounce, i.e. particle depolarizes.
      neutron.spin[2] = -neutron.spin[2]; //Flip the spin of the particle.
      return -1; //Indicate that depolarization has occurred.
    }
    return 0; //Particle did not depolarize.
  }

  if((int)rparams[neutron.region][14] == 1) {} //Future depolarization models go here.
  
  printf("+ERROR: Undefined Depolarization Model!\n");
  return -5;
}

int propinmed(double dt, double traj[3][3]) {
  
  int propmode;
  double t,D,xold,yold,zold;
  double costheta, sintheta, phi;
  double scattfuzz,tstep,v;
  double alpha,transm,scatt;
  double KEmed, KEvac;
  double kmed, kvac;
  double pdymhold[7];
  
  propmode = (int)rparams[neutron.region][11];
  
  //Absorbtion
  alpha = exp(-1.*rparams[neutron.region][9]*dt); //Calculate the probability of transmission through the medium.
  if(propmode == 1 || propmode == 4) {  //Propability of absorption is due to the medium alone.
    if(grn() < 1. - alpha) {  //Check to see if the particle is absorbed.
      if((int)rparams[neutron.region][10] == 0) return 8; //Particle was absorbed in a non-detector region.
      else return 9; //Particle was absorbed in a detector region.
    }
  }
  if(propmode == 2) {  //Probability of absorption is due to the medium and internal reflections (i.e. a thin foil).
    KEmed = 1./2.*nMASS*(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)); //This is the kinetic energy in the medium since if the particle has entered the
                                                                           //medium it has already had its energy shifted.
    KEvac = KEmed + rparams[neutron.region][7]*neV2J; //We assume vacuum on either side of a thin foil.
    if(KEmed < 0 || KEvac < 0) {  //Yipes!
      printf("+ERROR: Negative Kinetic Energy in 'propinmed' for particle number %d!\n",neutron.num);
      return -1;
    }
    kmed = sqrt(2.*nMASS/pow(hBAR,2)*KEmed); //Calculate the particle's wave number in the medium,
    kvac = sqrt(2.*nMASS/pow(hBAR,2)*KEvac); //Calculate the particle's wave number in vacuum.
    transm = 16.*alpha*pow(kmed,2)*pow(kvac,2)/(pow((kmed+kvac),4) - pow(alpha,2)*pow((kvac-kmed),4)); //Calculate the transmission probability.   
    if(transm < 0 || transm > 1) {  //Yipes!
      printf("+ERROR: Transmission probability for particle number %d not in range [0,1] in 'propinmed'!\n",neutron.num);
      return -1;
    }
    if(grn() > transm) {  //Check to see if the particle is transmitted.
      if((int)rparams[neutron.region][10] == 0) return 8; //Particle was absorbed in a non-detector region.
      else return 9; //Particle was absorbed in a detector region.
    }
  }
  
  //Scattering
  if(propmode == 3 || propmode == 4) {  
    scattfuzz = 1.0e-5; //This sets the localization of scattering, i.e. the particle's trajectory is divided up into intervals of (approximately) this
                     //length (in meters) and scattering is checked for on each interval.
    //Save initial dynmical information so that if there is no scattering we can reset the particle's initial dynamical informtion.
    pdymhold[0] = neutron.t;
    pdymhold[1] = neutron.x;
    pdymhold[2] = neutron.y;
    pdymhold[3] = neutron.z;
    pdymhold[4] = neutron.vx;
    pdymhold[5] = neutron.vy;
    pdymhold[6] = neutron.vz;
  
    v = sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)); //Calculate the particle's speed.
    tstep = scattfuzz/v; //The time step which will move the particle along the trajectory a distance equal to the defined interval.
    for(t=tstep ; t <= dt ; t = t+tstep) {  //Move the particle along its trajectory looking for scattering events.
      xold = neutron.x;
      yold = neutron.y;
      zold = neutron.z;
      if(move(t,traj,0) == -1) return -1; //Move the particle to the next time step and abandon it if it leaves the geometry.
      neutron.t = neutron.t - t + tstep; //'move' moves the particle from the point of last intersection so that we must move with incresingly large
                                        //temporal displacements in order to move along the trajectory. However, 'move' also updates the particle time each time
                                       //it is called, so we must remove the change from the last move since nothing happened, and increase the particle time by
                                      //the new temporal displacement.
      D = sqrt(pow(xold-neutron.x,2)+pow(yold-neutron.y,2)+pow(zold-neutron.z,2)); //Calculate the (straight-line) displacement of the particle.
      scatt = 1. - exp(-D/rparams[neutron.region][8]); //Calculate the probability of a scattering event.
      if(grn() < scatt) {  //The particle scattered.
        v = sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)); //Calculate the speed of the paticle at point of scatter. (It can change from the
                                                                        //previously calculated value due to gravity.)
        costheta = 1. - 2.*grn(); //Pick an isotropic theta in [0,Pi];
        sintheta = sqrt(1. - pow(costheta,2)); //Calculate the sin of the isotropic angle.
        phi = 2.*PI*grn(); //Pick a random phi in [0,2 Pi].
        if(EVENTS == "ON") event(17,0,0,0,0); //Write an event to record the particle's dynamical information just before the scatter.
        neutron.vx = v*sintheta*cos(phi);  //Calculate the new velocity vector
        neutron.vy = v*sintheta*sin(phi); //for the randomly chosen exit direction
        neutron.vz = v*costheta;         //assuming an elastic scatter.
        return 17; //Return the event code for a scattering event.
      }
    }
    //The particle did not scatter so reset the particle its original state.
    neutron.t = pdymhold[0];
    neutron.x = pdymhold[1];
    neutron.y = pdymhold[2];
    neutron.z = pdymhold[3];
    neutron.vx = pdymhold[4];
    neutron.vy = pdymhold[5];
    neutron.vz = pdymhold[6];
  }
    
  return -2; //Particle was not absorbed or scattered in the medium.
}

void normal(double *nx, double *ny, double *nz, double passnx, double passny, double passnz) {
  int face;
  double foff1,foff2,foff3,foff4;
  double nxp,nyp,nzp;
  double xp,yp,zp;
  double bpx,bpy,bpz;
  double angle,thetar,phir,psir;
  
  if (passnx != 0 && passny != 0 && passnz != 0) {
    *nx = passnx;
    *ny = passny;
    *nz = passnz;
    return;
  }
  
  *nx = 0;  *ny = 0; *nz = 0; //If there is an error, the normal vector (0,0,0) will be returned.
  if(neutron.xcode != -1 && neutron.xcode < 0) {  //'normal' should not be called unless there is a current intersection.
    printf("+ERROR: Call to 'normal' without a current intersecton!\n");
    return;
  }
  
  if(neutron.xcode == -1) {  //Intersection is with current region.
    //Get basepoint and orientation for the region.
    bpx = basepoints[neutron.region][0];
    bpy = basepoints[neutron.region][1];
    bpz = basepoints[neutron.region][2];
    thetar = regions[neutron.region][4];
    phir = regions[neutron.region][5];
    psir = regions[neutron.region][6];
    
    rotp(psir,thetar,phir,neutron.x-bpx,neutron.y-bpy,neutron.z-bpz,&xp,&yp,&zp); //Transform particle's position into the solving system.
    
    if((int)regions[neutron.region][0] == 1) {  //Region is a box-region.
      //Calculate the particle's offset from each face.
      foff1 = fabs(xp - regions[neutron.region][1]/2.); //+x face offset.
      foff2 = fabs(xp + regions[neutron.region][1]/2.); //-x face offset.
      foff3 = fabs(yp - regions[neutron.region][2]/2.); //+y face offset.
      foff4 = fabs(yp + regions[neutron.region][2]/2.); //-y face offset.
      
      //Determine with which face the particle has intersected. (Remember that +z and -z faces are cut-planes and so are checked elsewhere.)
      face = 0;
      if(foff1 < foff2 && foff1 < foff3 && foff1 < foff4) face = 1; //Face at +x-Dim.
      if(foff2 < foff1 && foff2 < foff3 && foff2 < foff4) face = 2; //Face at -x-Dim.
      if(foff3 < foff1 && foff3 < foff2 && foff3 < foff4) face = 3; //Face at +y-Dim.
      if(foff4 < foff1 && foff4 < foff2 && foff4 < foff3) face = 4; //Face at -y-Dim.
      if(face == 0) printf("+ERROR: Intersection face of box-region not identified in 'normal'!\n");
      
      //Calculate the normal vector in the global system based on the intersection face.
      if(face == 1) rota(psir,thetar,phir,-1.,0.,0.,nx,ny,nz); //Note that nx, ny, and nz are already pointers so no "&" operators are needed.
                                                         //+x-Dim face has a normal in solving system in -x direction.
      if(face == 2) rota(psir,thetar,phir,1.,0.,0.,nx,ny,nz);
      if(face == 3) rota(psir,thetar,phir,0.,-1.,0.,nx,ny,nz);
      if(face == 4) rota(psir,thetar,phir,0.,1.,0.,nx,ny,nz);
      return;
    }
    
    if((int)regions[neutron.region][0] == 2) {  //Region is a cylinder-region.
      //Check that the particle is indeed currently on the surface of the cylinder-region.
      if((pow(xp,2)+pow(yp,2)) > pow(regions[neutron.region][1]/2.,2) + GZERO || (pow(xp,2)+pow(yp,2)) < pow(regions[neutron.region][1]/2.,2) - GZERO) {
        printf("+ERROR: Normal sought for point not on surface of cylinder-region in 'normal'!\n");
        return;
      }
      angle = atan2(yp,xp); //Calculate the angle from the +x-axis that the particle's position makes in the solving system.
                           //Note that 'atan2(a,b)' calculates the inverse tan of a/b and correctly determines the quadrant based on the signs of a and b.
      //Calculate the components of the normal vector in the solving system.
      nxp = -cos(angle);
      nyp = -sin(angle);
      nzp = 0.;
      
      rota(psir,thetar,phir,nxp,nyp,nzp,nx,ny,nz); //Transform the normal vector back into the global coordinate system.
                                                  //Recall that rota when used as a passive transformation is the inverse of rotp.
                                                 //NB: nx, ny, and nz are already defined as pointers so that no "&" operators are required.
      return;
    }
    printf("+ERROR: Undefined region type in 'normal'!\n");
    return;  
  }
  
//  if (tjunc(neutron.region) == 1 && tjunc(neutron.xcode) == -1) { //particle is at a t-junction so we must find the normal vector of the correct cutplane
//    printf("\nRegion: %d. Xcode: %d.\n",neutron.region, neutron.xcode);
//  }
  
  //Intersection is with a cut-plane.
  //Get the basepoint for the region associated with the cut-plane and the orientation of the cutplane.
  bpx = basepoints[neutron.xcode][0];
  bpy = basepoints[neutron.xcode][1];
  bpz = basepoints[neutron.xcode][2];
  thetar = cplanes[neutron.xcode][0];
  phir = cplanes[neutron.xcode][1];
  psir = cplanes[neutron.xcode][7];
    
  rota(psir,thetar,phir,0.,0.,1.,nx,ny,nz); //In the solving system the normal to a cut-plane will always be in the +z-direction (into the region associated
                                           //with the cut plane) so all we must do is transform that vector back into the global coordinate system.
  return;
}

void gsys2bsys(double gcx,double gcy,double gcz,double *bcx,double *bcy,double *bcz, double passnx, double passny, double passnz) {
  double nx,ny,nz,lx,ly,lz,tx,ty,tz,mag;

  if(neutron.xcode != -1 && neutron.xcode < 0) {  //Particle does not have a current intersection.
    printf("+ERROR: No current particle intersection in 'gsys2bsys'!\n");
    *bcx = -999999;  *bcy = -999999;  *bcz = -999999; //If there is an error then we return something which should produce a geometry error.
    return;
  }
  
  normal(&nx,&ny,&nz,passnx,passny,passnz); //Calculate normal at current point of intersection. [Bounce system x-axis.]
  
//  if (passnx != 0 && passny != 0 && passnz != 0) {
//    printf("In gsys2bsys() we have a normal of (%f,%f,%f).\n\n",nx,ny,nz);
//  }
  
  //Generate a unit vector perpendicular to the unit normal. [Bounce system's z-axis.]
  if(timezero(ny) != 0. || timezero(nz) != 0.) {  //We may generate a perpendicular vector by crossing the normal with the global x-hat vector.
    lx = 0.;
    ly = nz;
    lz = -ny;
    mag = sqrt(pow(ly,2)+pow(lz,2));
    ly = ly/mag;
    lz = lz/mag;
  }
  if(timezero(nx) != 0. || timezero(nz) != 0.) {  //We may generate a perpendicular vector by crossing the normal with the global y-hat vector.
    lx = -nz;
    ly = 0.;
    lz = nx;
    mag = sqrt(pow(lx,2)+pow(lz,2));
    lx = lx/mag;
    lz = lz/mag;
  }

  //Take the cross product (n-hat x l-hat) to create a unit vector perpendicular to both the normal and "longitudinal" unit vectors. [Bounce system's y-axis.]
  tx = ny*lz - nz*ly;
  ty = -nx*lz + nz*lx;
  tz = nx*ly - ny*lx;
  
  //Calculate the components of the given vector in the bounce system by taking the vector's projections on the bounce system's basis vectors.
  *bcx = gcx*nx + gcy*ny + gcz*nz;
  *bcy = gcx*tx + gcy*ty + gcz*tz;
  *bcz = gcx*lx + gcy*ly + gcz*lz;

  return;
}

void bsys2gsys(double bcx,double bcy,double bcz,double *gcx,double *gcy,double *gcz, double passnx, double passny, double passnz) {
  double nx,ny,nz,lx,ly,lz,tx,ty,tz,mag;
  
  if(neutron.xcode != -1 && neutron.xcode < 0) {  //Particle does not have a current intersection.
    printf("+ERROR: No current particle intersection in 'gsys2bsys'!\n");
    *gcx = -999999;  *gcy = -999999;  *gcz = -999999;
    return;
  }
  
  normal(&nx,&ny,&nz,passnx,passny,passnz); //Calculate normal at current point of intersection. [Bounce system x-axis.]
  
  //Generate a unit vector perpendicular to the unit normal. [Bounce system's z-axis.]
  if(timezero(ny) != 0. || timezero(nz) != 0.) {  //We may generate a perpendicular vector by crossing the normal with the global x-hat vector.
    lx = 0.;
    ly = nz;
    lz = -ny;
    mag = sqrt(pow(ly,2)+pow(lz,2));
    ly = ly/mag;
    lz = lz/mag;
  }
  if(timezero(nx) != 0. || timezero(nz) != 0.) {  //We may generate a perpendicular vector by crossing the normal with the global y-hat vector.
    lx = -nz;
    ly = 0.;
    lz = nx;
    mag = sqrt(pow(lx,2)+pow(lz,2));
    lx = lx/mag;
    lz = lz/mag;
  }
  

  //Take the cross product (n-hat x l-hat) to create a unit vector perpendicular to both the normal and "longitudinal" unit vectors. [Bounce system's y-axis.]
  tx = ny*lz - nz*ly;
  ty = -nx*lz + nz*lx;
  tz = nx*ly - ny*lx;
  
  //Calculate the components of the given vector in the global system by taking the vector's projections on the global system's basis vectors.
  *gcx = bcx*nx + bcy*tx + bcz*lx;   //x-hat,bounce = nx x-hat + ny y-hat + nz z-hat
  *gcy = bcx*ny + bcy*ty + bcz*ly;  //y-hat,bounce = tx x-hat + ty y-hat + tz z-hat
  *gcz = bcx*nz + bcy*tz + bcz*lz; //z-hat,bounce = lx x-hat + ly y-hay + lz z-hat

  return;
}

double intersection(int *eventcode, double traj[3][3], int *intrsct) {
  int i;
  int sectflag,conreg;
  double dt,tempdt;

  //Calculate the intersections with all surfaces connected to the current region.
  dt = 999999;
  if(CHECK == "ON") printf("--------------------\n");
  if(CHECK == "ON") printf(">>Neutron: %d | Current Region: %d | Elapsed Time BEFORE this step: %f\n",neutron.num,neutron.region,neutron.t);
  if((int)regions[neutron.region][0] == 1) {  //If a box-region.
    tempdt = boxsect(neutron.region,traj); //Get intersection time.
    if(CHECK == "ON") printf(">>Box-Region Intersection Time: %f\n",tempdt);
    if(tempdt < 0 && CONSOLE == "ON") printf("+WARNING: Negative time intersection from 'boxsect'!\n");
    if(tempdt < dt && tempdt > 0) {
      sectflag = 1; //A valid intersection was found.
      dt = tempdt; //This intersection time is the smallest so far.
      *intrsct = -1; //Intersection with current region.
      *eventcode = 1;
    }
  }
  if((int)regions[neutron.region][0] == 2) {  //If a cylinder region.
    tempdt = cylindersect(neutron.region,traj); //Get intersection time.
    if(CHECK == "ON") printf(">>Cylinder-Region Intersection Time: %f\n",tempdt);
    if(tempdt < 0 && CONSOLE == "ON") printf("+WARNING: Negative time intersection from 'cylindersect'!\n");
    if(tempdt < dt && tempdt > 0) {
      sectflag = 1; //A valid intersection was found.
      dt = tempdt; //This intersection time is the smallest so far.
      *intrsct = -1; //Intersection is with current region.
      *eventcode = 1;
    }
  }
  for(i = 0 ; i < 6 ; i++) {  //Except for region 0, i=0 checks a cut-plane which is not connected to the current region since the first entry in 'connex' is
                             //the region that the current region connects to through its own cut-plane.
    conreg = connex[neutron.region][i]; //Holds connecting region number to save typing below.
    if(connex[neutron.region][i] == -1) break; //Intersections with cut-planes to all connecting regions have been checked.
    if(cplanes[(int)connex[neutron.region][i]][4] == 3 && cplanes[neutron.region][4] == 2) {  //We have a main region connecting to a type-3 T-ing region.
      if(CHECK == "ON") printf(">>Cut-plane checking for region %d skipped since it connects via a T-junction.\n",connex[neutron.region][i]);
    }
    else {  //The connection is not a T-junction and cut-plane connection simultaneously.
      tempdt = planesect(connex[neutron.region][i],traj); //Looks for an intersection with a cut-plane associated with a region connected to current region unless
                                                         //it is a T-ing region in which case checking is done separately.
      if(CHECK == "ON") printf(">>Region %d's Cut-Plane Intersection Time: %f\n",connex[neutron.region][i],tempdt);
      if(tempdt < 0 && CONSOLE == "ON") printf("+WARNING: Negative time intersection from 'planesect'!\n");
      if(tempdt < dt && tempdt > 0) {
        sectflag = 1; //A valid intersection was found.
        dt = tempdt; //This intersection time is the smallest so far.
        *intrsct = connex[neutron.region][i]; //Particle has intersected with another region's cut-plane so we return that region number.
        *eventcode = 2;
      }
    }
  }

  tempdt = planesect(neutron.region,traj); //Check for intersection with the current region's cut-plane.
  if(CHECK == "ON") printf(">>Current Region's Cut-Plane Intersection Time: %f\n",tempdt);
  if(tempdt < 0 && CONSOLE == "ON") printf("+WARNING: Negative time intersection from 'planesect'!\n");
  if(tempdt < dt && tempdt > 0) {
    sectflag = 1; //A valid intersection was found.
    dt = tempdt; //This intersection time is the smallest so far.
    *intrsct = neutron.region; //Particle has intersected with the current region's cut-plane so we return the current region number.
    *eventcode = 2;
  }
  
  if((int)cplanes[neutron.region][4] == 3 || (int)cplanes[neutron.region][4] == 4) {  //Particle is in a region which T's into another region, i.e. ends in the side
                                                                                     //of a main region, so we must check for intersections with that actual region.
    if (CHECK == "ON") printf("+WARNING: Handling transport from region T'ing into another region, see line 2060.\n");
                                //beware unanticipated consequences of changing the below line from 'conreg = connex[neutron.region][0]' to 'conreg = neutron.region'. -Erik Lutz
    if((int)cplanes[neutron.region][4] == 3) conreg = neutron.region; //The current region begins at a main section so that the region with which we
                                                                                //should check for intersections is the first in the current region's 'connex' list.
    if((int)cplanes[neutron.region][4] == 4) conreg = connex[neutron.region][1]; //The current region ends at a main section so that the region with which we should
                                                                                //check for intersections is the second in the current region's 'connex' list.
    if((int)regions[conreg][0] == 1) tempdt = boxsect(conreg,traj); //Check for intersection with region that is T'ed into (box).
    if((int)regions[conreg][0] == 2) tempdt = cylindersect(conreg,traj); //Check for intersection with region that is T'ed into (cylinder).
    if(CHECK == "ON") {
      printf(">>Current (T-ing) Region's Junction Intersection Time (with region %d): %f\n",conreg,tempdt);
      printf(">>Position of current (T-ing) region's cut-plane: %f",basepoints[neutron.region][0]);
      printf(",%f,%f (t-junc. flag = %d)",basepoints[neutron.region][1],basepoints[neutron.region][2],(int)cplanes[neutron.region][6]);
      printf(" | Length of T-ing region: %f\n",regions[neutron.region][2]);
    }
    if(tempdt < 0 && CONSOLE == "ON") printf("+WARNING: Negative time T-junction intersection! [%d]\n",neutron.num);
    if(tempdt < dt && timezero(tempdt) > 0) {  //Sometimes on passage into a T-ing region a particle ends up a very tiny distance from the T-junction so that
                                              //a very small positive time is returned which we must ignore since the particle should not intersect with the
                                             //T-junction again just after passing through it.
      sectflag = 1; //A valid intersection was found.
      dt = tempdt; //This intersection time is the smallest so far and the last we're checking so that if this is the smallest then it is THE intersection.
      *intrsct = conreg; //Return the intersection code for an intersection with the main region that the particle is trying to enter,
      *eventcode = 16; //Return cut-plane intersection code so that we'll end up in 'cplaneshandling'.
    }
  }

  if(sectflag != 1) {  //No intersection was found.
    if(CONSOLE == "ON") printf("+WARNING: Particle trajectory possessing no intersection with geometry!\n");
    *eventcode = -1;
    return 0;
  }

  if(CHECK == "ON") printf(">>Intersection Time: %f\n",dt);
  return dt;
}

double planesect(int reg, double traj[3][3]) {
  double t1,t2,D;
  double bpx,bpy,bpz,thetar,phir,psir;
  double trajp[3][3];

  //Get the basepoint for the cut-plane's region and the orientation of the cut-plane.
  bpx = basepoints[reg][0];
  bpy = basepoints[reg][1];
  bpz = basepoints[reg][2];
  thetar = cplanes[reg][0];
  phir = cplanes[reg][1];
  psir = cplanes[reg][7];

  //Transform the trajectory coefficient matrix to the solving system.
  rotp(psir,thetar,phir,traj[0][0],traj[1][0],traj[2][0],&trajp[0][0],&trajp[1][0],&trajp[2][0]); //Transform the acceleration vector.
  rotp(psir,thetar,phir,traj[0][1],traj[1][1],traj[2][1],&trajp[0][1],&trajp[1][1],&trajp[2][1]); //Transform the velocity vector.
  rotp(psir,thetar,phir,traj[0][2]-bpx,traj[1][2]-bpy,traj[2][2]-bpz,&trajp[0][2],&trajp[1][2],&trajp[2][2]); //Transform the shifted initial position vector.

  //Return the smallest non-zero intersection.
  return solve2(trajp[2][0],trajp[2][1],trajp[2][2]);
}

double boxsect(int reg, double traj[3][3]) {
  double t1,t2,t3,t4,t;
  double bpx,bpy,bpz,thetar,phir,psir;
  double trajp[3][3];
  
  //Get the basepoint and orientation for the region.
  bpx = basepoints[reg][0];
  bpy = basepoints[reg][1];
  bpz = basepoints[reg][2];
  thetar = regions[reg][4];
  phir = regions[reg][5];
  psir = regions[reg][6];
  
  //Transform the trajectory coefficient matrix to the solving system.
  rotp(psir,thetar,phir,traj[0][0],traj[1][0],traj[2][0],&trajp[0][0],&trajp[1][0],&trajp[2][0]); //Transform the acceleration vector.
  rotp(psir,thetar,phir,traj[0][1],traj[1][1],traj[2][1],&trajp[0][1],&trajp[1][1],&trajp[2][1]); //Transform the velocity vector.
  rotp(psir,thetar,phir,traj[0][2]-bpx,traj[1][2]-bpy,traj[2][2]-bpz,&trajp[0][2],&trajp[1][2],&trajp[2][2]); //Transform the shifted initial position vector.
  
   //For a box region, the region's cut-plane is taken to be -z-face of the box and the next region's cut-plane is taken to be the next region's cut-plane. Since
  //cut-plane intersections are checked for elsewhere, we don't need to check for those intersections here.
  
  //Look for intersection with +x-Dim and -x-Dim planes.
  t1 = solve2(trajp[0][0],trajp[0][1],trajp[0][2]-regions[reg][1]/2); // '- plane'
  t2 = solve2(trajp[0][0],trajp[0][1],trajp[0][2]+regions[reg][1]/2); // '+ plane'
  
  //Look for intersections with +y-Dim and -y-Dim planes.
  t3 = solve2(trajp[1][0],trajp[1][1],trajp[1][2]-regions[reg][2]/2); // '- plane'
  t4 = solve2(trajp[1][0],trajp[1][1],trajp[1][2]+regions[reg][2]/2); // '+ plane'
  
  //Return the smallest non-zero time or zero.
  if(t1 < 0 && t2 < 0 && t3 < 0 && t4 < 0) return 0; //No non-zero roots.
  t = 999999;
  if(t1 < t && t1 > 0) t = t1;
  if(t2 < t && t2 > 0) t = t2;
  if(t3 < t && t3 > 0) t = t3;
  if(t4 < t && t4 > 0) t = t4;
  return t;
}

double cylindersect(int reg, double traj[3][3]) {
  double bpx,bpy,bpz,thetar,phir,psir;
  double c4, c3, c2, c1, c0; //The coefficients of the quartic c4 t^4 + c3 t^3 + c2 t^2 + c1 t + c0 = 0.
  double trajp[3][3]; //Holds the trajectory coefficient matrix in the solving system.

  //Get the orientation and basepoint for the cylinder region.
  bpx = basepoints[reg][0];
  bpy = basepoints[reg][1];
  bpz = basepoints[reg][2];
  thetar = regions[reg][4];
  phir = regions[reg][5];
  psir = regions[reg][6];

   //Transform the trajectory coefficient matrix to the solving system.
  rotp(psir,thetar,phir,traj[0][0],traj[1][0],traj[2][0],&trajp[0][0],&trajp[1][0],&trajp[2][0]); //Transform the acceleration vector.
  rotp(psir,thetar,phir,traj[0][1],traj[1][1],traj[2][1],&trajp[0][1],&trajp[1][1],&trajp[2][1]); //Transform the velocity vector.
  rotp(psir,thetar,phir,traj[0][2]-bpx,traj[1][2]-bpy,traj[2][2]-bpz,&trajp[0][2],&trajp[1][2],&trajp[2][2]); //Transform the shifted initial position vector.

  //Calculate the quartic intersection equation coefficients.
  c4 = pow(trajp[0][0],2) + pow(trajp[1][0],2);
  c3 = 2.*trajp[0][0]*trajp[0][1] + 2.*trajp[1][0]*trajp[1][1];
  c2 = pow(trajp[0][1],2) + pow(trajp[1][1],2) + 2.*trajp[0][0]*trajp[0][2] + 2.*trajp[1][0]*trajp[1][2];
  c1 = 2.*trajp[0][1]*trajp[0][2] + 2.*trajp[1][1]*trajp[1][2];
  c0 = pow(trajp[0][2],2) + pow(trajp[1][2],2) - pow(regions[reg][1]/2.,2); //Note that the constant coefficient includes the -r^2 term.

  //Return the value of the smallest positive real solution.
  return solve4(c4,c3,c2,c1,c0);
}

int move(double dt, double traj[3][3], int gcheck) {
  double bpx,bpy,bpz,thetar,phir,psir;
  double ssx,ssy,ssz,gcx,gcy,gcz;
  double discrep,angle;
  
  if(gcheck != 0 && gcheck != 1) {
    printf("+ERROR: Checking code sent to 'move' not recognized!\n");
    return -1;
  }
  
  neutron.t = neutron.t + dt; //Increment the time by the elapsed time along trajectory.
  
//The trajectory equations have the current particle position corresponding to t=0s so plugging in elapsed time gives final absolute coordinates.
  neutron.x = traj[0][0]*pow(dt,2) + traj[0][1]*dt + traj[0][2];
  neutron.y = traj[1][0]*pow(dt,2) + traj[1][1]*dt + traj[1][2];
  neutron.z = traj[2][0]*pow(dt,2) + traj[2][1]*dt + traj[2][2];
  neutron.vx = 2.*traj[0][0]*dt + traj[0][1];
  neutron.vy = 2.*traj[1][0]*dt + traj[1][1];
  neutron.vz = 2.*traj[2][0]*dt + traj[2][1];
  
//Now we must do a reality check to ensure that the particle is still inside the geometry. If it is only a "little bit" outside it should be placed on the surface.  
  //Check to make sure the particle is inside the current region.
  //Get the basepoint and orientation information about the surface.
  bpx = basepoints[neutron.region][0];
  bpy = basepoints[neutron.region][1];
  bpz = basepoints[neutron.region][2];
  thetar = regions[neutron.region][4];
  phir = regions[neutron.region][5];
  psir = regions[neutron.region][6];
  
  //Transform the particle's new position into the current region's solving system.
  rotp(psir,thetar,phir,neutron.x-bpx,neutron.y-bpy,neutron.z-bpz,&ssx,&ssy,&ssz);
  
  if(regions[neutron.region][0] == 2) {  //The particle is currently in a cylinder-region.
    discrep = regions[neutron.region][1]/2. - sqrt(pow(ssx,2) + pow(ssy,2)); //Calculate the difference between the cylinder radius and the particle radius.
    
    //Check to see if the particle has left the sides of the region.
    if(discrep < -GZERO) {  //We'll consider the particle too far outside the geometry and abandon it.
      if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [discrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
	  event(-1,0,0,0,0);
      return -1;
    }
    
      //Check to see if the particle has left through the end. (If only a "little bit" on a cut-plane intersection it will be corrected below.)
     //Note that this is a rough check since a cylindrical region is only the specified length from basepoint to basepoint. Here we assume that
    //the "overhang" won't be larger than the diameter of the tube.
    if(ssz < -regions[neutron.region][2] - regions[neutron.region][1] || ssz > regions[neutron.region][2] + regions[neutron.region][1]) {  //Neutron has escaped
                                                                                                                                           //through a cut-plane.
      if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped through a cut-plane! [ssz = %e | maxssz = %e | region = %d]\n",neutron.num,ssz,regions[neutron.region][2] + regions[neutron.region][1],neutron.region);
      event(-1,0,0,0,0);
      return -1;
    }
    
    if(discrep < 0 && discrep >-GZERO) {  //Particle is just a "little bit" outside of geometry.
      angle = atan2(ssy,ssx); //Find the angle that the particle's position makes relative to the x-axis.
      //Place the particle on the surface.
      ssx = regions[neutron.region][1]/2.*cos(angle);
      ssy = regions[neutron.region][1]/2.*sin(angle);
    }
  }
  
  if(regions[neutron.region][0] == 1) {  //The particle is currently in a box-region.
    //Check +z-Dim.
    discrep = regions[neutron.region][3] - ssz;
    if(discrep < -GZERO) {  //We'll consider the particle too far outside the geometry and abandon it.
      if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [+zdiscrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
      event(-1,0,0,0,0);
      return -1;
    }
    if(discrep < 0 && discrep > -GZERO) ssz = regions[neutron.region][3]; //Place the particle on the z-face.
      
    //Check y-Dim.
    if(ssy > 0) {  //We should check the +y-Dim.
      discrep = regions[neutron.region][2]/2. - ssy;
      if(discrep < -GZERO) {  //we'll consider the particle too far outside the geometry and abandon it.
        if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [+ydiscrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
        event(-1,0,0,0,0);
        return -1;
      }
      if(discrep < 0 && discrep > -GZERO) ssy = regions[neutron.region][2]/2.; //Place the particle on the +y-face.
    }
    if(ssy < 0) {  //We should check the -y-Dim.
      discrep = regions[neutron.region][2]/2. + ssy;
      if(discrep < -GZERO) {  //we'll consider the particle too far outside the geometry and abandon it.
        if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [-ydiscrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
        event(-1,0,0,0,0);
        return -1;
      }
      if(discrep < 0 && discrep > -GZERO) ssy = -regions[neutron.region][2]/2.; //Place the particle on the -y-face.
    }
      
    //Check x-Dim.
    if(ssx > 0) {  //We should check the +x-Dim.
      discrep = regions[neutron.region][1]/2. - ssx;
      if(discrep < -GZERO) {  //we'll consider the particle too far outside the geometry and abandon it.
        if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [+xdiscrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
        event(-1,0,0,0,0);
        return -1;
      }
      if(discrep < 0 && discrep > -GZERO) ssx = regions[neutron.region][1]/2.; //Place the particle on the +x-face.
    }
    if(ssx < 0) {  //We should check the -x-Dim.
      discrep = regions[neutron.region][1]/2. + ssx;
      if(discrep < -GZERO) {  //we'll consider the particle too far outside the geometry and abandon it.
        if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the geometry! [-xdiscrep = %e | region = %d]\n",neutron.num,discrep,neutron.region);
        event(-1,0,0,0,0);
        return -1;
      }
      if(discrep < 0 && discrep > -GZERO) ssx = -regions[neutron.region][1]/2.; //Place the particle on the -x-face.  
    }
  }
  
  
  //Transform the adjusted position back to global coordinates, correcting the current particle position.
  rota(psir,thetar,phir,ssx,ssy,ssz,&gcx,&gcy,&gcz);
  neutron.x = gcx + bpx;
  neutron.y = gcy + bpy;
  neutron.z = gcz + bpz;
  
  
  if(gcheck == 1) {  //We must check for a valid cut-plane intersection.
    //Get basepoint and orientation information for the cut-plane.
    bpx = basepoints[neutron.xcode][0];
    bpy = basepoints[neutron.xcode][1];
    bpz = basepoints[neutron.xcode][2];
    thetar = cplanes[neutron.xcode][0];
    phir = cplanes[neutron.xcode][1];
    psir = cplanes[neutron.xcode][7];

    //Transform the particle's new position into the cut-plane's solving system.
    rotp(psir,thetar,phir,neutron.x-bpx,neutron.y-bpy,neutron.z-bpz,&ssx,&ssy,&ssz);
        
    discrep = ssz; //Particle should be on the cut-plane so solving-system z-coordinate should be zero.
    if(discrep < -GZERO || discrep > GZERO) {  //We'll consider the particle too far away from the cut-plane and abandon it.
      if(CONSOLE == "ON") printf("+WARNING: Particle number %d escaped from the current region! [discrep = %e | region = %d | xcode = %d]\n",neutron.num,discrep,neutron.region, neutron.xcode);
      event(-1,0,0,0,0);
      return -1;
    }
    if(discrep < GZERO && discrep > -GZERO) ssz = 0.; //Place the particle on the cut-plane.

    //Transform the adjusted position back to global coordinates, correcting the current particle position.
    rota(psir,thetar,phir,ssx,ssy,ssz,&gcx,&gcy,&gcz);
    neutron.x = gcx + bpx;
    neutron.y = gcy + bpy;
    neutron.z = gcz + bpz;
  }
  return 0;
}

void trajmatrix(double traj [3][3]) {

  //Assemble trajectory matrix (relative to global coordinate system) for a particle moving with accelerations given in 'regions'.
  traj[0][0] = -(1./2.)*neutron.spin[2]*(1./2.)*nGAMMA*hBAR*rparams[neutron.region][0]/nMASS; //1/2 ax
  traj[0][1] = neutron.vx; //vox
  traj[0][2] = neutron.x;  //xo
  traj[1][0] = -G/2.-(1./2.)*neutron.spin[2]*(1./2.)*nGAMMA*hBAR*rparams[neutron.region][1]/nMASS; //1/2 ay
  traj[1][1] = neutron.vy; //voy
  traj[1][2] = neutron.y; //yo
  traj[2][0] = -(1./2.)*neutron.spin[2]*(1./2.)*nGAMMA*hBAR*rparams[neutron.region][2]/nMASS; //1/2 az
  traj[2][1] = neutron.vz; //voz
  traj[2][2] = neutron.z; //zo
  
  return;
}

int tjunc(int tingreg) {
  int flag;
  double xp,yp,zp;
  double deltar;
  double thetar,phir,psir;

  flag = 0;
  thetar = regions[tingreg][4]; //Get theta-orientation of T-ing region.
  phir = regions[tingreg][5]; //Get phi-orientation of T-ing region.
  psir = regions[tingreg][6]; //Get psi-orientation of T-ing region.
  rotp(psir,thetar,phir,neutron.x-basepoints[tingreg][0],neutron.y-basepoints[tingreg][1],neutron.z-basepoints[tingreg][2],&xp,&yp,&zp); //Transform the particle's
                                                                                                                                        //position to the T-ing
                                                                                                                                       //region's solving system.
  if(cplanes[tingreg][4] == 4) {  //The T-ing region ENDS in a main region so that the T-ing region's length must be taken into account when checking for a
                                 //T-junction intersection from inside the main region. Note that in this case the T-ing region should be expanded into the main 
                                //region.
    if(regions[tingreg][0] == 1) {  //T-ing region is a box-region.
      if(fabs(xp) < regions[tingreg][1]/2. && fabs(yp) < regions[tingreg][2]/2. && zp < regions[tingreg][3]) return 1; //The particle is at the T-junction.
      else return 0; //The particle in not inside the bounds of the T-ing region and so is not at a T-junction.
    }
    if(regions[tingreg][0] == 2) {  //T-ing region is a cylinder-region.
      deltar = sqrt(pow(xp,2) + pow(yp,2)); //Distance of the particle from the axis of the T-ing region.
      if(deltar < regions[tingreg][1]/2. && zp < regions[tingreg][2]) return 1; //The particle is inside the bounds of the T-ing region.
      else return 0; //The particle is not inside the bounds of the T-ing region and so is not at a T-junction.
    }
  }
  
  if(cplanes[tingreg][4] == 3) {  //The T-ing region BEGINS in a main region so that the length of the T-ing region should not be taken into account when checking
                                 //for a T-junction intersection from inside the main region. Note that in this case the T-ing region will not be expanded into the
                                //main region since the particle is currently in the main region.
   if(regions[tingreg][0] == 1) {  //T-ing region is a box-region.
      if(fabs(xp) < regions[tingreg][1]/2. && fabs(yp) < regions[tingreg][2]/2. && zp > -1.*cplanes[tingreg][5]) return 1; //The particle is at the T-junction.
      else return 0; //The particle in not inside the bounds of the T-ing region and so is not at a T-junction.
    }
    if(regions[tingreg][0] == 2) {  //T-ing region is a cylinder-region.
      deltar = sqrt(pow(xp,2) + pow(yp,2)); //Distance of the particle from the axis of the T-ing region.
      if(deltar < regions[tingreg][1]/2. && zp > -1.*cplanes[tingreg][5]) return 1; //The particle is inside the bounds of the T-ing region.
      else return 0; //The particle is not inside the bounds of the T-ing region and so is not at a T-junction.
    }
  }
  
  return -1;
}

void cplaneshift(int mvcode, int cutplane) {
  double ssdeltaz,deltax,deltay,deltaz;

  if(mvcode == 1 && (int)cplanes[cutplane][6] == 1) {  //A region should never be expanded twice in a row.
    printf("+ERROR: Attemp to expand region %d twice in a row using 'cplaneshift'!\n",cutplane);
    return;
  }
  
  if(mvcode == -1 && (int)cplanes[cutplane][6] == -1) {  //A region should never be shrunk twice in a row.
    printf("+ERROR: Attemp to shrink region %d twice in a row using 'cplaneshift'!\n",cutplane);
    return;
  }
  
  if((int)cplanes[cutplane][4] == 4) {  //T-ing region ends in a main section so we just need to change the region's length.
    if((int)regions[cutplane][0] == 1) regions[cutplane][3] = regions[cutplane][3] + mvcode*cplanes[cutplane][5]; //Change z-dim for a box-region.
    if((int)regions[cutplane][0] == 2) regions[cutplane][2] = regions[cutplane][2] + mvcode*cplanes[cutplane][5]; //Change length for a cylinder-region.
    cplanes[cutplane][6] = mvcode; //Set the region-shift flag appropriately.
    return;
  }
  
  if((int)cplanes[cutplane][4] == 3) {  //T-ing region begins in a main section so we need to adjust basepoint and change length.
    //Transform the basepoint offset along solving system z axis to the global system and adjust the basepoint.
    ssdeltaz = -1.*mvcode*cplanes[cutplane][5]; //Translate.
    rota(cplanes[cutplane][7],cplanes[cutplane][0],cplanes[cutplane][1],0.,0.,ssdeltaz,&deltax,&deltay,&deltaz); //Transform shift.
  
    //Apply Shift,
    basepoints[cutplane][0] = basepoints[cutplane][0] + deltax;
    basepoints[cutplane][1] = basepoints[cutplane][1] + deltay;
    basepoints[cutplane][2] = basepoints[cutplane][2] + deltaz;
  
    if(regions[cutplane][0] == 1) regions[cutplane][3] = regions[cutplane][3] + mvcode*cplanes[cutplane][5]; //Alter length (z-dim) of box region to compensate
                                                                                                            //for shift in basepoint.
    if(regions[cutplane][0] == 2) regions[cutplane][2] = regions[cutplane][2] + mvcode*cplanes[cutplane][5]; //Alter length of cylinder region to compensate
                                                                                                            //for shift in basepoint.
    cplanes[cutplane][6] = mvcode; //Set the region-shift flag appropriately.
    return;
  }
  
  printf("+ERROR: 'cplaneshift' called without a valid internal cut-plane special handling flag!\n");
  return;
}

int cplanehandling(void) {
  int regf,cpcode,icpcode,icrcode,jcode,ecode,bcheck;
  double pangle,flipeff,radius,psir,thetar,phir,pangleloc;
  double deltakE, minusdeltakEneV;
  double vn,vyb,vzb,bpx,bpy,bpz,xp,yp,zp,vxloc,vyloc,vzloc;
  double traj[3][3],dt;
  double xb,yb,zb,raperature,rparticle;
  double passnx,passny,passnz;
  double vninitial, vnfinal, vyi,vzi,vyf,vzf,nx,ny,nz,nxperturb,nyperturb,nzperturb;
  int cyclecount;
  int i;
  double vnlocalinital,vnlocalfinal;
  int roughnessmodel, passthroughregf;
  double passthroughx,passthroughy,passthroughz,passthroughthickness;
  
  if(neutron.xcode < 0) {  //Current intersection is with a region or there is no current intersection.
    printf("+ERROR: 'cplanehandling' called without a cut-plane intersection! Particle number %d abandoned!\n",neutron.num);
    return -1;
  }

  //Get the special-handling codes (if any) and determine which region the particle should end up in if it is passed through the cut-plane.
  cpcode = (int)cplanes[neutron.xcode][2]; //User specified special-handling code.
  icpcode = (int)cplanes[neutron.xcode][4]; //Internal special-handling code for the cut-plane.
  if(neutron.region == neutron.xcode) regf = connex[neutron.region][0]; //Particle is leaving current region through that region's cut-plane so the
                                                                       //particle will end end up in the region indicated first in 'connex' for current region.
  if(neutron.region != neutron.xcode) regf = neutron.xcode; //Particle is entering the region associated with the cut plane so it will end up in the
                                                           //cut-plane's region (or is passing through a T-juntion).
  
  if(icpcode == 1) {  //Cut-plane is connecting a box region and a cylinder region. Note that if no lip is encountered, code passes on to check for special-handling.
    if (neutron.region != neutron.xcode) {
        bpx = basepoints[regf][0];
        bpy = basepoints[regf][1];
        bpz = basepoints[regf][2];
        }
    if (neutron.region == neutron.xcode) {
        bpx = basepoints[regf][0];
        bpy = basepoints[regf][1];
        bpz = basepoints[regf][2] + regions[regf][2]; // VERY UNROBUST! only works with cylinders pointed directly in the z direction. Otherwise the math breaks
                                                     // didnt feel like figuring out Euler angles just yet.
        }
    psir = cplanes[neutron.xcode][7];
    thetar = cplanes[neutron.xcode][0];
    phir = cplanes[neutron.xcode][1];
    
    if(regions[neutron.region][0] == 1 && regions[regf][0] == 2) {  //We are going from a box-region into a cylinder-region.
      //Calculate the distance from the axis of the guides to the particle's location (which is somewhere on the cut-plane and so is the particle's radius).
      radius = sqrt(pow(bpx-neutron.x,2)+pow(bpy-neutron.y,2)+pow(bpz-neutron.z,2));
      if(radius >= regions[regf][1]/2.) {  //The particle is outside the radius of the cylindrical guide and must be bounced.
        if(neutron.xcode == neutron.region) bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],1,0,0,0); //Bounce is being called for reflection off the current
                                                                                                       //region's cut-plane. The models/values and potential will
                                                                                                      //be those associated with current (box) region.
        else bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],-1,0,0,0); //The normal used in 'bounce' must be reversed since the particle has hit a cut-plane
                                                                         //associated with an adjoining region. Note that the models/values and potential used will
                                                                        //be those associated with the cylinder-region.
        if(bcheck == 0) return 0; //The particle bounced off the lip and has been given a new bounced velocity. Note that no further special handling is done.
        if(bcheck == -1) return -1; //An error occurred during the bounce so that the particle has been lost unphysically.
        if(bcheck == 1 || bcheck == 2) return 1; //The particle was lost physically.
        printf("+ERROR: Unrecognized code from 'bounce'! [Particle#: %d]\n",neutron.num);
        return -1;
      }   
    }
    if(regions[neutron.region][0] == 2 && regions[regf][0] == 1) {  //We are going from a cylinder-region into a box-region.
      rotp(psir,thetar,phir,neutron.x-bpx,neutron.y-bpy,neutron.z-bpz,&xp,&yp,&zp); //Transform particle's position into the solving system.
      if(xp >= regions[regf][1]/2. || xp <= -regions[regf][1]/2. || yp >= regions[regf][2]/2. || yp <= -regions[regf][2]/2.) {  //Particle must be bounced.
        if(neutron.xcode == neutron.region) bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],1,0,0,0); //Bounce is being called for reflection off the current
                                                                                                       //region's cut-plane. The models/values and potential will
                                                                                                      //be those associated with the cylinder region.
        else bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],-1,0,0,0); //The normal used in 'bounce' must be reversed since the particle has hit a cut-plane
                                                                         //associated with an adjoining region. Note that the models/values and potential used will
                                                                        //be those associated with the box-region.
        if(bcheck == 0) return 0; //The particle bounced off the lip and has been given a new bounced velocity. Note that no further special handling is done.
        if(bcheck == -1) return -1; //An error occurred during the bounce so that the particle has been lost unphysically.
        if(bcheck == 1 || bcheck == 2) return 1; //The particle was lost physically.
        printf("+ERROR: Unrecognized code from 'bounce'! [Particle#: %d]\n",neutron.num);
        return -1;
      }   
    }
  }
  
  if(cpcode == 1) {  //Cut-plane is 100% absorptive.
    if(EVENTS == "ON") event(8,0,0,0,0); //Log an event.
    if((int)rparams[neutron.xcode][10] != 0) return 9; // if region is defined to be a detector, log a detector event
    return 1; //Return the physical loss code so that the particle will be dropped.
  }

  if(cpcode == 3) {  //Spin-flip resonance.
    flipeff = 0.996; //Set the spin-flip efficiency.
    if(grn() < flipeff) neutron.spin[2] = -neutron.spin[2]; //Check for spin-flip and flip the spin if there is one.
    if(EVENTS == "ON") event(10,0,0,0,0); //Log an event.
    cpcode = -1; //Set the cut-plane code to -1 so that the particle will be passed through below.
  }

//  if(cpcode == 4) {  //Cut-plane is a detector that also records angular information.
//    pangle = atan2(sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)),neutron.vz); //Calculate particle's polar angle (relative to global system) at the cut-plane.
//    event(9,neutron.xcode,pangle,sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2)),0); //Record a detection event keyed to cut-plane number.
//    cpcode = -1; //Set the cut-plane code to -1 so that the particle will be passed through below. (Cut-plane detectors are passive detectors!)
//  }
if(cpcode == 4) {  //Cut-plane is a detector that also records angular information.
  gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vn,&vyb,&vzb,0,0,0);
	thetar = regions[neutron.region][4];
	phir = regions[neutron.region][5];
	psir = regions[neutron.region][6];
	rotp(psir,thetar,phir,neutron.vx,neutron.vy,neutron.vz,&vxloc,&vyloc,&vzloc);
	pangle = atan2(sqrt(pow(neutron.vx,2)+pow(neutron.vy,2)),neutron.vz);
    pangleloc = atan2(sqrt(pow(vxloc,2)+pow(vyloc,2)),vzloc); //Calculate particle's polar angle (relative to global system) at the cut-plane.
    //The 'detectors' array is passed to this functions so that cut-plane detections may be written to channels in that file if desired.
      //Increment the correct bin in 'detectors[][0]' for the measured polar angle.
      //for(i=0 ; i < 3000 ; i++) if(pangle >= PI/(2.*3000.) * i && pangle < PI/(2.*3000.) * (i+1)) detectors[i][0]++;
      //Increment the correct bin in 'detectors[][1]' for the corresponding time of detection.
      //for(i=0 ; i < 3000 ; i++) if(neutron.t >= i/10. && neutron.t < (i+1)/10.) detectors[i][1]++;
    event(18,neutron.xcode,cos(pangleloc),.5*nMASS*(pow(neutron.vx,2)+pow(neutron.vy,2)+pow(neutron.vz,2))/neV2J,.5*nMASS*pow(vn,2)/neV2J); //Record a detection event keyed to cut-plane number.
    cpcode = -1; //Set the cut-plane code to -1 so that the particle will be passed through below. (Cut-plane detectors are passive detectors!)
  }
  
  if(cpcode == 5) {  //Cut-plane is at a lip between two dissimilar-sized same-type guides.
    if(regions[neutron.region][0] != regions[regf][0]) {
      printf("+WARNING: Special-handling code 5 used at a connection involving two regions of different type!\n");
    }
    if((regions[neutron.region][4] != regions[regf][4] || regions[neutron.region][5] != regions[regf][5] || regions[neutron.region][6] != regions[regf][6]) && CONSOLE == "ON") {
      printf("+WARNING: A lip (special-handling code 5) has been specified between two regions that have diferent orientations!\nNeutron is leaving region %d and entering region %d.\n", neutron.region, regf);
    }
    if(regions[regf][0] == 2) {  //The connection is between cylinder-regions.
      if(regions[neutron.region][1] < regions[regf][1]) cpcode = -1; //Nothing needs to be done if the particle is going from the smaller to the larger region.
      else {  //Now we must check whether to bounce or pass the particle.
        //Calculate the distance from the axis of the guides to the particle's location (which is somewhere on the cut-plane and so is the particle's radius).
        radius = sqrt(pow(basepoints[neutron.xcode][0]-neutron.x,2)+pow(basepoints[neutron.xcode][1]-neutron.y,2)+pow(basepoints[neutron.xcode][2]-neutron.z,2));
        if(radius < regions[regf][1]/2.) cpcode = -1; //The particle is whithin the radius of the smaller guide and may pass through.
        else {  //The particle is outside the radius of the smaller guide and must be bounced
          if(neutron.xcode == neutron.region) bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],1,0,0,0); //Bounce is being called for reflection off the current
                                                                                                         //region's cut-plane. The models and values for a bounce
                                                                                                        //will be those associated with the cut-plane's region.
          else bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],-1,0,0,0); //The normal used in 'bounce' must be reversed since the particle has hit a cut-plane
                                                                           //associated with an adjoining region. Note that the models and values for the bounce will
                                                                          //still be those associated with the cut-plane's region.
          if(bcheck == 0) return 0; //The particle bounced off the lip and has been given a new bounced velocity.
          if(bcheck == -1) return -1;  //An error occurred during the bounce so that the particle has been lost unphysically.
          if(bcheck == 1 || bcheck == 2) return 1; //The particle was lost physically.
          printf("+ERROR: Unrecognized code from 'bounce'! [Particle#: %d\n",neutron.num);
          return -1;
        }   
      }
    }
    if(regions[regf][0] == 1) {  //The connection is between box-regions.
      bpx = basepoints[neutron.xcode][0];
      bpy = basepoints[neutron.xcode][1];
      bpz = basepoints[neutron.xcode][2];
      psir = cplanes[neutron.xcode][7];
      thetar = cplanes[neutron.xcode][0];
      phir = cplanes[neutron.xcode][1];
      rotp(psir,thetar,phir,neutron.x-bpx,neutron.y-bpy,neutron.z-bpz,&xp,&yp,&zp); //Transform particle's position into the solving system.
      if(xp >= regions[regf][1]/2. || xp <= -regions[regf][1]/2. || yp >= regions[regf][2]/2. || yp <= -regions[regf][2]/2.) {  //Particle must be bounced.
        if(neutron.xcode == neutron.region) bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],1,0,0,0); //Bounce is being called for reflection off the current
                                                                                                       //region's cut-plane. The models/values and potential will
                                                                                                      //be those associated with the cylinder region.
        else bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],-1,0,0,0); //The normal used in 'bounce' must be reversed since the particle has hit a cut-plane
                                                                         //associated with an adjoining region. Note that the models/values and potential used will
                                                                        //be those associated with the box-region.
        if(bcheck == 0) return 0; //The particle bounced off the lip and has been given a new bounced velocity.
        if(bcheck == -1) return -1; //An error occurred during the bounce so that the particle has been lost unphysically.
        if(bcheck == 1 || bcheck == 2) return 1; //The particle was lost physically.
        printf("+ERROR: Unrecognized code from 'bounce'! [Particle#: %d]\n",neutron.num);
        return -1;
      }
      cpcode = -1; //The particle is whithin the guide that it is entering and so may be passed through.
    }
  }
  
  if (cpcode == 7) { //cutplane has a spherical aperature (must have this positioned before cpcode = 6 and -1 so that they can be called depending on particle position). Made by Erik Lutz
    
    //printf("We got a sevener!\n");
    
    //printf("(%f,%f,%f)\n",neutron.x,neutron.y,neutron.z);
    
    xb = neutron.x - basepoints[neutron.xcode][0]; //translate to basepoint
    yb = neutron.y - basepoints[neutron.xcode][1];
    zb = neutron.z - basepoints[neutron.xcode][2];
    gsys2bsys(xb,yb,zb,&xb,&yb,&zb,0,0,0); //transform to bounce frame
    
    //printf("(%f,%f,%f)\n",xb,yb,zb);
    
    raperature = cplanes[neutron.xcode][3]; //get radius of aperature
    rparticle = sqrt(pow(yb,2)+pow(zb,2)); //calculate distance from center of cutplane
    
    //printf("(%f,%f)\n",raperature,rparticle);
    
    if (2*raperature > regions[neutron.xcode][1]) raperature = regions[neutron.xcode][1]/2.; //aperature cannot be bigger than region
    
    if (rparticle <  raperature) cpcode = -1; //particle is within aperature so pass it through
    if (rparticle >= raperature) cpcode = 6; //particle is outside of aperature so bounce it back
    
    if ((rparticle <  raperature) && (rparticle >= raperature)){ //Error
      if(CONSOLE == "ON") printf("Particle lost unphysically on an aperature cutplane");
      return -1;
    }
      
  }
  
  if (cpcode == 6) { //The cutplane should behave just like the walls of its associated region. Made by Erik Lutz
    if(neutron.xcode == neutron.region) bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],1,0,0,0); //Bounce is being called for reflection off the current
    //region's cut-plane. The models/values and potential will
    //be those associated with current (box) region.
    else bcheck = bounce(neutron.xcode,rparams[neutron.xcode][6],-1,0,0,0); //The normal used in 'bounce' must be reversed since the particle has hit a cut-plane
    //associated with an adjoining region. Note that the models/values and potential used will
    //be those associated with the cylinder-region.
    if(bcheck == 0) return 0; //The particle bounced off the lip and has been given a new bounced velocity. Note that no further special handling is done.
    if(bcheck == -1) return -1; //An error occurred during the bounce so that the particle has been lost unphysically.
    if(bcheck == 1 || bcheck == 2) return 1; //The particle was lost physically.
  }
  
  if (cpcode == 8) { //cutplane becomes 100% absorptive after neutron passes through it
    
    if (cplanes[neutron.region][3] == neutron.num){
        if(EVENTS == "ON") event(8,0,0,0,0); //Log an event.
        if((int)rparams[neutron.xcode][10] != 0) return 9; // if region is defined to be a detector, log a detector event
        return 1; //Return the physical loss code so that the particle will be dropped.
    }
    else {
        cplanes[neutron.region][3] = neutron.num;
        cpcode = -1;
    }
    
  }


  if(cpcode == -1 || cpcode == 2 || cpcode == 9) {  //There is no special handling, but we must take care of any changes in Fermi potential when passing from one region to the next.
                                    //OR the surface has a defined surface roughness, implemented by Erik Lutz
                                    //OR the plane is skewd off-perpendicular from the region, implemented by Erik Lutz, use with care, not robust at all, no trapping for unphysical or unrealistic angles
    if(neutron.xcode == 0) {  //Particle is incident on the start-region's cut-plane, whose behavior should have been specified by a special-handling code.
      printf("+WARNING: Particle number %d escaped through the start-region's cut-plane since no definite handling instruction was given!\n",neutron.num);
      return -1; //Return the error code so that particle will be dropped.
    }
    if(rparams[neutron.region][7] == rparams[regf][7]) {  //There is no difference in Fermi potential between the regions so simply pass the particle through.
      neutron.region = regf; //Set the particle's region so that intersections with the correct parts of the geometry will be sought.
      return 0;
    }
    else {
      
      passnx = 0;
      passny = 0;
      passnz = 0;
      vn     = 0;
      
      deltakE = (rparams[neutron.region][7] - rparams[regf][7])*neV2J; //The difference in Fermi potential between the region the particle is leaving (defined
                                                                      //at its cut-plane) and the region it is entering (defined at its cut-plane) gives the
                                                                     //required energy shift, i.e deltakE = PEi - PEf
        
      minusdeltakEneV = -(rparams[neutron.region][7] - rparams[regf][7]);

      gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vninitial,&vyi,&vzi,0,0,0);
      //if (cpcode == 2) printf("\nInitial velocity in bounce frame = (%f, %f, %f).\n",vninitial,vyi,vzi);
      
      
      if(neutron.xcode == neutron.region){ //particle is interacting with the cutplane of the region its in, normal velocity should be negative
        
        if (cpcode == 2) {//the surface has a defined roughness
          cyclecount = 0;
          while (vn >= 0) {
            roughnessmodel = perturbednormal(&passnx,&passny,&passnz);
            gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vn,&vyb,&vzb,passnx,passny,passnz); //Transform the particle's velocity into the bounce system (where, since the intersection
                                                                                           //is with a cut-plane, the +x-axis will be normal to the cut-plane and point into the
                                                                                          //cut-plane's region) so that normal velocity component shifts may be calculated below.
            cyclecount++;
            if (cyclecount == 1000) {
              if (CONSOLE == "ON") printf("+ERROR: Tried unsuccessfully to find a suitable random normal.\n");
              return -1;
            }
          }
        }
        else gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vn,&vyb,&vzb,passnx,passny,passnz);
        bcheck = bounce(regf,minusdeltakEneV,1,passnx,passny,passnz); //Bounce is being called for reflection off the current region's cut-plane.
                                                                     //The models/values used for the bounce will be those of the region into which
                                                                    //the particle is heading (with the material potential taken as that region's
                                                                   //bulk potential).
      }
      
      
      
      else { //particle is interacting with a connecting region's cutplane, normal velocity should be positive
        
        if (cpcode == 2) {//the surface has a defined roughness
          gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vninitial,&vyi,&vzi,0,0,0);
          cyclecount = 0;
          while (vn <= 0) {
            roughnessmodel = perturbednormal(&passnx,&passny,&passnz);
            gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vn,&vyb,&vzb,passnx,passny,passnz); //Transform the particle's velocity into the bounce system (where, since the intersection
                                                                                           //is with a cut-plane, the +x-axis will be normal to the cut-plane and point into the
                                                                                          //cut-plane's region) so that normal velocity component shifts may be calculated below.
            //printf("cycling to find suitable velocity.\n");
            cyclecount++;
            if (cyclecount == 1000) {
              if (CONSOLE == "ON") printf("+ERROR: Tried unsuccessfully to find a suitable random normal.\n");
              return -1;
            }
          }
        }
        if (cpcode == -1) gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vn,&vyb,&vzb,passnx,passny,passnz);
        bcheck = bounce(regf,minusdeltakEneV,-1,passnx,passny,passnz); //The normal used in 'bounce' must be reversed since particle may reflect
                                                                      //off a cut-plane associated with an
                                                                     //adjoining region. The models/values for the region the particle is trying to enter will be used for the
                                                                    //bounce (with the material potential taken as that region's bulk potential).
        
      }
      
      
      gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vnfinal,&vyf,&vzf,0,0,0);

      
      if(bcheck == 0) { //The particle bounced off the cut-plane and has been given a new bounced velocity.
        if(EVENTS == "ON") event(19,0,0,0,0);
        if ((vninitial > 0 && vnfinal > 0) || (vninitial < 0 && vnfinal < 0)) { //even though particle is bounced it still goes into next region because of roughness
            
            if (roughnessmodel == 0 || roughnessmodel == 1) { //passthrough's are not allowed so return 2 so that particle will reinteract with surface.
                return 2;
            }
        
            else if ((roughnessmodel == 2 || roughnessmodel == 3) && rparams[neutron.region][7] == 0.) { //we must pass the particle through the next region
            
                passthroughregf = neutron.region + 2*(regf - neutron.region); //define region to pass neutron to
                if (neutron.region == neutron.xcode) {passthroughthickness = -1*regions[regf][2];} //define thickness of region to pass through
                else {passthroughthickness = regions[regf][2];}
                
                
                //calculate change in particle position as is travels through pass through region, in bounce frame
                passthroughx = passthroughthickness;
                passthroughy = (passthroughthickness/vnfinal)*vyf;
                passthroughz = (passthroughthickness/vnfinal)*vzf;
                
                //transform neutron postion into bounce system
                xb = neutron.x - basepoints[neutron.xcode][0]; //translate to basepoint
                yb = neutron.y - basepoints[neutron.xcode][1];
                zb = neutron.z - basepoints[neutron.xcode][2];
                gsys2bsys(xb,yb,zb,&xb,&yb,&zb,0,0,0); //transform to bounce frame
                
                //move neutron
                xb += passthroughx;
                yb += passthroughy;
                zb += passthroughz;
                bsys2gsys(xb,yb,zb,&xb,&yb,&zb,0,0,0); //transform to global frame
                
                //translate back from basepoint
                xb += basepoints[neutron.xcode][0];
                yb += basepoints[neutron.xcode][1];
                zb += basepoints[neutron.xcode][2];
                
                
                //set neutron position
                neutron.x = xb;
                neutron.y = yb;
                neutron.z = zb;
                
                neutron.region = passthroughregf; //set new region
                
                return 0;

            }
            
        else return 2;

        }
        return 0;
      }
      
      if(bcheck == 1) {  //The particle penetrated into the region and so it should be passed through the cut-plane with an energy shift.
        if(EVENTS == "ON") event(20,0,0,0,0);
        if((2.*deltakE/nMASS + pow(vn,2)) < 0) {  //After the energy shift the particle speed will be imaginary!
          printf("+ERROR: Particle number %d arrived inside a medium with an imaginary velocity component! [Reg: %d, xcode: %d]\n",neutron.num,neutron.region,neutron.xcode);
          return -1; //Abandon the particle.
        }
        
        
        vnlocalinital = vn;
        vn = vn/fabs(vn)*sqrt(2.*deltakE/nMASS + pow(vn,2)); //Calculate the particle's new normal velocity component.
        vnlocalfinal = vn;
        bsys2gsys(vn,vyb,vzb,&neutron.vx,&neutron.vy,&neutron.vz,passnx,passny,passnz); //Transform the particle's new velocity back into the global system.
        
        gsys2bsys(neutron.vx,neutron.vy,neutron.vz,&vnfinal,&vyf,&vzf,0,0,0);
        
        if ((vninitial > 0 && vnfinal < 0) || (vninitial < 0 && vnfinal > 0)) {
            neutron.region = regf;
            return 2;
        }
        
        else neutron.region = regf; //Pass the particle through the cut-plane.

        if(EVENTS == "ON") event(15,0,0,0,0); //Log an event.
        return 0;
      }
      if(bcheck == 2) return 1; //The particle was lost physically during the bounce.
      if(bcheck == -1) return -1; //The particle was lost unphysically during the bounce.
      printf("+ERROR: Eroneous code received from 'bounce'! [Particle# = %d\n]",neutron.num);
      return -1;

      if(CONSOLE == "ON") printf("+WARNING: Error during particle %d's encounter with a bulk medium!\n",neutron.num);
      return -1; //Drop particle since error in determination of particle interaction with medium.

	}
	}        
  printf("+ERROR: 'cplanehandling' called with no result. Particle number %d abandoned!\n",neutron.num);
  return -1;
}

void recdet(int detectors[5000][10], int *counts) {
  int i;
    
  fprintf(detectorsfp,"Total Counts:\n");
  for(i=0 ; i < 10 ; i++) fprintf(detectorsfp,"Detector %d -> %d\n",(i+1),counts[i]);
  
  fprintf(detectorsfp,"\n\nTime      |   Counts D1   |   Counts D2   |   Counts D3   |   Counts D4   |   Counts D5   |   Counts D6   |   Counts D7   |   Counts D8   |   Counts D9   |   Counts D10 \n");
  for(i=0 ; i < 5000 ; i++) {
    fprintf(detectorsfp,"%lf ",i/10.);
    fprintf(detectorsfp,"      %d               %d               %d",detectors[i][0],detectors[i][1],detectors[i][2]);
    fprintf(detectorsfp,"               %d               %d",detectors[i][3],detectors[i][4]);
    fprintf(detectorsfp,"               %d               %d               %d",detectors[i][5],detectors[i][6],detectors[i][7]);
    fprintf(detectorsfp,"               %d               %d\n",detectors[i][8],detectors[i][9]);
  }
  return;
}

int paramadjust(char action[5],char pname[10],int startreg,int endreg,double pvalue) {
  int i;
  
   //Note that 'strcmp' is a function defined in the <string.h> header which compares two strings. Format is: int strcmp(char *s1,cahr *s2).
  //The function returns -1 if s1 < s2. 0 if s1 = s2, and 1 if s1 > s2.
  if(strcmp(pname,"gradBx") == 0) {  //'gradBx' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'gradBx' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][0]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][0] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'gradBx' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][0] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"gradBy") == 0) {  //'gradBy' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'gradBy' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][1]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][1] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'gradBy' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][1] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"gradBz") == 0) {  //'gradBz' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'gradBz' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][2]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][2] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'gradBz' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][2] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"spec") == 0) {  //'spec' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'spec' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][3]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][3] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'spec' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][3] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"loss") == 0) {  //'loss' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'loss' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][4]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][4] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'loss' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][4] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"depol") == 0) {  //'depol' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'depol' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][5]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][5] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'depol' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][5] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"wpot") == 0) {  //'wpot' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'wpot' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][6]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][6] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'wpot' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][6] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"bpot") == 0) {  //'bpot' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'bpot' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][7]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][7] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'bpot' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][7] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"scat") == 0) {  //'scat' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'scat' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][8]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][8] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'scat' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][8] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"abs") == 0) {  //'abs' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'abs' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][9]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][9] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'abs' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][9] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"PM") == 0) {  //'PM' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'PM' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][11]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][11] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'PM' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][11] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"SM") == 0) {  //'SM' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'SM' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][12]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][12] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'SM' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][12] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"LM") == 0) {  //'LM' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'LM' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][13]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][13] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'LM' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][13] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  if(strcmp(pname,"DM") == 0) {  //'DM' to be adjusted.
    if(strcmp(action,"set") == 0) {  //'DM' is to be set to a new value.
      for(i=startreg ; i <= endreg ; i++) bparamhold[i] = rparams[i][14]; //Save the old set of parameters.
      for(i=startreg ; i <= endreg ; i++) rparams[i][14] = pvalue; //Load the new set of parameters.
      return 0; //Parameter name recognized and adjusted... return OK value.
    }
    if(strcmp(action,"reset") == 0)  {  //'DM' is to be reset to original value.
      for(i=startreg ; i <= endreg ; i++) rparams[i][14] = bparamhold[i]; //Reset to the original set of parameters.
      for(i=0 ; i < 100 ; i++) bparamhold[i] = 0; //Zero the array for safety.
      return 0; //Parameter name recognized and reset... return OK value.
    }
    printf("+ERROR: Unrecognized action code '%s' in 'paramadjust'!\n",action);
    return -1; //Return the error code.
  }
  printf("+ERROR: 'batch' file calls for a change in the unrecognized parameter '%s'.\n",pname);
  return -1; //Return the error code.
}

void event(int eventcode,double data1,double data2,double data3, double data4) {
  
  fprintf(eventsfp,"%-11d %-7d %-6d %-9d",neutron.num,eventcode,neutron.region,neutron.xcode);
  fprintf(eventsfp,"%-10f %-11f %-11f %-12f",neutron.t,neutron.x,neutron.y,neutron.z);
  fprintf(eventsfp,"%-11f %-12f %-12f",neutron.vx,neutron.vy,neutron.vz);
  fprintf(eventsfp,"%-11f %-11f %-11f %-11f %-12f\n",neutron.spin[2],data1,data2,data3,data4);
  return;
}

int perturbednormal(double *nx,double *ny,double *nz){ //created by Erik Lutz
  double inx,iny,inz;
  double surfparamexp, surfparamlin;
  double theta, phi, randthetaseed, randtheta, randphiseed, randphi;
  double thetaf, phif;
  int modeltype;
    
  normal(&inx,&iny,&inz,0,0,0);
  
  randthetaseed = grn();
  randphiseed = grn();
  
  theta = acos(inx);   //define current inclination angle
  phi = atan2(inz,iny); //define current inclination angle
    
  if ((cplanes[neutron.xcode][3] <= 1.) || ((cplanes[neutron.xcode][3] > 2.) && (cplanes[neutron.xcode][3] <= 3.))) { //use exponential model
  
    if (cplanes[neutron.xcode][3] <= 1.) {surfparamexp = 1.0/cplanes[neutron.xcode][3];}    //grab user supplied surface roughness parameter for exponential model}
    else {surfparamexp = 1.0/(cplanes[neutron.xcode][3]-2);}
    
    randtheta = (1./2.)*acos(1.-2.*pow(randthetaseed,surfparamexp));         //generate a weighted random inclination angle seed between 0 and Pi/2
    randphi = (PI/2.)*pow(randphiseed,surfparamexp);          //generate an weighted random azimuthal angle between 0 and Pi/2
  }
  
  if (((cplanes[neutron.xcode][3] > 1.) && (cplanes[neutron.xcode][3] <= 2.)) || ((cplanes[neutron.xcode][3] > 3.) && (cplanes[neutron.xcode][3] <= 4.))) { //use linear model
  
    if ((cplanes[neutron.xcode][3] > 1.) && (cplanes[neutron.xcode][3] <= 2.)) {surfparamlin = cplanes[neutron.xcode][3]-1.;}    //grab user supplied surface roughness parameter for linear model
    else {surfparamlin = cplanes[neutron.xcode][3]-3.;}
  
    randtheta = (1./2.)*acos(1.-2.*surfparamlin*randthetaseed);
    randphi = (PI/2.)*surfparamlin*randphiseed;
  }
  
  if(grn() < 0.5) randtheta = -randtheta;   //give random theta a random sign, now is a weighted random number in [-Pi/2,Pi/2]
  if(grn() < 0.5) randphi = -randphi;   //give random phi a random sign, now is a weighted random number in [-Pi/2,Pi/2]
  
  thetaf = theta + 2*randtheta;
  phif = phi + 2*randphi;
  
  *nx = cos(thetaf);
  *ny = sin(thetaf)*cos(phif);
  *nz = sin(thetaf)*sin(phif);
  
  if ((cplanes[neutron.xcode][3] > 0.) && (cplanes[neutron.xcode][3] <= 1.)) { return 0;} //exponential model, no pass through
  else if ((cplanes[neutron.xcode][3] > 1.) && (cplanes[neutron.xcode][3] <= 2.)) { return 1;} //linear model, no pass through
  else if ((cplanes[neutron.xcode][3] > 2.) && (cplanes[neutron.xcode][3] <= 3.)) { return 2;} //exponential model, pass through
  else if ((cplanes[neutron.xcode][3] > 3.) && (cplanes[neutron.xcode][3] <= 4.)) { return 3;} // linear model, pass through
  
  return -1;
}
