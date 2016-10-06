#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
////#include <gsl/gsl_rng.h>
////#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "T2toolkit.h"
//#include "cpgplot.h"
//#include <omp.h>
#define WINSIZE 2.5
#define PRECISION 0.005

typedef struct acfStruct {
	int n; // number of dynamic spectrum

	double cFlux;  // pulsar flux density
	double whiteLevel;  // white noise level

	double cFreq; // observing central frequency
	double bw; // observing bandwidth
	double f0; // scintillation bandwidth
	double tint; // integration time
	double t0; // scintillation time-scale
	int nchn;
	int nsubint;

	double **dynSpec; // dynamic spectrum 

	float *dynPlot; // dynamic spectrum for pgplot

	float probability;
} acfStruct;

typedef struct noiseStruct {
	int n; // number of dynamic spectrum
	int nchn;
	int nsubint;

	float *noisePlot; // dynamic spectrum for pgplot
	double whiteLevel;  // white noise level
	double detection; // detection threshold
} noiseStruct;

typedef struct controlStruct {
	int n; // number of dynamic spectrum
	//char oname[1024]; // output file name
	char dname[1024]; // graphics device

	double T; // total integration time
	double BW; // total BW
	double cFreq;  // observing central frequency
	double tsub;  // subintegration time
	int nsub;  // number of subintegrations
	double chanBW;  // subchannel bandwidth
	int nchan; // number of subchannels

	double scint_freqbw;   // scintillation bandwidth
	double scint_ts;       // scintillation timescale

	double whiteLevel;   // white noise level, mJy
	double whiteLevel0;   // white noise level, mJy
	double cFlux;        // flux density of pulsars, mJy
	double cFlux0;        // flux density of pulsars, mJy
	double cFlux1;        // flux density of pulsars, mJy
	//double fluxStep;      

	int noplot;

	double precision;
}controlStruct;

void deallocateMemory (acfStruct *acfStructure);
void allocateMemory (acfStruct *acfStructure);
int simDynSpec (acfStruct *acfStructure, long seed);
int winDynSpec (acfStruct *acfStructure, long seed);
int calculateScintScale (acfStruct *acfStructure, controlStruct *control);
int calculateNDynSpec (acfStruct *acfStructure, controlStruct *control, noiseStruct *noiseStructure);
float find_peak_value (int n, float *s);
int readParams(char *fname, char *dname, int n, controlStruct *control);
//int readParams(char *fname, char *oname, char *dname, int n, controlStruct *control);
void initialiseControl(controlStruct *control);

//void heatMap (acfStruct *acfStructure, char *dname);
//void palett(int TYPE, float CONTRA, float BRIGHT);
//int plotDynSpec (char *pname, char *dname);

//int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control);
float chiSquare (float *data, int n, float noise);
float moduIndex (float *data, int n);
float variance (float *data, int n);
int histogram (float *data, int n, float *x, float *val, float low, float up, int step);
float find_max_value (int n, float *s);
float find_min_value (int n, float *s);

int calNoise (noiseStruct *noiseStructure, controlStruct *control);
int simNoise (noiseStruct *noiseStructure, long seed);
int calThreshold (noiseStruct *noiseStructure, float *var_n);

void readDiss (char *Tname, double *tdiss);
int readDissNum (char *Tname);

int calNoise (noiseStruct *noiseStructure, controlStruct *control)
{
	int nchn, nsubint;
	long seed;
	int i;
	int n = control->n;

	float *var_n;
	var_n = (float*)malloc(sizeof(float)*n);

	noiseStructure->n = control->n; 
	noiseStructure->whiteLevel = control->whiteLevel; // mJy
	noiseStructure->nchn = control->nchan; 
	noiseStructure->nsubint = control->nsub; 

	nchn = noiseStructure->nchn;
	nsubint = noiseStructure->nsubint;

	// allocate memory
	noiseStructure->noisePlot = (float *)malloc(sizeof(float)*nsubint*nchn);

	// simulate noise
	for (i=0; i<noiseStructure->n; i++)
	{
		seed = TKsetSeed();
		simNoise (noiseStructure, seed);

		var_n[i] = variance (noiseStructure->noisePlot, nsubint*nchn);
	}

	calThreshold (noiseStructure, var_n);

	// deallocate memory
	free(noiseStructure->noisePlot);
	free(var_n);

	return 0;
}

int calThreshold (noiseStruct *noiseStructure, float *var_n)
{
	int i, nMax;
	int n = noiseStructure->n;

	float max, min;
	float percent;
	double threshold;
	int num;

	max = find_max_value (n, var_n);
	min = find_min_value (n, var_n);
	//printf ("Results: %f %f %f %f\n", meanVar, varVar, meanVarN, varVar_n);
	
	percent = 1.0;
	nMax = 0;
	//threshold = max;
	//while (percent >= 0.95)
	while (fabs(percent-0.95) >= 0.001 && nMax <= 100)
	{
		//threshold = threshold - 0.01;
		threshold = min+(max-min)/2.0;

		num = 0;
		for (i=0; i<n; i++)
		{
			if (var_n[i] <= threshold)
			{
				num++;
			}
		}
		percent = (float)(num)/n;
		
		if (percent > 0.95)
		{
			max = threshold;
		}
		else
		{
			min = threshold;
		}
		nMax++;

		//printf ("%.6lf %.6lf %d\n", threshold, percent, num);
	}

	noiseStructure->detection = threshold;

	return 0;
}

int calculateNDynSpec (acfStruct *acfStructure, controlStruct *control, noiseStruct *noiseStructure)
{
	long seed;
	int i;

	int n = acfStructure->n;
	int nsub = control->nsub;
	int nchan = control->nchan;

	int num;
	
	float var;
	//float *var;
	//var = (float*)malloc(sizeof(float)*n);

	acfStructure->cFlux = control->cFlux; // mJy

	num = 0;
	for (i=0; i<acfStructure->n; i++)
	{
		seed = TKsetSeed();
	
		simDynSpec (acfStructure, seed);

		winDynSpec (acfStructure, seed);
		//printf ("Make DynSpec %d\n", i);
		
		// calculate the variance
		var = variance (acfStructure->dynPlot, nsub*nchan);
		//var[i] = variance (acfStructure->dynPlot, nsub*nchan);

		//printf ("Var %f %f\n", var, noiseStructure->detection);
		if (var >= noiseStructure->detection)
		//if (var[i] >= noiseStructure->detection)
		{
			num++;
		}
	}

	//printf ("%d\n", num);
	acfStructure->probability = (float)(num)/n;

	//free(var);

	return 0;
}

int calculateScintScale (acfStruct *acfStructure, controlStruct *control)
//int calculateScintScale (acfStruct *acfStructure, controlStruct *control, long seed)
{
	//FILE *fin;
	int nchn, nsubint;

	//printf ("Starting simulating dynamic spectrum\n");
	// moved to preAllocateMemory
	acfStructure->n = control->n; 
	acfStructure->whiteLevel = control->whiteLevel; // mJy

	acfStructure->cFreq = control->cFreq; // MHz
	acfStructure->bw = fabs(control->chanBW*control->nchan); // MHz
	acfStructure->tint = control->nsub*control->tsub;  // s

	acfStructure->f0 = control->scint_freqbw;  // MHz
	acfStructure->t0 = control->scint_ts; // s

	acfStructure->nchn = control->nchan;
	acfStructure->nsubint = control->nsub;

	nchn = acfStructure->nchn;
	nsubint = acfStructure->nsubint;
	//printf ("Scintillation bandwidth: %lf (MHz)\n", acfStructure->f0);
	//printf ("Scintillation time-scale: %lf (s)\n", acfStructure->t0);

	allocateMemory (acfStructure);

	//for (i=0; i<acfStructure->n; i++)
	//{
	//	seed = TKsetSeed();
	//	//acfStructure->phaseGradient = TKgaussDev(&seed);
	//	acfStructure->phaseGradient = 0.0;
	//	//printf ("Phase gradient: %lf\n", acfStructure->phaseGradient);

	//	winDynSpec (acfStructure, seed, i);
	//}

	/*
	if (acfStructure->n == 1)
	{
		printf ("Dynamic spectrum is output into %s\n", control->oname);

		if ((fin=fopen(control->oname,"w"))==NULL)
		{
			printf ("Can't open output file!\n");
			exit(1);
		}

		fprintf(fin,"INFO nsub nchn bandwidth tint cFreq\n");
		fprintf(fin,"START %d %d %lf %lf %lf\n",acfStructure->nsubint,acfStructure->nchn,acfStructure->bw,acfStructure->tint,acfStructure->cFreq);

		for (i=0;i<acfStructure->nchn;i++)
		{
			for (j=0;j<acfStructure->nsubint;j++)
			{
				fprintf(fin,"%d %d %lf\n", i, j, acfStructure->dynPlot[0][i*acfStructure->nsubint+j]);
			}
		}

		if (fclose(fin))
		{
			printf ("Can't close output file!\n");
			exit(1);
		}
	}
	else
	{
		printf ("%d dynamic spectra are simulated\n", acfStructure->n);
	}
	*/

	return 0;
}

void allocateMemory (acfStruct *acfStructure)
{
	int i;
	int nsub = acfStructure->nsubint;
	int nchn = acfStructure->nchn;
	
	acfStructure->dynSpec = (double **)malloc(sizeof(double *)*nchn);

	for (i = 0; i < nchn; i++)
	{
		acfStructure->dynSpec[i] = (double *)malloc(sizeof(double)*nsub);
	}
	
	acfStructure->dynPlot = (float *)malloc(sizeof(float)*nsub*nchn);
}

void deallocateMemory (acfStruct *acfStructure)
{
	//int n = acfStructure->n; // number of dynamic spectrum
	int nchn = acfStructure->nchn;
	//int nchn = acfStructure->nchn;
	
	int i;
	for (i = 0; i < nchn; i++)
	{
		free(acfStructure->dynSpec[i]);
	}

	//for (i = 0; i < nchn; i++)
	//{
	//	free(acfStructure->dynSpecWindow[i]);
	//}

	free(acfStructure->dynPlot);
}

int simNoise (noiseStruct *noiseStructure, long seed)
{
	int nchn = noiseStructure->nchn;
	int nsubint = noiseStructure->nsubint;

	int i, j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsubint; j++)
		{
			noiseStructure->noisePlot[i*nsubint+j] = (float)(noiseStructure->whiteLevel*TKgaussDev(&seed));   // create noise image pixels
		}
	}

	return 0;
}

int simDynSpec (acfStruct *acfStructure, long seed)
{
	int nchn = acfStructure->nchn;
	int nsub = acfStructure->nsubint;

	int i;
	int j;

	int n = 0;
	double sum = 0.0;

	//printf ("Simulating dynamic spectrum\n");
	//seed = TKsetSeed();
	//printf ("seed %ld\n",seed);

	/////////////////////////////////////////////////////////////////////////////////

	// form the matrix and normalize
	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsub; j++)
		{
			acfStructure->dynSpec[i][j] = pow(TKgaussDev(&seed),2.0)+pow(TKgaussDev(&seed),2.0);
			//fprintf (fp, "%lf  ", acfStructure->dynSpec[i][j]);
			sum += acfStructure->dynSpec[i][j];
			n++;
		}
		//fprintf (fp, "\n");
	}

	sum = sum/n;

	//printf ("Normalization %.10lf\n",sum);
	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsub; j++)
		{
			acfStructure->dynSpec[i][j] = acfStructure->dynSpec[i][j]/sum;
			//printf ("%d %d\n", i, j);
		}
	}

	return 0;
}

int winDynSpec (acfStruct *acfStructure, long seed)
{
	int nchn = acfStructure->nchn;
	int nsubint = acfStructure->nsubint;

	int i;
	int j;

	double dynSpecWindow;

	//printf ("Simulating dynamic spectrum\n");
	//seed = TKsetSeed();
	//printf ("seed %ld\n",seed);

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsubint; j++)
		{
			dynSpecWindow = acfStructure->dynSpec[i][j];
			//acfStructure->dynSpecWindow[i][j] = acfStructure->dynSpec[i+nf0][j+ns0];
			acfStructure->dynPlot[i*nsubint+j] = (float)(dynSpecWindow*acfStructure->cFlux+acfStructure->whiteLevel*TKgaussDev(&seed));   // add in noise here
			//acfStructure->dynPlot[i*nsubint+j] = (float)(acfStructure->dynSpecWindow[i][j]*acfStructure->cFlux+acfStructure->whiteLevel*TKgaussDev(&seed));   // add in noise here
			//printf ("noise rand %lf\n",TKgaussDev(&seed));
			//acfStructure->dynPlot[i*nsubint+j] = (float)(acfStructure->dynSpecWindow[i][j]);
			//fprintf (fp, "%.10lf  ", acfStructure->dynSpec[i][j]/sum);
		}
	}

	return 0;
}

float find_peak_value (int n, float *s)
{
	int i;
	float temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	return temp[n-1];
}
						
int readParams(char *fname, char *dname, int n, controlStruct *control)
//int readParams(char *fname, char *oname, char *dname, int n, controlStruct *control)
{
	FILE *fin;
	char param[1024];
	int endit=-1;
	int finished=0;

	// define the output file name
	//strcpy(control->oname,oname);
	strcpy(control->dname,dname);

	control->n = n;

	///////////////////////////////////////
	if ((fin=fopen(fname,"r"))==NULL)
	{
		printf ("Can't open file!\n");
		exit(1);
	}

	//printf("Reading parameters...\n");

      	// Find the start observation
	while (!feof(fin))
	{
		if (fscanf(fin,"%s",param)==1)
		{
			if (strcasecmp(param,"START_OBS")==0)
			{
				endit=0;
				break;
			}
		}
		else 
			return 1;
	}

	if (endit==-1)
		return 1;

	do
	{
		fscanf(fin,"%s",param);
		if (strcasecmp(param,"END_OBS")==0)
			endit=1;
		else
		{
			//if (strcasecmp(param,"SCINT_TS")==0)
			//	fscanf(fin,"%lf",&(control->scint_ts));
			//else if (strcasecmp(param,"SCINT_FREQBW")==0)
			//	fscanf(fin,"%lf",&(control->scint_freqbw));	  
			if (strcasecmp(param,"T")==0)
      				fscanf(fin,"%lf",&(control->T));
			else if (strcasecmp(param,"BW")==0)
      				fscanf(fin,"%lf",&(control->BW));
			else if (strcasecmp(param,"CFREQ")==0)
      				fscanf(fin,"%lf",&(control->cFreq));
			//else if (strcasecmp(param,"CHAN_BW")==0)
      			//	fscanf(fin,"%lf",&(control->chanBW));
			//else if (strcasecmp(param,"NCHAN")==0)
			//     	fscanf(fin,"%d",&(control->nchan));
			//else if (strcasecmp(param,"NSUB")==0)
      			//	fscanf(fin,"%d",&(control->nsub));
			else if (strcasecmp(param,"WHITE_LEVEL")==0)    // noise level of the survey
      				fscanf(fin,"%lf",&(control->whiteLevel0));
			else if (strcasecmp(param,"CFLUX0")==0)
			      	fscanf(fin,"%lf",&(control->cFlux0));
			else if (strcasecmp(param,"CFLUX1")==0)
			      	fscanf(fin,"%lf",&(control->cFlux1));
			else if (strcasecmp(param,"PRECISION")==0)
			      	fscanf(fin,"%lf",&(control->precision));
		}
	} while (endit==0);

	//////////////////////////////////////////////////////////////////////
	//if (control->tsys != 0.0 && control->tsky != 0.0 && control->gain != 0.0 && control->whiteLevel == 0)
	//{
	//	control->radioNoise = (control->tsys+control->tsky)/(control->gain)/sqrt(2.0*(control->tsubRequested/control->nbin)*(fabs(control->obsBW)/control->nchan));
	//}
	//else if (control->tsys == 0.0 && control->tsky == 0.0 && control->gain == 0.0 && control->whiteLevel != 0)
	//{
	//	control->radioNoise = control->whiteLevel;
	//}
	//else 
	//{
	//	printf ("Double definiation of radio-meter noise!\n");
	//	exit (1);
	//}
	//printf ("Nchan: %d; Tsys: %lf; Tsky: %lf; Gain: %lf; Radio-meter noise: %lf mJy\n", control->nchan, control->tsys, control->tsky, control->gain, control->radioNoise);

	if (fclose(fin))
	{
		printf ("Can't close file.\n");
		exit(1);
	}

	return finished;
}

void initialiseControl(controlStruct *control)
{
	//strcpy(control->primaryHeaderParams,"UNKNOWN");
	//strcpy(control->exact_ephemeris,"UNKNOWN");
	//strcpy(control->src,"UNKNOWN");
	//strcpy(control->oname,"UNKNOWN");
	strcpy(control->dname,"1/xs");
	
	control->n = 1; // simulate 1 dynamic spectrum by default

	// Standard defaults
	control->T = 1000.0;
	control->BW = 1000.0;
	control->nchan = 10;
	control->nsub = 10;
	control->cFreq = 1400.0;
	control->chanBW = 100;   // MHz
	control->tsub = 100;  // second
	control->whiteLevel = 0;   // mJy
	control->whiteLevel0 = 0;   // mJy
	control->scint_ts  = 50.0;  // second
	control->scint_freqbw = 50.0;   // MHz
	//control->scint_freqbw0 = 0.0;   // MHz
	//control->scint_freqbw1 = 0.0;   // MHz
	//control->scint_f_step = 0.0;   // MHz
	control->cFlux = 0.0;   // mJy
	control->cFlux0 = 0.0;   // mJy
	control->cFlux1 = 0.0;   // mJy

	control->precision = 0.05;   
}

//void heatMap (acfStruct *acfStructure, char *dname)
//{
//	//int i,j;                     
//	//int dimx = acfStructure.ns;
//	//int dimy = acfStructure.nf; // dimensions 
//	//float tab[dimx*dimy];       // value
//	char caption[1024];
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", acfStructure->cFreq, "MHz", "BW:", acfStructure->bw, "MHz", "Length:", acfStructure->tint, "s");
//  
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	int dimx = acfStructure->nsubint;
//	int dimy = acfStructure->nchn;
//	double bw = acfStructure->bw;
//  
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//
//	double f1 = acfStructure->cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = acfStructure->cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//printf ("f1 f2: %lf %lf\n", f1, f2);
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, acfStructure->dynPlot[0]);
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
//
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	//cpgswin(0,dimx,164,132); // set window
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//	//cpgbox("BCTSIN",10,5,"BCTSIN",50,5);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//      	//cpglab("Subintegration","Frequency (MHz)","Freq: 150.0 MHz BW: -32.000 MHz Length: 960.0 s");
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(acfStructure->dynPlot[0],dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(acfStructure->dynPlot,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//} 
//
//void palett(int TYPE, float CONTRA, float BRIGHT)
//{
////-----------------------------------------------------------------------
//// Set a "palette" of colors in the range of color indices used by
//// PGIMAG.
////-----------------------------------------------------------------------
//	float GL[] = {0.0, 1.0};
//	float GR[] = {0.0, 1.0};
//	float GG[] = {0.0, 1.0};
//	float GB[] = {0.0, 1.0};
//	float RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
//	float RR[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
//	float RG[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
//	float RB[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
//	float HL[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float HR[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float HG[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float HB[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	float WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
//	float WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
//	float WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
//	float WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
//	float AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
//	float AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//	float AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
//	float AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//      
//	if (TYPE == 1)
//	{   
//		//-- gray scale
//		cpgctab(GL, GR, GG, GB, 2, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 2) 
//	{
//		//-- rainbow
//		cpgctab(RL, RR, RG, RB, 9, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 3) 
//	{
//		//-- heat
//		cpgctab(HL, HR, HG, HB, 5, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 4) 
//	{
//		//-- weird IRAF
//		cpgctab(WL, WR, WG, WB, 10, CONTRA, BRIGHT);
//	}
//	else if (TYPE == 5) 
//	{
//		//-- AIPS
//		cpgctab(AL, AR, AG, AB, 20, CONTRA, BRIGHT);
//	}
//}

//int plotDynSpec (char *pname, char *dname)
//{
//	FILE *fin;
//	char start[128];
//	double tint,bw,cFreq;
//	int nsub,nchn;
//	float val;
//	int n1, n2;
//	float *dynSpec;
//
//	int i;
//
//	int dimx;
//	int dimy;
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	char caption[1024];
//
//	if ((fin=fopen(pname,"r"))==NULL)
//	{
//		printf ("Can't open dynamic spectrum!\n");
//		exit(1);
//	}
//
//      	// Find the start of dynamic spectrum, which contains basic info
//	while (!feof(fin))
//	{
//		if (fscanf(fin,"%s %d %d %lf %lf %lf",start,&nsub,&nchn,&bw,&tint,&cFreq)==6)
//		{
//			if (strcasecmp(start,"START")==0)
//			{
//				break;
//			}
//		}
//		//else 
//		//	return 1;
//	}
//
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", cFreq, "MHz", "BW:", bw, "MHz", "Length:", tint, "s");
//
//	dynSpec = (float*)malloc(sizeof(float)*nsub*nchn);
//
//	i = 0;
//	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
//	{
//		dynSpec[i] = val;
//		i++;
//	}
//
//	if (fclose(fin))
//	{
//		printf ("Can't close dynamic spectrum!\n");
//		exit(1);
//	}
//
//	dimx = nsub;
//	dimy = nchn;
//  
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, dynSpec);
//
//	double f1 = cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
// 
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,dname,1,1);
//	//cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//
//	free(dynSpec);
//
//	return 0;
//} 

//int qualifyVar (acfStruct *acfStructure, noiseStruct *noiseStructure, controlStruct *control)
//{
//	int i;
//	int n = acfStructure->n;
//	int n_n = noiseStructure->n;
//
//	float m0;
//	//float *m;
//
//	float *var; 
//	//float *var_n;
//	//float varVar, meanVar;
//	//float varVar_n, meanVar_n;
//
//	int nsub = control->nsub;
//	int nchan = control->nchan;
//
//	int num;
//
//	//m = (float*)malloc(sizeof(float)*n);
//	var = (float*)malloc(sizeof(float)*n);
//	//var_n = (float*)malloc(sizeof(float)*n_n);
//
//	for (i=0; i<n; i++)
//	{
//		//m[i] = moduIndex (acfStructure->dynPlot[i], nsub*nchan);
//		var[i] = variance (acfStructure->dynPlot[i], nsub*nchan);
//		//printf ("%f\n", var[i]);
//		//printf ("%f \n", m[i]);
//	}
//
//	num = 0;
//	for (i=0; i<n; i++)
//	{
//		if (var[i] >= noiseStructure->detection)
//		{
//			num++;
//		}
//	}
//	acfStructure->probability = (float)(num)/n;
//	//printf ("%d %d %f %f\n", num, n, (float)(num)/n, acfStructure->probability);
//
//	//for (i=0; i<n_n; i++)
//	//{
//	//	var_n[i] = variance (noiseStructure->noisePlot[i], nsub*nchan);
//	//	//fprintf (fin2, "%f\n", var_n[i]);
//	//}
//
//	/*
//	if (fclose(fin1))
//	{
//		printf ("Can't close file...\n");
//		exit(1);
//	}
//
//	if (fclose(fin2))
//	{
//		printf ("Can't close file...\n");
//		exit(1);
//	}
//	*/
//
//	//m0 = 0.0;
//	//meanVar = 0.0;
//	//for (i=0; i<n; i++)
//	//{
//	//	m0 += m[i];
//	//	meanVar += var[i];
//	//}
//	//m0 = m0/n;
//	//meanVar = meanVar/n;
//
//	//meanVar_n = 0.0;
//	//for (i=0; i<n_n; i++)
//	//{
//	//	meanVar_n += var_n[i];
//	//}
//	//meanVar_n = meanVar_n/n_n;
//
//	//varVar = variance (var, n);
//	//varVar_n = variance (var_n, n_n);
//
//	//printf ("Results: %f %f %f %f %f\n", m0, meanVar, varVar, meanVar_n, varVar_n);
//
//	/*
//	maxV = find_max_value(n,var);
//	minV = find_min_value(n,var);
//	maxVN = find_max_value(n_n,var_n);
//	minVN = find_min_value(n_n,var_n);
//	////////////////////////////
//	// make histogram
//	if (control->noplot != 1)
//	{
//		cpgbeg(0,control->dname,1,1);
//
//		xHis = (float*)malloc(sizeof(float)*step);
//		val = (float*)malloc(sizeof(float)*step);
//		xHisN = (float*)malloc(sizeof(float)*step);
//		valN = (float*)malloc(sizeof(float)*step);
//		histogram (var, n, xHis, val, minV-2*(maxV-minV)/step, maxV+2*(maxV-minV)/step, step);
//		histogram (var_n, n_n, xHisN, valN, minVN-2*(maxVN-minVN)/step, maxVN+2*(maxVN-minVN)/step, step);
//
//      		cpgsch(1); // set character height
//      		cpgscf(1); // set character font
//
//		// find the max and min
//		max1 = find_max_value(step,val);
//		max2 = find_max_value(step,valN);
//		max = (max1 >= max2 ? max1 : max2);
//      		//cpgenv(-5,5,0,4500,0,1); // set window and viewport and draw labeled frame
//      		cpgenv(minVN-1, maxV+1, 0, max+0.1*max, 0, 1); // set window and viewport and draw labeled frame
//
//		sprintf(caption, "%s", "Variance histogram");
//      		cpglab("Variance","Number",caption);
//		cpgbin(step,xHis,val,0);
//		cpgsci(2);
//		cpgbin(step,xHisN,valN,0);
//		cpgsci(3);
//		cpgline(100, xThreshold, yThreshold);
//		///////////////////////////////////////////////////////
//		cpgend();
//
//		free(xHis);
//		free(val);
//		free(xHisN);
//		free(valN);
//	}
//	*/
//
//	//free(m);
//	free(var);
//	//free(var_n);
//	//free(flux0);
//	//free(flux);
//
//	return 0;
//}

float chiSquare (float *data, int n, float noise)
{
	int i;

	float ave;
	float chiS;
	
	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	chiS = 0.0;
	for (i=0; i<n; i++)
	{
		chiS += pow(data[i]-ave,2)/pow(noise,2);
	}

	return chiS/(n-1);
}

float moduIndex (float *data, int n)
{
	int i;

	float ave, devi;
	float m;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = sqrt(devi/n);

	m = devi/ave;

	return m;
}

float variance (float *data, int n)
{
	int i;

	float ave, devi;

	ave = 0.0;
	for (i=0; i<n; i++)
	{
		ave += data[i];
	}
	ave = ave/n;

	devi = 0.0;
	for (i=0; i<n; i++)
	{
		devi += pow(data[i]-ave,2);
	}
	devi = devi/n;

	return devi;
}

int histogram (float *data, int n, float *x, float *val, float low, float up, int step)
{
	int i,j,count;
	float width;
	float *temp;

	temp = (float*)malloc(sizeof(float)*(step+1));

	width = (up-low)/step;
	for (i=0; i<step; i++)
	{
		x[i] = low + i*width + width/2.0;
	}

	for (i=0; i<=step; i++)
	{
		temp[i] = low + i*width;
	}

	for (i=0; i<step; i++)
	{
		count = 0;
		for (j=0; j<n; j++)
		{
			if (data[j]>=temp[i] && data[j]<temp[i+1])
			{
				count += 1;
			}
		}
		//val [i] = count;
		val[i] = (float)(count)/n;
	}

	free(temp);
	return 0;
}

float find_max_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

float find_min_value (int n, float *s)
{
	int i;
	float *temp;

	temp = (float *)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a <= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	c = temp[n-1];
	free(temp);

	return c;
}

int readDissNum (char *Tname)
{
	int nt;
	FILE *fint;
	double temp;

	if ((fint=fopen(Tname, "r"))==NULL)
	{
		printf ("Can't open scintillation time-scale...\n");
		exit(1);
	}

	nt = 0;
	while (fscanf(fint, "%lf", &temp) == 1)
	{
		nt++;
	}

	if (fclose(fint))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}

	return nt;
}



void readDiss (char *Tname, double *tdiss)
{
	int i;
	FILE *fint;
	double temp;

	///////////////////////////////////////////////////////////////////////
	if ((fint=fopen(Tname, "r"))==NULL)
	{
		printf ("Can't open scintillation time-scale...\n");
		exit(1);
	}

	i = 0;
	while (fscanf(fint, "%lf", &temp) == 1)
	{
		tdiss[i] = temp;
		i++;
	}

	if (fclose(fint))
	{
		printf ("Can't close scintillation time-scale...\n");
		exit(1);
	}
}
