// Sensitivity map of matched filter as a function of scintillation time-scale and bandwidth for a given survey
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "detection.h"
//#include <mpi.h>

int main (int argc, char* argv[])
{
	//int id;  //  process rank
	//int p;   //  number of processes
	//int num;
	long seed;
	int nMax;

	//MPI_Init (&argc, &argv);
	//MPI_Comm_rank (MPI_COMM_WORLD, &id);
	//MPI_Comm_size (MPI_COMM_WORLD, &p);

	//FILE *fin;
	controlStruct control;
	acfStruct acfStructure;
	noiseStruct noiseStructure;

	char fname[1024];   // read in parameter file
	char Tname[1024];   // read in scintillation time-scales
	char Fname[1024];   // read in scintillation bandwidth

	int nt;
	double *tdiss;
	int nf;
	double *fdiss;

	//char oname[1024];   // output file name
	//char pname[1024];   // file to plot
	char dname[1024];   // graphics device

	int i, j;
	int n = 1;  // number of dynamic spectrum to simulate

	double flux0;
	double flux1;

	int plotMode = 0;  // plot existing dynamic spectrum, off by default
	control.noplot = 1;  // don't show dynamic spectrum while simulationg, on by default

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
			strcpy(Tname,argv[++i]);
			strcpy(Fname,argv[++i]);
			printf ("Parameters are in %s\n", fname);
		}
		//else if (strcmp(argv[i],"-o")==0)
		//{
		//	strcpy(oname,argv[++i]);
		//	//printf ("Dynamic spectrum is output into %s\n", oname);
		//}
		else if (strcmp(argv[i],"-n")==0)
		{
			n = atoi (argv[++i]);
			printf ("Number of dynamic spectrum to simulate %d\n", n);
		}
		//else if (strcmp(argv[i],"-p")==0) // just plot 
		//{
		//	strcpy(pname,argv[++i]);
		//	plotMode = 1;
		//	printf ("Plotting exiting dynamic spectrum.\n");
		//}
		//else if (strcmp(argv[i],"-noplot")==0)
		//{
		//	control.noplot = 1; // Don't show dynamic spetrum while simulationg
		//}
		//else if (strcmp(argv[i],"-dev")==0)
		//{
		//	strcpy(dname,argv[++i]);
		//	printf ("Graphic device: %s\n", dname);
		//}
	}

	nt = readDissNum (Tname);
	nf = readDissNum (Fname);
	tdiss = (double *)malloc(sizeof(double)*nt);
	fdiss = (double *)malloc(sizeof(double)*nf);
	readDiss (Tname, tdiss);
	readDiss (Fname, fdiss);

	sprintf(dname, "%s", "1/xs");
	/*
	if ((fin=fopen(oname, "w"))==NULL)
	{
		printf ("Can't open file...\n");
		exit(1);
	}
	*/

	if (plotMode==0)
	{
		// Simulate dynamic spectrum
		// read parameters
		initialiseControl(&control);
		readParams (fname,dname,n,&control);
		//readParams (fname,oname,dname,n,&control);
		printf ("Finished reading parameters.\n");

		//preAllocateMemory (&acfStructure, &control);
		//allocateMemory (&acfStructure);

		//for (i=id; i<nt; i+=p)
		for (i=0; i<nt; i++)
		{
			//control.T = 1000.0;
			//control.nsub = tdiss[i];
			//control.tsub = control.T/(double)(control.nsub);
			//control.scint_ts = control.tsub/2.0;
			control.scint_ts = pow(10.0, tdiss[i]);
			control.tsub = control.scint_ts*2.0;
			control.nsub = (int)(control.T/control.tsub);    // nsub >= 2

			for (j=0; j<nf; j++)
			{
				//control.BW = 1000.0;
				//control.chanBW = 100.0;
				//control.scint_freqbw = 50.0;
				//control.nchan = 10;
				control.scint_freqbw = pow(10.0, fdiss[j]);
				control.chanBW = control.scint_freqbw*2.0;
				control.nchan = (int)(control.BW/control.chanBW);  // nchan >= 2

				control.whiteLevel = control.whiteLevel0*sqrt(control.nsub*control.nchan);  

				//printf ("%d %d %lf\n", control.nchan, control.nsub, control.whiteLevel);
		
				calNoise (&noiseStructure, &control);
				//printf ("%lf\n", control.whiteLevel);

				flux0 = control.cFlux0;
				flux1 = control.cFlux1;
				
				acfStructure.probability = 1.0;
					
				// initialize the parameters
				calculateScintScale (&acfStructure, &control);

				// simulate n dynamic spectrum and noise
				seed = TKsetSeed();
				simDynSpec (&acfStructure, seed);

				nMax = 0;
				while (fabs(acfStructure.probability-0.8) >= control.precision && nMax <= 100)
				{
					control.cFlux = flux0+(flux1-flux0)/2.0;
					//printf ("%lf %lf %.8lf %.3f\n", tdiff, fdiff, control.cFlux, acfStructure.probability);
		
					// simulate dynamic spectra
					calculateNDynSpec (&acfStructure, &control, &noiseStructure);

					//if (control.noplot==0 && n == 1)
					//{
					//	// plot while simulating
					//	heatMap (&acfStructure, dname);
					//}

					// merged into calculateNDynSpec
					//qualifyVar (&acfStructure, &noiseStructure, &control);

					if (acfStructure.probability>0.8)
					{
						flux1 = control.cFlux;
					}
					else 
					{
						flux0 = control.cFlux;
					}
					nMax++;
					//printf ("%lf %f %d\n", control.cFlux, acfStructure.probability, nMax);
				}
					
				//printf ("%d %lf %lf %lf %f %d\n", control.nsub, control.whiteLevel, noiseStructure.detection, control.cFlux, acfStructure.probability, nMax);
				printf ("%lf %lf %lf %lf %f %lf %d\n", control.scint_ts, control.scint_freqbw, control.whiteLevel, noiseStructure.detection, control.cFlux, acfStructure.probability, nMax);
				fflush (stdout);

				deallocateMemory (&acfStructure);
			}
		}

		//MPI_Finalize ();

		printf ("Threshold: %f \n", noiseStructure.detection);
		// deallocate memory
	}
	//else
	//{
	//	plotDynSpec(pname, dname);
	//}

	/*
	if (fclose(fin))
	{
		printf ("Can't close file...\n");
		exit(1);
	}
	*/

	free(tdiss);
	free(fdiss);
	return 0;
}

