#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dot(a,b) ((a[0])*(b[0]) + (a[1])*(b[1]) + (a[2])*(b[2]))
#define lengthSquared(a) dot(a,a)
#define length(a) sqrt(lengthSquared(a))

#define cross(c,a,b) c[0] = (a[1])*(b[2]) - (a[2])*(b[1]); \
                     c[1] = (a[2])*(b[0]) - (a[0])*(b[2]); \
                     c[2] = (a[0])*(b[1]) - (a[1])*(b[0]);

#define kMaxParticles 2048
#define kNumDims 3

double masses[kMaxParticles],
  positions[kMaxParticles][kNumDims],
  lf_positions[kMaxParticles][kNumDims],
  lf_velocities[kMaxParticles][kNumDims],
  previousPositions[kMaxParticles][kNumDims],
  previousVelocities[kMaxParticles][kNumDims],
  forces[kMaxParticles][kNumDims],
  lf_forces[kMaxParticles][kNumDims],
  forceDots[kMaxParticles][kNumDims];

double currentTime = 0.0,
  nextTime = 0.0,
  timeStep,
  finalTime,
  eta,
  epsilonSquared;

int nSteps,
  currStep;

/* Gravitational constant */

double g = 6.67384e-11;

int numSteps = 0,
  numParticles;

void readParameters();
void lf_readParticles();
void initializeParticles();
void lf_initializeParticles();
void lf_outputData();
void updateForces();
void updateVelocities();
void updatePositions();
void writeAllParticleData();

FILE * inputFile;
FILE * lf_inputFile;
FILE * outputFile;

int main()
{
  int currStep;

  lf_inputFile = fopen("initc.data", "r");
  readParameters();
  lf_readParticles();
  fclose(lf_inputFile);

  outputFile = fopen("sphere.data", "w");
  lf_initializeParticles();

  nSteps = 10; /* (int) finalTime / timeStep; */

  for (currStep = 0; currStep < nSteps; currStep++) {

    printf("\n           steps: %14i\n", currStep);
    lf_outputData();
    updateVelocities();
    updatePositions();
    updateForces();
  }

  fclose(outputFile);

  writeAllParticleData();

  return 0;
}

void readParameters()
{
    printf("Enter numParticles, eta, timeStep, finalTime and epsilonSquared: ");
    fscanf(lf_inputFile, "%i", &numParticles);
    if(numParticles > kMaxParticles)
    {
        fprintf(stderr, "The program currently supports up to %i particles.\n"
                        "If you reqire more particles, please increase\n"
                        "kMaxParticles in nbody.c and recompile.\n", kMaxParticles);
        exit(1);
    }
    fscanf(lf_inputFile, "%lf", &eta);
    fscanf(lf_inputFile, "%lf", &timeStep);
    fscanf(lf_inputFile, "%lf", &finalTime);
    fscanf(lf_inputFile, "%lf", &epsilonSquared);
    printf("\n");
}

void lf_readParticles()
{
    int i, j;
    for(i=0; i<numParticles; i++)
    {
        fscanf(lf_inputFile, "%lf", &masses[i]);
        for(j=0; j<kNumDims; j++)
            fscanf(lf_inputFile, "%lf", &lf_positions[i][j]);
        for(j=0; j<kNumDims; j++)
            fscanf(lf_inputFile, "%lf", &lf_velocities[i][j]);
    }
}

void lf_initializeParticles()
{
    int i,j,k;
    double deltaPos[kNumDims];
    double a, b;

    printf("Initializing particles\n");
    for(i=0; i<numParticles; i++)
    {
        for(k=0; k<kNumDims; k++)
        {
            lf_forces[i][k] = 0.0;
        }
        for(j=0; j<numParticles; j++)
        {
            if(j!=i)
            {
                for(k=0; k<kNumDims; k++)
                {
                    deltaPos[k] = lf_positions[i][k] - lf_positions[j][k];
                }
                a = 1.0/(lengthSquared(deltaPos) + epsilonSquared);
                b = -g * masses[i] * masses[j] * a * sqrt(a);
                for(k=0; k<kNumDims; k++)
                {
                    lf_forces[i][k] += b * deltaPos[k];
                }
            }
        }
    }
    printf("\n");
}

void lf_outputData()
{
    double totalEnergy = 0.0,
      totalMass = 0.0;

    double deltaPos[kNumDims],
      velocity[kNumDims];

    int i,j,k;

    /* kinetic energy */
    for(i=0; i<numParticles; i++) {
      for(k = 0; k < kNumDims; k++)
	velocity[k] = lf_velocities[i][k];
      totalEnergy += 0.5 * masses[i] * lengthSquared(velocity);
      totalMass += masses[i];
      fprintf(outputFile, "%14.7g %14.7g %14.7g %14.7g \n",
	      masses[i], positions[i][0], positions[i][1], positions[i][2]);
    }

    /* potential energy */
    for(i=0; i<numParticles; i++) {
      for(j=0; j<i; j++) {
	for(k=0; k<kNumDims; k++)
	  deltaPos[k] = positions[i][k] - positions[j][k];
	totalEnergy -= masses[i]*masses[j] / sqrt(lengthSquared(deltaPos) + epsilonSquared);
      }
    }

    printf("          energy: %14.7g\n",
           totalEnergy);
}

void updateForces() {
  int i, j, k;
  double deltaPos[kNumDims];
  double a, b;

  for(i=0; i<numParticles; i++) {
    for(k=0; k<kNumDims; k++) {
      lf_forces[i][k] = 0.0;
    }
    for(j=0; j<numParticles; j++) {
      if(j!=i) {
	for(k=0; k<kNumDims; k++) {
	  deltaPos[k] = lf_positions[i][k] - lf_positions[j][k];
	}
	a = 1.0/(lengthSquared(deltaPos) + epsilonSquared);
	b = -g * masses[i] * masses[j] * a * sqrt(a);
	for(k=0; k<kNumDims; k++) {
	  lf_forces[i][k] += b * deltaPos[k];
	}
      }
    }
  }
}

void updateVelocities() {
  int i, k;
  for (i = 0; i < numParticles; i++) {
    for (k = 0; k < kNumDims; k++) {
      lf_velocities[i][k] += lf_forces[i][k] * timeStep / masses[i];
    }
  }
}

void updatePositions() {
  int i, k;
  for (i = 0; i < numParticles; i++) {
    for (k = 0; k < kNumDims; k++) {
      lf_positions[i][k] += lf_velocities[i][k] * timeStep;
    }
  }
}

void writeAllParticleData() {
  FILE * file;
  int i;
  file = fopen("particlesLocal.data", "w");
  for(i=0; i<numParticles; i++) {
    fprintf(file, "%i\n"
	    "      f = <%16.10e %16.10e %16.10e>\n"
	    "    pos = <%16.10e %16.10e %16.10e>\n"
	    "    vel = <%16.10e %16.10e %16.10e>\n",
	    i,
	    lf_forces[i][0], lf_forces[i][1], lf_forces[i][2],
	    lf_positions[i][0], lf_positions[i][1], lf_positions[i][2],
	    lf_velocities[i][0], lf_velocities[i][1], lf_velocities[i][2]
	    );
      }
  fclose(file);
}
