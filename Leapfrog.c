#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dot(a,b) ((a[0])*(b[0]) + (a[1])*(b[1]) + (a[2])*(b[2]))
#define lengthSquared(a) dot(a,a)
#define length(a) sqrt(lengthSquared(a))

#define kNumParticles 3
#define kNumDims 3

double masses[kNumParticles] =
 {
   5.9722e24, 1.8986e27, 1.9889e30
 };
double lf_positions[kNumParticles][kNumDims] =
  {
     1.470983e11, 0.000000, 0.000000,
    -7.405736e11, 0.000000, 0.000000,
     0.000000e00, 0.000000, 0.000000
  };
double lf_velocities[kNumParticles][kNumDims] =
  {
      2.694354528161541e03,   3.016946927465788e04,  0.000000,
     -1.0965244901087316e02, -1.3710001990210707e04, 0.000000,
      0.000000,               0.000000,              0.000000
  };
double lf_forces[kNumParticles][kNumDims];

double timeStep = 864000.0;
double epsilonSquared = 0.250000;

int nSteps,
  currStep;

/* Gravitational constant */

double g = 6.67384e-11;

int numParticles = kNumParticles;

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

  lf_initializeParticles();

  nSteps = 1000000;

  for (currStep = 0; currStep < nSteps; currStep++) {
    updateVelocities();
    updatePositions();
    updateForces();
  }

  writeAllParticleData();
  printf("\n");
  lf_outputData();
  return 0;
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
    }
    printf("          kinetic energy: %14.7g\n",
           totalEnergy);

    /* potential energy */
    for(i=0; i<numParticles; i++) {
      for(j=0; j<i; j++) {
	for(k=0; k<kNumDims; k++)
	  deltaPos[k] = lf_positions[i][k] - lf_positions[j][k];
	totalEnergy -= masses[i]*masses[j] / sqrt(lengthSquared(deltaPos) + epsilonSquared);
      }
    }

    printf("          total energy: %14.7g\n",
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
  for(i=0; i<numParticles; i++) {
    printf("%i\n"
	   "      f = <%16.10e %16.10e %16.10e>\n"
	   "    pos = <%16.10e %16.10e %16.10e>\n"
	   "    vel = <%16.10e %16.10e %16.10e>\n",
	   i,
	   lf_forces[i][0], lf_forces[i][1], lf_forces[i][2],
	   lf_positions[i][0], lf_positions[i][1], lf_positions[i][2],
	   lf_velocities[i][0], lf_velocities[i][1], lf_velocities[i][2]
	   );
      }
}
