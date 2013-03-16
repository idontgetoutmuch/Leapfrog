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

/* per-particle physical and dynamic properties */

double masses[kMaxParticles],
  positions[kMaxParticles][kNumDims],
  lf_positions[kMaxParticles][kNumDims],
  lf_velocities[kMaxParticles][kNumDims],
  previousPositions[kMaxParticles][kNumDims],
  previousVelocities[kMaxParticles][kNumDims],
  forces[kMaxParticles][kNumDims],
  lf_forces[kMaxParticles][kNumDims],
  forceDots[kMaxParticles][kNumDims];

/* other per-particle state variables */

double steps[kMaxParticles],
  t0[kMaxParticles],
  t1[kMaxParticles],
  t2[kMaxParticles],
  t3[kMaxParticles],
  d1[kMaxParticles][kNumDims],
  d2[kMaxParticles][kNumDims],
  d3[kMaxParticles][kNumDims];

double currentTime = 0.0,
  nextTime = 0.0,
  timeStep,
  finalTime,
  eta,
  epsilonSquared;

int nSteps;

/* Gravitational constant */

double g = 6.67384e-11;

int numSteps = 0,
  numParticles;

void readParameters();
void readParticles();
void lf_readParticles();
void initializeParticles();
void lf_initializeParticles();
void outputData();
void advanceParticles();
void updateForce();
void updateVelocities();
void updatePositions();
void writeAllParticleData();

FILE * inputFile;
FILE * lf_inputFile;
FILE * outputFile;

int main()
{
  int i;
  inputFile = fopen("initc.data", "r");
  readParameters();
  readParticles();
  fclose(inputFile);

  lf_inputFile = fopen("initc.data", "r");
  readParameters();
  lf_readParticles();
  fclose(inputFile);

  outputFile = fopen("sphere.data", "w");
  lf_initializeParticles();

  nSteps = (int) finalTime / timeStep;
  printf("Num of steps %14.7g\n", (int) finalTime / timeStep);

  for (i = 0; i < nSteps; i++) {

    for (i = 0; i < numParticles; i++) {
      printf("F_x = %14.7g, F_y = %14.7g, F_z = %14.7g\n",
	     lf_forces[i][0],
	     lf_forces[i][1],
	     lf_forces[i][2]
	     );
    }
    for (i = 0; i < numParticles; i++) {
      printf("V_x = %14.7g, V_y = %14.7g, V_z = %14.7g\n",
	     lf_velocities[i][0],
	     lf_velocities[i][1],
	     lf_velocities[i][2]
	     );
    }

    updateVelocities();

    for (i = 0; i < numParticles; i++) {
      printf("V_x = %14.7g, V_y = %14.7g, V_z = %14.7g\n",
	     lf_velocities[i][0],
	     lf_velocities[i][1],
	     lf_velocities[i][2]
	     );
    }

    for (i = 0; i < numParticles; i++) {
      printf("P_x = %14.7g, P_y = %14.7g, P_z = %14.7g\n",
	     lf_positions[i][0],
	     lf_positions[i][1],
	     lf_positions[i][2]
	     );
    }

    updatePositions();

    for (i = 0; i < numParticles; i++) {
      printf("P_x = %14.7g, P_y = %14.7g, P_z = %14.7g\n",
	     lf_positions[i][0],
	     lf_positions[i][1],
	     lf_positions[i][2]
	     );
    }
  }

  fclose(outputFile);

  return 0;
}


void readParameters()
{
    printf("Enter numParticles, eta, timeStep, finalTime and epsilonSquared: ");
    fscanf(inputFile, "%i", &numParticles);
    if(numParticles > kMaxParticles)
    {
        fprintf(stderr, "The program currently supports up to %i particles.\n"
                        "If you reqire more particles, please increase\n"
                        "kMaxParticles in nbody.c and recompile.\n", kMaxParticles);
        exit(1);
    }
    fscanf(inputFile, "%lf", &eta);
    fscanf(inputFile, "%lf", &timeStep);
    fscanf(inputFile, "%lf", &finalTime);
    fscanf(inputFile, "%lf", &epsilonSquared);
    printf("\n");
}

void readParticles()
{
    int i, j;
    for(i=0; i<numParticles; i++)
    {
        fscanf(inputFile, "%lf", &masses[i]);
        for(j=0; j<kNumDims; j++)
            fscanf(inputFile, "%lf", &previousPositions[i][j]);
        for(j=0; j<kNumDims; j++)
            fscanf(inputFile, "%lf", &previousVelocities[i][j]);
    }
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

    printf("Initializing particles %i %i\n", numParticles, kNumDims);
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

void outputData()
{
    double totalEnergy = 0.0, totalMass = 0.0;
    double centerOfMass[kNumDims], totalMomentum[kNumDims],
        totalAngularMomentum[kNumDims];

    double refPos[kNumDims], temp[kNumDims];

    double f2Dot[kNumDims], velocity[kNumDims], deltaPos[kNumDims];

    int i,j,k;

    double dt;

    /* kinetic energy */

    for(i=0; i<numParticles; i++)
    {
        dt = nextTime - t0[i];
        for(k=0; k<kNumDims; k++)
        {
            f2Dot[k] = d3[i][k] * ((t0[i]-t1[i]) + (t0[i]-t2[i])) + d2[i][k];
            positions[i][k] = ((((0.05*d3[i][k]*dt + f2Dot[k]/12.0)*dt + forceDots[i][k])
                *dt + forces[i][k])*dt + previousVelocities[i][k])*dt + previousPositions[i][k];
            velocity[k] = (((0.25*d3[i][k]*dt + f2Dot[k]/3.0)*dt + 3.0*forceDots[i][k])
                *dt + 2.0*forces[i][k])*dt + previousVelocities[i][k];

        }

        /* write the particle to the file */

        totalEnergy += 0.5 * masses[i] * lengthSquared(velocity);

        totalMass += masses[i];

        fprintf(outputFile, "%14.7g %14.7g %14.7g %14.7g \n",
	  masses[i], positions[i][0], positions[i][1], positions[i][2]);

    }

    /* potential energy */
    for(i=0; i<numParticles; i++)
    {
        for(j=0; j<i; j++)
        {
            for(k=0; k<kNumDims; k++)
                deltaPos[k] = positions[i][k] - positions[j][k];
            totalEnergy -= masses[i]*masses[j] / sqrt(lengthSquared(deltaPos) + epsilonSquared);
        }
    }


    printf("\n"
           "            time: %14.7g\n"
           "           steps: %14i\n"
           "          energy: %14.7g\n",
           currentTime,
           numSteps,
           totalEnergy);
}

void updateForce() {
  double d, d2, d2_softened, d3_2;
  int i, j, k;
  for (i = 0; i < numParticles; i++) {
    for (j = 0; j < numParticles; j++) {
      d2 = 0;
      for (k = 0; k < kNumDims; k++) {
	d   = lf_positions[i][k] - lf_positions[j][k];
	d2 += d * d;
      }
      d2_softened = d2 + epsilonSquared;
      d3_2 = d2_softened * sqrt(d2_softened);
      for (k = 0; k < kNumDims; k++) {
	lf_forces[i][k] += -g * masses[i] * masses[j] / d3_2;
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

void advanceParticles()
{
    double s, dt, dt1,dt2, dt3, t1pr, t2pr, t3pr;

    int i,j,k;

    double deltaPos[kNumDims], f1Dot[kNumDims], f2Dot[kNumDims], f3Dot[kNumDims],
        f1[kNumDims];
    double a,b, temp1, temp2, temp3, temp4;

    printf("\n"
           "Advancing Particles to goal time: %10.7f\n", nextTime);
    printf("                            steps       time\n");
    /* this is a hack for printing out the current state */
    printf("                                            ");

    while(currentTime < nextTime)
    {
        currentTime = 1e10;
        for(j=0; j<numParticles; j++)
        {
            if(currentTime > t0[j] + steps[j])
            {
                i = j;
                currentTime = t0[j] + steps[j];
            }
        }

        /* predict all coordinates to first order in force derivative */

        for(j=0; j<numParticles; j++)
        {
            s = currentTime - t0[j];
            for(k=0;k<kNumDims; k++)
                positions[j][k] = ((forceDots[j][k]*s + forces[j][k])*s
                    + previousVelocities[j][k])*s + previousPositions[j][k];
        }

        /* include 2nd and 3rd order and obtain the velocity */

        dt = currentTime - t0[i];

        for (k=0; k<kNumDims; k++)
        {
            f2Dot[k] = d3[i][k] * ((t0[i]-t1[i]) + (t0[i]-t2[i])) + d2[i][k];
            positions[i][k] = (d3[i][k]*0.05*dt + f2Dot[k]/12.0)
                *(dt*dt*dt*dt) + positions[i][k];
            previousVelocities[i][k] = (((d3[i][k]*0.25*dt + f2Dot[k]/3.0)*dt
                + forceDots[i][k]*3.0)*dt + forces[i][k]*2.0)*dt + previousVelocities[i][k];
            f1[k] = 0.0;
        }

        /* obtain current force on i-th body */

        for (j=0; j<numParticles; j++)
        {
            if (j != i)
            {
                for(k=0; k<kNumDims; k++)
                    deltaPos[k] = positions[j][k] - positions[i][k];
                a = 1. / (lengthSquared(deltaPos) + epsilonSquared);
                b = masses[j] * a * sqrt(a);

                for(k=0; k<kNumDims; k++)
                    f1[k] += b*deltaPos[k];
            }
        }

        /* set time intervals for new difference and update the times */

        dt1 = currentTime - t1[i];
        dt2 = currentTime - t2[i];
        dt3 = currentTime - t3[i];
        t1pr = t0[i] - t1[i];
        t2pr = t0[i] - t2[i];
        t3pr = t0[i] - t3[i];
        t3[i] = t2[i];
        t2[i] = t1[i];
        t1[i] = t0[i];
        t0[i] = currentTime;

        /* form new differences and include 4th order semi-iterative */

        for (k=0; k<kNumDims; k++)
        {
            temp1 = (f1[k] - 2.0*forces[i][k]) / dt;
            temp2 = (temp1 - d1[i][k]) / dt1;
            temp3 = (temp2 - d2[i][k]) / dt2;
            temp4 = (temp3 - d3[i][k]) / dt3;
            d1[i][k] = temp1;
            d2[i][k] = temp2;
            d3[i][k] = temp3;
            f1Dot[k] = t1pr * t2pr * t3pr * temp4;
            f2Dot[k] = (t1pr * t2pr + t3pr * (t1pr + t2pr)) * temp4;
            f3Dot[k] = (t1pr + t2pr + t3pr) * temp4;
            previousPositions[i][k] = (((temp4 * dt / 30.0 + f3Dot[k] * 0.05) * dt
                + f2Dot[k] / 12.0) * dt + f1Dot[k] / 6.0) * (dt*dt*dt) + positions[i][k];
            previousVelocities[i][k] = (((temp4 * .2 * dt + f3Dot[k] * 0.25) * dt
                + f2Dot[k] / 3.0) * dt + f1Dot[k] * 0.5) * (dt*dt) + previousVelocities[i][k];
        }

        /* scale F and FDOT by factorials and set new integration step */

        for (k=0; k<kNumDims; k++)
        {
            forces[i][k] = f1[k] * 0.5;
            forceDots[i][k] = ((d3[i][k] * dt1 + d2[i][k]) * dt + d1[i][k]) / 6.0;
            f2Dot[k] = (d3[i][k] * (dt + dt1) + d2[i][k]) * 2.0;
        }
        steps[i] = sqrt(eta * sqrt(lengthSquared(f1) / lengthSquared(f2Dot)));
        numSteps++;
        if(numSteps%100 == 0)
        {
            for(k=0; k<21; k++)
                printf("\b");
            printf("%10i %10.7f", numSteps, currentTime);
            fflush(stdout);
        }
    }
    for(k=0; k<21; k++)
        printf("\b");
    printf("%10i %10.7f", numSteps, currentTime);
    fflush(stdout);
    printf("\n");
}

void writeAllParticleData()
{
    FILE * file;
    int i;
    file = fopen("particlesLocal.data", "w");
    for(i=0; i<numParticles; i++)
    {
        fprintf(file, "%i\n"
                      "    pos = <%f %f %f>\n"
                      "      f = <%f %f %f>\n"
                      "   fdot = <%f %f %f>\n"
                      "prevPos = <%f %f %f>\n"
                      "prevVel = <%f %f %f>\n"
                      "     t0 = %f\n"
                      "     t1 = %f\n"
                      "     t2 = %f\n"
                      "     t3 = %f\n"
                      "   step = %f\n"
                      "     d1 = <%f %f %f>\n"
                      "     d2 = <%f %f %f>\n"
                      "     d3 = <%f %f %f>\n\n", i,
                      positions[i][0], positions[i][1], positions[i][2],
                      forces[i][0], forces[i][1], forces[i][2],
                      forceDots[i][0], forceDots[i][1], forceDots[i][2],
                      previousPositions[i][0], previousPositions[i][1], previousPositions[i][2],
                      previousVelocities[i][0], previousVelocities[i][1], previousVelocities[i][2],
                      t0[i], t1[i], t2[i], t3[i], steps[i],
                      d1[i][0], d1[i][1], d1[i][2],
                      d2[i][0], d2[i][1], d2[i][2],
                      d3[i][0], d3[i][1], d3[i][2]);
    }
    fclose(file);
}
