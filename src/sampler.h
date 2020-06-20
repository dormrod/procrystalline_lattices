#ifndef PROCRYSTAL_SAMPLER_H
#define PROCRYSTAL_SAMPLER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "lattice.h"

class Sampler {
    //Generates procrystal samples

private:

    //Procrystal parameters
    string latticeType;
    bool pattern;
    VecF<int> patternR, patternL;
    int latticeDim, nodeCnd;

    //Monte Carlo parameters
    int randomSeed,numSamples;
    double temperature;

    //Output parameters
    string outputPrefix;
    bool writeSamples,writeRDFs,writeEnvs,writeSk,writeRingAn,writeChainAn;
    int skMaxN;
    double rdfDelta,skDelta;


public:

    //Constructor and setters
    Sampler();
    void setProcrystal(string latType, int cnd, bool pat, VecF<int> patR, VecF<int> patL, int latDim);
    void setMonteCarlo(int seed, double temp, int samples);
    void setOutput(string outPrefix, int write, int envs, int rdf, double delta,  int sk, double sdelta, int smaxn);

    //Start sampler
    void sample(Logfile &logfile);
    void bruteForce(Logfile &logfile);
};


#endif //PROCRYSTAL_SAMPLER_H
