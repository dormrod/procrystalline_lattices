#ifndef PROCRYSTAL_LATTICE_H
#define PROCRYSTAL_LATTICE_H

#include <map>
#include <random>
#include <string>
#include <map>
#include <unordered_map>
#include <sstream>
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "node.h"
#include "ring.h"
#include "chain.h"
#include "outputfile.h"

using namespace std;

class Lattice {

private:

    //Data members
    int latticeCnd,nodeCnd;
    bool pattern;
    VecF<int> patternR,patternL;
    VecR<Node> nodes;
    VecR<Ring> rings;
    VecR<Chain> chains;
    mt19937 mtGen;
    uniform_int_distribution<int> randRL;
    uniform_real_distribution<double> rand01;
    bool converged;
    VecF<double> periodicity;
    bool calcEnvs;
    string latticeCode;
    bool calcRings,calcChains;

    //Member functions
    void initialiseSqLattice(int dim, bool restart=false);
    void initialiseTriLattice(int dim, bool restart=false);
    void initialiseSnubSqLattice(int dim, bool restart=false);
    void initialiseIsoSnubQuadLattice(int dim, bool restart=false);
    void initialiseTriHexLattice(int dim, bool restart=false);
    void initialiseHexLattice(int dim, bool restart=false);
    void initialiseProLattice();
    int findRings();
    int findChains();
    void findEnvironments();
    void calculateRDF(VecF<double> &x, VecF<double> &y, VecF<double> &rdf, double delta);
    void calculateSk(VecF<double> &x, VecF<double> &y, VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax);
    VecF<double> periodicCentreOfMass(VecF<double> &x, VecF<double> &y);

public:

    //Constructor
    Lattice();
    Lattice(int seed);

    //Member functions
    void initialise(string latType, int cnd, int latDim, bool pat, VecF<int> patR, VecF<int> patL, string prefix, bool restart=false, bool crystal=false, bool envs=false);
    VecF<int> generate(int maxIterations, double temperature);
    int generate(VecR<int> &pairA, VecR<int> &pairB);
    VecF<double> networkAnalysis(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk);
    void chainAnalysis(VecR<int> &chainLengths);
    VecF<double> calculateNetworkProperties(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk);
    void rdfAnalysis(VecF<double> &latticeRDF, VecF<double> &dualRDF, double delta);
    void nodeAnalysis(VecF<int> &pt, VecF< VecF<int> > &pst, VecF< VecF< VecF<int> > > &prst, OutputFile &testFile);
    void skAnalysis(VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax);
    int getNumNodes();
    int getNumActiveRings();
    VecF<double> getPeriodicity();
    VecF<int> getEnvironments();
    string getEnvironmentCode();
    void writeCrds(string prefix);
    void writeNetwork(string prefix, int index);
    void getEdgeCombinations(VecR<int> &pairA, VecR<int> &pairB, int &n, int &r);
};


#endif //PROCRYSTAL_LATTICE_H
