#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "sampler.h"

using namespace std;

int main() {

    //Set up logfile
    Logfile logfile("./procrystal.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("Procrystal Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input parameters
    logfile.write("Reading input file parameters");
    ifstream inputFile("./procrystal.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file procrystal.inpt");
    string skip,line;
    ++logfile.currIndent;
    //Procrystal parameters
    string latticeType; //square or triangular lattice
    int latticeCnd; //underlying lattice coordination
    int nodeCnd; //coordination number of nodes (after MC)
    string linkType; //whether links random or patterned
    bool linkPattern=false; //whether pattern selected
    VecF<int> patternR, patternL; //right and left patterns
    int latticeDim; //lattice dimension
    logfile.write("Reading procrystal parameters");
    ++logfile.currIndent;
    for(int i=0; i<3; ++i) getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>latticeType;
    if(latticeType=="sq" or latticeType=="trihex") latticeCnd=4;
    else if(latticeType=="snub" or latticeType=="isosnub") latticeCnd=5;
    else if(latticeType=="tri") latticeCnd=6;
    else logfile.criticalError("Lattice type not valid");
    logfile.write("Lattice type:",latticeType);
    getline(inputFile,line);
    istringstream(line)>>nodeCnd;
    logfile.write("Node coordination:",nodeCnd);
    getline(inputFile,line);
    istringstream(line)>>linkType;
    getline(inputFile,line);
    if(linkType.substr(0,7)=="pattern"){
        linkPattern=true;
        patternR=VecF<int>(latticeCnd);
        patternL=VecF<int>(latticeCnd);
        istringstream ss(line);
        for(int i=0; i<latticeCnd; ++i) ss>>patternR[i];
        if(vSum(patternR)!=nodeCnd) logfile.criticalError("Pattern does not match node coordination");
        if(linkType=="pattern2"){
            for(int i=0, j=latticeCnd-1; i<latticeCnd; ++i, --j) patternL[i]=patternR[j];
        }
        else patternL=patternR;
        logfile.write("Link pattern selected");
    }
    else logfile.write("Random links selected");
    getline(inputFile,line);
    istringstream(line)>>latticeDim;
    logfile.write("Lattice dimension:",latticeDim);
    --logfile.currIndent;
    //Monte Carlo parameters
    int randomSeed; //random seed
    double temperature; //MC temperature
    int samples; //number of samples to generate
    logfile.write("Reading Monte Carlo parameters");
    ++logfile.currIndent;
    for(int i=0; i<2; ++i) getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>randomSeed;
    logfile.write("Random seed: ",randomSeed);
    getline(inputFile,line);
    istringstream(line)>>temperature;
    logfile.write("Temperature: ",temperature);
    getline(inputFile,line);
    istringstream(line)>>samples;
    logfile.write("Samples: ",samples);
    --logfile.currIndent;
    //Output
    logfile.write("Reading output parameters");
    string outPrefix; //prefix for output files
    int writeSamples; //flag to write sample coordinates and connectivities
    int writeEnvironments; //flag to write node environments
    int writeRDFs; //flag to calculate and write radial distribution functions
    double rdfDelta; //histogram width for rdf
    int writeSk; //flag to calculate structure factor
    int skMaxN; //maximum n value for structure factor
    double skDelta; //width for sk
    ++logfile.currIndent;
    for(int i=0; i<2; ++i) getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>outPrefix;
    logfile.write("Output prefix: ",outPrefix);
    getline(inputFile,line);
    istringstream(line)>>writeSamples;
    logfile.write("Write samples: ",writeSamples);
    getline(inputFile,line);
    istringstream(line)>>writeEnvironments;
    logfile.write("Write clusters: ",writeEnvironments);
    getline(inputFile,line);
    istringstream(line)>>writeRDFs;
    logfile.write("Calculate and write RDFs: ",writeRDFs);
    getline(inputFile,line);
    istringstream(line)>>rdfDelta;
    logfile.write("RDF delta: ",rdfDelta);
    getline(inputFile,line);
    istringstream(line)>>writeSk;
    logfile.write("Calculate and write structure factor: ",writeSk);
    getline(inputFile,line);
    istringstream(line)>>skDelta;
    logfile.write("Sk delta ",skDelta);
    getline(inputFile,line);
    istringstream(line)>>skMaxN;
    logfile.write("Sk max n ",skMaxN);

    logfile.currIndent-=2;
    logfile.separator();

    //Initialise sampler
    Sampler sampler;
    sampler.setProcrystal(latticeType,nodeCnd,linkPattern,patternR,patternL,latticeDim);
    sampler.setMonteCarlo(randomSeed,temperature,samples);
    sampler.setOutput(outPrefix,writeSamples,writeEnvironments,writeRDFs,rdfDelta,writeSk,skDelta,skMaxN);

    //Run sampler
//    sampler.bruteForce(logfile);
    sampler.sample(logfile);

    return 0;
}