#include "sampler.h"

Sampler::Sampler() {
    //Default constructor
}


void Sampler::setProcrystal(string latType, int cnd, bool pat, VecF<int> patR, VecF<int> patL, int latDim) {
    //Set procrystal parameters

    latticeType=latType;
    nodeCnd=cnd;
    pattern=pat;
    patternR=patR;
    patternL=patL;
    latticeDim=latDim;
}


void Sampler::setMonteCarlo(int seed, double temp, int samples) {
    //Set monte carlo parameters

    randomSeed=seed;
    temperature=temp;
    numSamples=samples;
}


void Sampler::setOutput(string outPrefix, int write, int envs, int rdf, double delta, int sk, double sdelta, int smaxn) {
    //Set output parameters

    outputPrefix=outPrefix;
    if(write==1) writeSamples=true;
    else writeSamples=false;
    if(envs==1) writeEnvs=true;
    else writeEnvs=false;
    if(rdf==1){
        writeRDFs=true;
        rdfDelta=delta;
    }
    else writeRDFs=false;
    if(sk==1){
        writeSk=true;
        skDelta=sdelta;
        skMaxN=smaxn;
    }
    else writeSk=false;
    if(nodeCnd>2){
        writeRingAn=true;
        writeChainAn=false;
    }
    else{
        writeRingAn=false;
        writeChainAn=true;
    }
}


void Sampler::sample(Logfile &logfile) {
    //Sample procrystalline lattices

    logfile.write("Procrystal Generation - Monte Carlo");
    ++logfile.currIndent;

    //Initialise procrystalline lattice
    Lattice procrystal(randomSeed);
    procrystal.initialise(latticeType,nodeCnd,latticeDim,pattern,patternR,patternL,outputPrefix,false,false,writeEnvs);
    int maxIterations = pow(procrystal.getNumNodes(),3);
    VecF<double> pbc=procrystal.getPeriodicity();
//    maxIterations=10000;
    cout<<procrystal.getNumNodes()<<endl;
    int numNodes=procrystal.getNumNodes();

    //Setup analysis containers and output file
    int maxK=2*numNodes;
    VecF<int> k(maxK);
    VecF<int> pkTotal(maxK);
    VecF< VecF<int> > ejkTotal(maxK);
    for(int i=0; i<maxK; ++i) k[i]=i;
    for(int i=0; i<maxK; ++i) ejkTotal[i]=VecF<int>(maxK);
    OutputFile netFile(outputPrefix+"_net.dat");
    netFile.initVariables(10,4,60,20);
    OutputFile rdfFile(outputPrefix + "_rdf.dat");
    int maxBin=0;
    if(writeRDFs) {
        if (pbc[0] <= pbc[1]) maxBin = floor(pbc[0] / (2 * rdfDelta)) + 1;
        else maxBin = floor(pbc[1] / (2 * rdfDelta)) + 1;
    }
    VecF<double> latticeRDF(maxBin), dualRDF(maxBin);
    VecF< VecF< map<int,int> > >clusters;
    VecF<int> clusters2;
    OutputFile skFile(outputPrefix + "_sk.dat");
    int maxBinSk=0;
    if(writeSk) {
        maxBinSk = floor((2 * M_PI / pbc[0]) * skMaxN * (2 * M_PI / pbc[1]) * skMaxN / skDelta) + 1;
    }
    VecF<double> dualSk(maxBinSk);
    VecF<int> multiplicity(maxBinSk);
    OutputFile envFile(outputPrefix+"_env.dat");
    envFile.initVariables(10,4,60,20);
    VecR<string> environments(0,numSamples);
    OutputFile chainFile(outputPrefix + "_chain.dat");
    VecR<int> chainLengths(0,numSamples*numNodes/2);

    //Generate required number of samples
    int convergedSamples=0;
    int attemptedSamples=0;
    int averageIterations=0;
    for(;;){
        VecF<int> opt=procrystal.generate(maxIterations,temperature);
        if(opt[0]==1){
            VecF<int> pk;
            VecF< VecF<int> > ejk;
            if(writeRingAn){
                VecF<double> netAnalysis=procrystal.networkAnalysis(k,pk,ejk);
                netFile.writeRowVector(netAnalysis);
                pkTotal+=pk;
                for(int i=0; i<maxK; ++i) ejkTotal[i]+=ejk[i];
            }
            if(writeChainAn){
                procrystal.chainAnalysis(chainLengths);
            }
            if(writeRDFs) procrystal.rdfAnalysis(latticeRDF,dualRDF,rdfDelta);
            if(writeSk) procrystal.skAnalysis(dualSk,multiplicity,skDelta,skMaxN);
            if(writeSamples) procrystal.writeNetwork(outputPrefix,convergedSamples);
            if(writeEnvs) environments.addValue(procrystal.getEnvironmentCode());
            logfile.write("Sample "+to_string(attemptedSamples)+" generated in "+to_string(opt[1])+" iterations");
            ++convergedSamples;
            averageIterations+=opt[1];
        }
        else if(opt[0]==0){
            logfile.write("Sample "+to_string(attemptedSamples)+" did not converge");
        }
        else if(opt[0]==2){
            logfile.write("Sample "+to_string(attemptedSamples)+" converged but sample split");
        }
        else if(opt[0]==3){
            logfile.write("Sample "+to_string(attemptedSamples)+" converged but failed with error");
        }
        ++attemptedSamples;
        cout<<convergedSamples<<"/"<<attemptedSamples<<endl;
        if(convergedSamples==numSamples) break;
        procrystal.initialise(latticeType,nodeCnd,latticeDim,pattern,patternR,patternL,outputPrefix,true,false,writeEnvs);
    }
    cout<<double(averageIterations)/numSamples<<endl;

    //Write combined analysis
    if(writeRingAn) {
        VecF<double> netAnalysis = procrystal.calculateNetworkProperties(k, pkTotal, ejkTotal);
        netFile.writeRowVector(netAnalysis);
        OutputFile pkFile(outputPrefix + "_pk.dat");
        pkFile.initVariables(10, 4, 60, 20);
        double norm = vSum(pkTotal);
        for (int i = 0; i < maxK; ++i) pkFile.write(pkTotal[i] / norm);
    }
    if(writeChainAn) {
        for(int i=0; i<chainLengths.n; ++i) chainFile.write(chainLengths[i]);
    }

    //Normalise RDFs and write
    if(writeRDFs) {
        VecF<double> bins(latticeRDF.n);
        for (int i = 0; i < bins.n; ++i) bins[i] = rdfDelta * (i + 0.5);
        int numNodes = procrystal.getNumNodes();
        int numDual = procrystal.getNumActiveRings();
        double normLat, normDual;
        normLat = numNodes * (numNodes / (pbc[0] * pbc[1]) * M_PI * numSamples);
        normDual = numDual * (numDual / (pbc[0] * pbc[1]) * M_PI * numSamples);
        for (int i = 0; i < bins.n; ++i) {
            latticeRDF[i] /= normLat * (pow((i + 1) * rdfDelta, 2) - pow(i * rdfDelta, 2));
            dualRDF[i] /= normDual * (pow((i + 1) * rdfDelta, 2) - pow(i * rdfDelta, 2));
        }
        for (int i = 0; i < bins.n; ++i) rdfFile.write(bins[i], latticeRDF[i], dualRDF[i]);
    }

    //Normalise structure factor and write
    if(writeSk){
        VecF<double> bins(dualSk.n);
        for(int i=0; i<bins.n; ++i) bins[i]=skDelta*(i+0.5);
        dualSk /= numSamples;
        for(int i=0; i<dualSk.n; ++i){
            if(multiplicity[i]>0) dualSk[i]/=multiplicity[i];
        }
        for (int i = 0; i < bins.n; ++i){
            if(multiplicity[i]>0) skFile.write(bins[i], dualSk[i]);
        }
    }

    //Write environments
    if(writeEnvs){
        for(int i=0; i<environments.n; ++i) envFile.write(environments[i]);
    }

    --logfile.currIndent;
    logfile.separator();
}


void Sampler::bruteForce(Logfile &logfile) {
    //Find all combinations of procrystalline lattices

    logfile.write("Procrystal Generation - Brute Force");
    ++logfile.currIndent;

    //Initialise procrystalline lattice
    Lattice procrystal(randomSeed);
    procrystal.initialise(latticeType,nodeCnd,latticeDim,pattern,patternR,patternL,outputPrefix,false,true,writeEnvs);
    int numNodes=procrystal.getNumNodes();
    cout<<numNodes<<endl;

    //Setup analysis containers and output file
    int maxK=4*latticeDim;
    VecF<int> k(maxK);
    VecF<int> pkTotal(maxK);
    VecF< VecF<int> > ejkTotal(maxK);
    for(int i=0; i<maxK; ++i) k[i]=i;
    for(int i=0; i<maxK; ++i) ejkTotal[i]=VecF<int>(maxK);
    OutputFile netFile(outputPrefix+"_net.dat");
    netFile.initVariables(10,4,60,20);

    //Get edge pairs
    int n,r;
    VecR<int> pairA,pairB;
    procrystal.getEdgeCombinations(pairA,pairB,n,r);
    for(int i=0; i<n; ++i){
        cout<<i<<" "<<pairA[i]<<" "<<pairB[i]<<endl;
    }

    //Brute force all combinations of pairs
    vector<bool> combination(n);
    fill(combination.begin(), combination.begin() + r, true);
    VecF<int> nodeCount(numNodes);
    int bfCount=0,acceptCount=0;
    do {
        nodeCount=0;
        for (int i = 0; i < n; ++i) {
            if (combination[i]) {
                ++nodeCount[pairA[i]];
                ++nodeCount[pairB[i]];
            }
        }
        bool accept=true;
        for(int i=0; i<numNodes; ++i){
            if(nodeCount[i]!=1) accept = false;
        }
        //Generate procrystalline lattice
        if(accept){
            VecR<int> breakPairA(0,r),breakPairB(0,r);
            for (int i = 0; i < n; ++i) {
                if (combination[i]) {
                    breakPairA.addValue(pairA[i]);
                    breakPairB.addValue(pairB[i]);
                }
            }
            procrystal.initialise(latticeType,nodeCnd,latticeDim,pattern,patternR,patternL,outputPrefix,true,true,writeEnvs);
            procrystal.generate(breakPairA,breakPairB);
            VecF<int> pk;
            VecF< VecF<int> > ejk;
            VecF<double> netAnalysis=procrystal.networkAnalysis(k,pk,ejk);
            pkTotal+=pk;
            for(int i=0; i<maxK; ++i) ejkTotal[i]+=ejk[i];
            netFile.writeRowVector(netAnalysis);
            if(writeSamples) procrystal.writeNetwork(outputPrefix,acceptCount);
            logfile.write("Configuration found");
            procrystal.initialise(latticeType,nodeCnd,latticeDim,pattern,patternR,patternL,outputPrefix,true,writeEnvs);
            ++acceptCount;
        }
        ++bfCount;
        if(bfCount%1000000==0) cout<<bfCount<<endl;
    } while (prev_permutation(combination.begin(), combination.end()));
    logfile.write("Total found / attempted",acceptCount,bfCount);
    cout<<acceptCount<<endl;
    cout<<bfCount<<endl;
    //Write combined analysis
    VecF<double> netAnalysis=procrystal.calculateNetworkProperties(k,pkTotal,ejkTotal);
    netFile.writeRowVector(netAnalysis);
    OutputFile pkFile(outputPrefix+"_pk.dat");
    pkFile.initVariables(10,4,60,20);
    double norm=vSum(pkTotal);
    for(int i=0; i<maxK; ++i) pkFile.write(pkTotal[i]/norm);

    --logfile.currIndent;
    logfile.separator();
}
