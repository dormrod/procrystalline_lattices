#ifndef PROCRYSTAL_CHAIN_H
#define PROCRYSTAL_CHAIN_H

#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

using namespace std;
class Chain {
    //Chain of 2-coordinate nodes in lattice

private:

public:

    //Static variables
    static int autoId;

    //Data members
    int id;
    VecR<int> nodes;

    //Constructor
    Chain();
    Chain(int maxNodes);

    //Member functions
    void addNode(int nId);
};


#endif //PROCRYSTAL_CHAIN_H
