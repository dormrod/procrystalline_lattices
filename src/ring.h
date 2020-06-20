#ifndef PROCRYSTAL_RING_H
#define PROCRYSTAL_RING_H

#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

using namespace std;

class Ring {
    //Ring in lattice with connectivity information and associated nodes

private:

public:

    //Static variables
    static int autoId,totalActive;

    //Data members
    int id;
    int environment;
    bool active;
    bool selfInteracting;
    VecR<int> cnxs,nodes,arcCnxs,arcNodes;
    VecF<double> crd;

    //Constructor
    Ring();
    Ring(int maxCnxs);

    //Member functions
    void addNode(int nId);
    void addCnx(int cId);
    int getNumCnxs();
    void setCrd(VecF<double> c);
    VecF<string> getEdges();
//    Ring operator + (Ring& source);
    Ring merge(Ring& other, int nid0, int nid1);
    bool checkEdge(int nid0, int nid1);
    void clear();
};


#endif //PROCRYSTAL_RING_H
