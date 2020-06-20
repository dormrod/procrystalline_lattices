#ifndef PROCRYSTAL_NODE_H
#define PROCRYSTAL_NODE_H

#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "ring.h"

using namespace std;

class Node {
    //Node in lattice with connectivity information and associated rings

private:

public:

    //Static variables
    static int autoId;

    //Data members
    int id;
    int environment,cluster;
    VecR<int> allowedCnxs,cnxs,rings;
    VecF<double> crd;

    //Constructors
    Node();
    Node(int maxCnxs);

    //Member functions
    void resetRings();
    void setCrd(VecF<double> c);
    void addCnx(int cId);
    void addRing(int rId);
    void randomCnxs(int numCnxs, mt19937 &mtGen);
    void randomCnxs(int numCnxs, VecF<int> pattern, mt19937 &mtGen);
    void setCnxs(VecR<int> ids);
    void breakCnx(int id);
    void maximiseCnxs();
    void swapRing(int delId,int addId);
    VecR<int> getVacantCnxs();
    VecR<int> getReciprocalCnxs(VecR<Node> &others);
    VecR<int> getSharedRings(Node &other);
    VecR<int> getSharedCnxs(Ring &ring);
};


#endif //PROCRYSTAL_NODE_H
