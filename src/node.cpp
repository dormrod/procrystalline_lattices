#include "node.h"

int Node::autoId = 0;

Node::Node() {}


Node::Node(int maxCnxs) {
    //Construct with maximum allowed nodes

    //Set incremental id
    id=Node::autoId;
    ++Node::autoId;

    //Initialise containters
    allowedCnxs=VecR<int>(0,maxCnxs);
    cnxs=VecR<int>(0,maxCnxs);
    rings=VecR<int>(0,maxCnxs);
    crd=VecF<double>(2);
    environment=-1;
    cluster=-1;
}


void Node::resetRings() {
    //Clear ring vector

    rings=VecR<int>(0,allowedCnxs.n);
}


void Node::setCrd(VecF<double> c) {
    //Set coordinate

    crd=c;
}


void Node::addCnx(int cId) {
    //Add allowed connection id

    allowedCnxs.addValue(cId);
}


void Node::addRing(int rId) {
    //Add ring id

    rings.addValue(rId);
}

void Node::randomCnxs(int numCnxs, mt19937 &mtGen) {
    //Shuffle available connections and pick given number

    cnxs=vShuffle(allowedCnxs,mtGen);
    cnxs.setSize(numCnxs);
}

void Node::randomCnxs(int numCnxs, VecF<int> pattern, mt19937 &mtGen) {
    //Pick connections in orientation of supplied pattern

    int n=allowedCnxs.n;
    uniform_int_distribution<int> randPos(0,n);
    int patternStart=randPos(mtGen);
    int j=0;
    for(int i=0; i<n; ++i){
        if(pattern[i]==1){
            cnxs[j]=allowedCnxs[(patternStart+i)%n];
            ++j;
        }
    }
    cnxs.setSize(numCnxs);
}

void Node::setCnxs(VecR<int> ids) {
    //Set connections explicitly

    cnxs=ids;
}


void Node::maximiseCnxs() {
    //Set connections to all allowed

    cnxs=allowedCnxs;
}


void Node::breakCnx(int id) {
    //Break connection

    cnxs.delValue(id);
}


void Node::swapRing(int delId, int addId) {
    //Swap all instances of ring id

    for(int i=0; i<rings.n; ++i){
        if(rings[i]==delId) rings[i]=addId;
    }
    rings=vUnique(rings);
}


VecR<int> Node::getReciprocalCnxs(VecR<Node> &others) {
    //Find connections which are reciprocated by other nodes

    VecR<int> rCnxs(0,cnxs.n);
    for(int i=0; i<cnxs.n; ++i){
        if(vContains(others[cnxs[i]].cnxs,id)) rCnxs.addValue(cnxs[i]);
    }
    return rCnxs;
}


VecR<int> Node::getVacantCnxs() {
    //Get allowed connections which are not formed

    VecR<int> vacantCnxs(0, allowedCnxs.n);
    for(int i=0; i<allowedCnxs.n; ++i){
        if(!vContains(cnxs,allowedCnxs[i])) vacantCnxs.addValue(allowedCnxs[i]);
    }
    return vacantCnxs;
}


VecR<int> Node::getSharedRings(Node &other) {
    //Get rings shared with other node

    VecR<int> shared(0,rings.n);
    for(int i=0; i<rings.n; ++i){
        if(vContains(other.rings,rings[i])) shared.addValue(rings[i]);
    }
    return shared;
}

VecR<int> Node::getSharedCnxs(Ring &ring) {
    //Get shared nodes with ring

    VecR<int> shared(0,cnxs.n);
    for(int i=0; i<cnxs.n; ++i){
        if(vContains(ring.nodes,cnxs[i])) shared.addValue(cnxs[i]);
    }
    return shared;
}