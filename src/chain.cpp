#include "chain.h"

int Chain::autoId=0;


Chain::Chain() {}


Chain::Chain(int maxNodes) {

    //Set incremental id
    id=Chain::autoId;
    ++Chain::autoId;

    //Setup containers
    nodes=VecR<int>(0,maxNodes);
}


void Chain::addNode(int nId) {
    //Add node id

    nodes.addValue(nId);
}
