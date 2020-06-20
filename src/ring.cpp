#include "ring.h"

int Ring::autoId = 0;
int Ring::totalActive = 0;

Ring::Ring() {}


Ring::Ring(int maxCnxs) {

    //Set incremental id
    active=true;
    selfInteracting=false;
    id=Ring::autoId;
    ++Ring::autoId;
    ++Ring::totalActive;
    environment=-1;

    //Setup containers
    active=true;
    cnxs=VecR<int>(0,maxCnxs);
    nodes=VecR<int>(0,maxCnxs);
    crd=VecF<double>(2);
}


void Ring::addNode(int nId) {
    //Add node id

    nodes.addValue(nId);
}


void Ring::addCnx(int cId) {
    //Add connection id

    cnxs.addValue(cId);
}

int Ring::getNumCnxs() {
    //Get number of connections

    return cnxs.n;
}

void Ring::setCrd(VecF<double> c) {
    //Set coordinate

    crd=c;
}

VecF<string> Ring::getEdges() {
    //Get edge codes

    VecF<string> edges(nodes.n);
    int nid0,nid1;
    for(int i=0; i<nodes.n-1; ++i){
        nid0=nodes[i];
        nid1=nodes[i+1];
        if(nid0<nid1) edges[i]="#"+to_string(nid0)+"#"+to_string(nid1);
        else edges[i]="#"+to_string(nid1)+"#"+to_string(nid0);
    }
    nid0=nodes[0];
    nid1=nodes[nodes.n-1];
    if(nid0<nid1) edges[nodes.n-1]="#"+to_string(nid0)+"#"+to_string(nid1);
    else edges[nodes.n-1]="#"+to_string(nid1)+"#"+to_string(nid0);

    return edges;
}

//Ring Ring::operator+(Ring &source) {
//    //New ring with combined nodes
//
//    Ring ring(this->cnxs.nMax);
//    VecR<int> nonuniqueNodes(0,this->nodes.n+source.nodes.n);
//    for(int i=0; i<this->nodes.n; ++i) nonuniqueNodes.addValue(this->nodes[i]);
//    for(int i=0; i<source.nodes.n; ++i) nonuniqueNodes.addValue(source.nodes[i]);
//    VecR<int> uniqueNodes=vUnique(nonuniqueNodes);
//    for(int i=0; i<uniqueNodes.n; ++i) ring.addNode(uniqueNodes[i]);
//
//    return ring;
//}

Ring Ring::merge(Ring &other, int nid0, int nid1) {
    //Combine two rings to form new ring

    Ring ring(cnxs.nMax);
    int posa0=-1,posa1=-1,posb0=-1,posb1=-1;
    for(int i=0; i<nodes.n; ++i){
        if(nodes[i]==nid0) posa0=i;
        else if(nodes[i]==nid1) posa1=i;
    }
    for(int i=0; i<other.nodes.n; ++i){
        if(other.nodes[i]==nid0) posb0=i;
        else if(other.nodes[i]==nid1) posb1=i;
    }
    bool aForwards,bForwards;
    if(abs(posa1-posa0)==1){
        if(posa1-posa0<0) aForwards=true;
        else aForwards=false;
    }
    else{
        if(posa0==0) aForwards=true;
        else aForwards=false;
    }
    if(abs(posb1-posb0)==1){
        if(posb1-posb0<0) bForwards=false;
        else bForwards=true;
    }
    else{
        if(posb0==0) bForwards=false;
        else bForwards=true;
    }
    if(aForwards){
        for(int i=0; i<nodes.n-1; ++i) ring.addNode(nodes[(posa0+i)%nodes.n]);
    }
    else{
        for(int i=0; i<nodes.n-1; ++i) ring.addNode(nodes[(posa0+nodes.n-i)%nodes.n]);
    }
    if(bForwards){
        for(int i=0; i<other.nodes.n-1; ++i) ring.addNode(other.nodes[(posb1+i)%other.nodes.n]);
    }
    else{
        for(int i=0; i<other.nodes.n-1; ++i) ring.addNode(other.nodes[(posb1+other.nodes.n-i)%other.nodes.n]);
    }
    ring.crd = crd;

    return ring;
}


bool Ring::checkEdge(int nid0, int nid1) {
    //Check if given node ids adjacent in ring

    VecR<int> pos0(0,4);
    for(int i=0; i<nodes.n; ++i){
        if(nodes[i]==nid0) pos0.addValue(i);
    }

    if(pos0.n==0) return false;
    for(int i=0; i<pos0.n; ++i){
        if(nodes[(pos0[i]+1)%nodes.n]==nid1) return true;
        else if(nodes[(pos0[i]+nodes.n-1)%nodes.n]==nid1) return true;
    }
    return false;
}


void Ring::clear() {
    //Deactivate and archive arrays

    active=false;
    arcCnxs = cnxs;
    arcNodes = nodes;
    cnxs=VecR<int>(0);
    nodes=VecR<int>(0);
    --Ring::totalActive;
}