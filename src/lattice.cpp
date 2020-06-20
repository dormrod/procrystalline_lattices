#include "lattice.h"

Lattice::Lattice() {
    //Default constructor

    //Set random number generator
    mtGen.seed(0);
}


Lattice::Lattice(int seed) {
    //Constructor with random seed

    //Set random number generator
    mtGen.seed(seed);
    rand01=uniform_real_distribution<double>(0,1);
    randRL=uniform_int_distribution<int>(0,2);
}


void Lattice::initialise(string latType, int cnd, int latDim, bool pat, VecF<int> patR, VecF<int> patL, string prefix,
                         bool restart, bool crystal, bool envs) {
    //Initialise lattice

    //Set variables
    nodeCnd = cnd;
    pattern = pat;
    patternR = patR;
    patternL = patL;
    Node::autoId = 0;
    Ring::autoId = 0;
    Ring::totalActive = 0;
    Chain::autoId = 0;
    converged=false;
    periodicity=VecF<double>(2);
    calcEnvs=envs;
    latticeCode = latType+to_string(nodeCnd);
    if(nodeCnd>2){
        calcRings=true;
        calcChains=false;
    }
    else{
        calcRings=false;
        calcChains=true;
    }

    //Make regular lattice
    if(latType=="sq") initialiseSqLattice(latDim,restart);
    else if(latType=="tri") initialiseTriLattice(latDim,restart);
    else if(latType=="snub") initialiseSnubSqLattice(latDim,restart);
    else if(latType=="isosnub") initialiseIsoSnubQuadLattice(latDim,restart);
    else if(latType=="trihex") initialiseTriHexLattice(latDim,restart);
    else if(latType=="hex") initialiseHexLattice(latDim,restart);

    //Calculate ring coordinates
    for(int i=0; i<rings.n; ++i){
        VecF<double> x(rings[i].nodes.n),y(rings[i].nodes.n);
        for(int j=0; j<rings[i].nodes.n; ++j){
            x[j] = nodes[rings[i].nodes[j]].crd[0];
            y[j] = nodes[rings[i].nodes[j]].crd[1];
        }
        VecF<double> com=periodicCentreOfMass(x,y);
        rings[i].setCrd(com);
    }

    //Write crystal
    if(!restart){
        for(int i=0; i<nodes.n; ++i) nodes[i].maximiseCnxs();
        writeCrds(prefix);
        writeNetwork(prefix,-1);
    }

    //Make initial procrystalline lattice
    if(!crystal) initialiseProLattice();
    //writeNetwork("./output/mc",-1);
}


void Lattice::initialiseSqLattice(int dim, bool restart) {
    //Initialise square lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 4;
    int dimSq = dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,dimSq+(latticeCnd-nodeCnd)*dimSq);
    for(int i=0; i<dimSq; ++i) rings.addValue(Ring(4*dim));
    periodicity = dim;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim; ++y) {
            crd[1]=y;
            for(int x=0; x<dim; ++x) {
                crd[0]=x;
                nodes[id].setCrd(crd);
                ++id;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx = y * dim + (id + dim - 1) % dim;
                nodes[id].addCnx(cnx);
                cnx = (id + dim) % dimSq;
                nodes[id].addCnx(cnx);
                cnx = y * dim + (id + 1) % dim;
                nodes[id].addCnx(cnx);
                cnx = (id + dimSq - dim) % dimSq;
                nodes[id].addCnx(cnx);
                id += 1;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            ring = id;
            nodes[id].addRing(ring);
            ring = y * dim + (id + 1) % dim;
            nodes[id].addRing(ring);
            ring=((y * dim + (id + 1) % dim) + dimSq - dim) % dimSq;
            nodes[id].addRing(ring);
            ring=(id + dimSq - dim) % dimSq;
            nodes[id].addRing(ring);
            id += 1;
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseTriLattice(int dim, bool restart) {
    //Initialise triangular lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 6;
    int dimSq = dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,2*dimSq+(latticeCnd-nodeCnd)*2*dimSq);
    for(int i=0; i<2*dimSq; ++i) rings.addValue(Ring(2*4*dim));
    periodicity[0] = dim;
    periodicity[1] = dim*sqrt(3)/2;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        double dy=sqrt(3.0)*0.5;
        for(int y=0; y<dim; ++y) {
            crd[1]=y*dy;
            for(int x=0; x<dim; ++x) {
                crd[0]=0.5*(y%2)+x;
                nodes[id].setCrd(crd);
                ++id;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx=y*dim+(id+dim-1)%dim;
                nodes[id].addCnx(cnx);
                if(y%2==0){
                    cnx=(cnx+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(id+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else{
                    cnx=(id+dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=((y+1)*dim+(id+1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                cnx=y*dim+(id+1)%dim;
                nodes[id].addCnx(cnx);
                if(y%2==0){
                    cnx=(id+dimSq-dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dim+(id+dim-1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else{
                    cnx=(dimSq+(y-1)*dim+(id+dim+1)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dim+(id+dim)%dim)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    int dimSq2=2*dimSq;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dim; ++x){
            ring=(2*y*dim+(id+dim-1)%dim);
            nodes[id].addRing(ring);
            ring=((2*y+1)*dim+(id)%dim);
            nodes[id].addRing(ring);
            ring=(2*y*dim+id%dim);
            nodes[id].addRing(ring);
            if(y%2==0){
                ring=(2*y*dim+id%dim-dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+(id+dim-1)%dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+(id+dim-1)%dim+dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
            }
            else{
                ring=(2*y*dim+(id+1)%dim-dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+id%dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
                ring=(2*(y-1)*dim+id%dim+dim+dimSq2)%dimSq2;
                nodes[id].addRing(ring);
            }
            ++id;
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseSnubSqLattice(int dim, bool restart) {
    //Initialise snub-square lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 5;
    int xDim=dim,yDim=dim*2;
    int dimSq4=4*xDim*yDim;
    int dim4=4*xDim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq4);
        for(int i=0; i<dimSq4; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,3*dimSq4/2+(latticeCnd-nodeCnd)*3*dimSq4/2);
    for(int i=0; i<3*dimSq4/2; ++i) rings.addValue(Ring(2*4*dim));
    double dxa=sqrt(3.0)/2;
    double dxe=1.0;
    double dya=0.5;
    double dye=0.5+sqrt(3.0)/2;
    periodicity[0] = xDim*(dxa*2+dxe);
    periodicity[1] = yDim*dye;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<yDim; ++y){
            crd[1]=dye*(y+0.5);
            for(int x=0; x<xDim; ++x){
                crd[0]=x*(2*dxa+dxe)+(y%2)*(dxe/2+dxa);
                nodes[y*dim4+x*4+0].setCrd(crd);
                crd[0]+=dxa;
                crd[1]+=dya;
                nodes[y*dim4+x*4+1].setCrd(crd);
                crd[1]-=2*dya;
                nodes[y*dim4+x*4+2].setCrd(crd);
                crd[0]+=dxa;
                crd[1]+=dya;
                nodes[y*dim4+x*4+3].setCrd(crd);
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<yDim; ++y){
            for(int x=0; x<xDim; ++x){
                //environment 0
                cnx=y*dim4+(id+dim4-1)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-2)+dim4)%dim4;
                else cnx=(id+dim4+2)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                cnx=id+2;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-3)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=id-dim4+1;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 1
                cnx=id-1;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-2)+dim4)%dim4;
                else cnx=(id+dim4+2)%dimSq4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-1)+dim4)%dim4;
                else cnx=((y+1)*dim4+(id+3)%dim4)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=id+2;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 2
                cnx=id-2;
                nodes[id].addCnx(cnx);
                cnx=id-1;
                nodes[id].addCnx(cnx);
                cnx=id+1;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-2)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+2)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-3)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+1)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                ++id;
                //environment 3
                cnx=id-1;
                nodes[id].addCnx(cnx);
                cnx=id-2;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(y+1)*dim4+((id%dim4-1)+dim4)%dim4;
                else cnx=((y+1)*dim4+((id%dim4+3)+dim4)%dim4)%dimSq4;
                nodes[id].addCnx(cnx);
                cnx=y*dim4+(id+dim4+1)%dim4;
                nodes[id].addCnx(cnx);
                if(y%2==0) cnx=(((y-1)*dim4+((id%dim4-2)+dim4)%dim4)+dimSq4)%dimSq4;
                else cnx=(y-1)*dim4+((id%dim4+2)+dim4)%dim4;
                nodes[id].addCnx(cnx);
                ++id;
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    map<string,int> ringCodes;
    int id0,id1,id2,id3;
    int ringId=0;
    for(int i=0; i<nodes.n; ++i){
        id0=i;
        for(int j=0; j<5; ++j){
            id1=nodes[i].cnxs[j];
            id2=nodes[i].cnxs[(j+1)%5];
            VecR<int> ringPath(0,4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if(!vContains(nodes[id1].cnxs,id2)){
                VecR<int> common=vCommonValues(nodes[id1].cnxs,nodes[id2].cnxs);
                common.delValue(id0);
                id3=common[0];
                ringPath.addValue(id3);
            }
            ringPath=vSort(ringPath);
            string rCode="";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            if(ringCodes.count(rCode)==0){
                ringCodes[rCode]=ringId;
                ++ringId;
            }
        }
    }
    for(int i=0; i<nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < 5; ++j) {
            id1 = nodes[i].cnxs[j];
            id2 = nodes[i].cnxs[(j + 1) % 5];
            VecR<int> ringPath(0, 4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if (!vContains(nodes[id1].cnxs, id2)) {
                VecR<int> common = vCommonValues(nodes[id1].cnxs, nodes[id2].cnxs);
                common.delValue(id0);
                id3 = common[0];
                ringPath.addValue(id3);
            }
            ringPath = vSort(ringPath);
            string rCode= "";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            ringId=ringCodes.at(rCode);
            nodes[i].addRing(ringId);
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseIsoSnubQuadLattice(int dim, bool restart) {
    //Initialise isosnub-quadrilateral lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 5;
    int dimSq=dim*dim;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,3*dimSq/2+(latticeCnd-nodeCnd)*3*dimSq/2);
    for(int i=0; i<3*dimSq/2; ++i) rings.addValue(Ring(2*4*dim));
    periodicity[0] = dim;
    periodicity[1] = dim/2+dim*sqrt(3)/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                crd[0]=x;
                if((y%4)>1) crd[0]+=0.5;
                nodes[id].setCrd(crd);
                ++id;
            }
            if(y%2==0) crd[1]+=1;
            else crd[1]+=sqrt(3)/2;
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim; ++y){
            for(int x=0; x<dim; ++x){
                cnx = y * dim + (id + dim - 1) % dim;
                nodes[id].addCnx(cnx);
                if(y%4==1){
                    cnx = ((id + dimSq + dim - 1) % dimSq);
                    if(x==0) cnx+=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = (id + dim) % dimSq;
                nodes[id].addCnx(cnx);
                if(y%4==3){
                    cnx = (id + dimSq + dim + 1) % dimSq;
                    if(x==dim-1) cnx-=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = y * dim + (id + 1) % dim;
                nodes[id].addCnx(cnx);
                if(y%4==2){
                    cnx = id - dim + 1;
                    if(x==dim-1) cnx-=dim;
                    nodes[id].addCnx(cnx);
                }
                cnx = (id + dimSq - dim) % dimSq;
                nodes[id].addCnx(cnx);
                if(y%4==0){
                    cnx = ((id + dimSq - dim - 1) % dimSq);
                    if(x==0) cnx+=dim;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    map<string,int> ringCodes;
    int id0,id1,id2,id3;
    int ringId=0;
    for(int i=0; i<nodes.n; ++i){
        id0=i;
        for(int j=0; j<5; ++j){
            id1=nodes[i].cnxs[j];
            id2=nodes[i].cnxs[(j+1)%5];
            VecR<int> ringPath(0,4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if(!vContains(nodes[id1].cnxs,id2)){
                VecR<int> common=vCommonValues(nodes[id1].cnxs,nodes[id2].cnxs);
                common.delValue(id0);
                id3=common[0];
                ringPath.addValue(id3);
            }
            ringPath=vSort(ringPath);
            string rCode="";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            if(ringCodes.count(rCode)==0){
                ringCodes[rCode]=ringId;
                ++ringId;
            }
        }
    }
    for(int i=0; i<nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < 5; ++j) {
            id1 = nodes[i].cnxs[j];
            id2 = nodes[i].cnxs[(j + 1) % 5];
            VecR<int> ringPath(0, 4);
            ringPath.addValue(i);
            ringPath.addValue(id1);
            ringPath.addValue(id2);
            if (!vContains(nodes[id1].cnxs, id2)) {
                VecR<int> common = vCommonValues(nodes[id1].cnxs, nodes[id2].cnxs);
                common.delValue(id0);
                id3 = common[0];
                ringPath.addValue(id3);
            }
            ringPath = vSort(ringPath);
            string rCode= "";
            for(int j=0; j<ringPath.n; ++j) rCode += "#" + to_string(ringPath[j]);
            ringId=ringCodes.at(rCode);
            nodes[i].addRing(ringId);
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseTriHexLattice(int dim, bool restart) {
    //Initialise trihex lattice, either from scratch or restarting previous run

    //Initialise nodes and rings
    latticeCnd = 4;
    int dim3=3*dim;
    int dim_2=dim/2;
    int dimSq3_4=dim*dim3/4;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq3_4);
        for(int i=0; i<dimSq3_4; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,dimSq3_4+(latticeCnd-nodeCnd)*dimSq3_4);
    for(int i=0; i<dimSq3_4; ++i) rings.addValue(Ring(4*dim));
    periodicity[0] = dim;
    periodicity[1] = 2*sqrt(3)*dim/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        double dy=sqrt(3)/2;
        VecF<double> crd(2);
        int id=0;
        for(int y=0; y<dim/4; ++y){
            for(int yy=0; yy<4; ++yy) {
                if (yy % 4 == 0 || yy % 4 == 2) {
                    for (int x = 0; x < dim; ++x) {
                        crd[0] = x;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                else if (yy % 4 == 1) {
                    for (int x = 0; x < dim_2; ++x) {
                        crd[0] = 2*x+0.5;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                else if (yy % 4 == 3) {
                    for (int x = 0; x < dim_2; ++x) {
                        crd[0] = 2*x+1.5;
                        nodes[id].setCrd(crd);
                        ++id;
                    }
                }
                crd[1] += dy;
            }
        }

        //Make node connections
        id=0;
        int cnx;
        for(int y=0; y<dim/4; ++y){
            for(int yy=0; yy<4; ++yy){
                if(yy%4==0){
                    for(int x=0; x<dim; ++x){
                        cnx = y*dim3+(x+dim-1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+int(floor(x/2.0));
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+(x+1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = (y*dim3-dim_2+(dim+int(ceil(x/2.0))-1)%dim_2+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==1){
                    for(int x=0; x<dim_2; ++x){
                        cnx = y*dim3+dim+dim_2+2*x;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+2*x+1;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+2*x+1;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+2*x;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==2){
                    for(int x=0; x<dim; ++x){
                        cnx = y*dim3+dim+dim_2+(x+dim-1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+dim+(int(ceil(x/2.0))-1+dim_2)%dim_2;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+dim_2+(x+1)%dim;
                        nodes[id].addCnx(cnx);
                        cnx = y*dim3+dim+int(floor(x/2.0));
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
                if(yy%4==3){
                    for(int x=0; x<dim_2; ++x){
                        cnx = ((y+1)*dim3+2*x+1)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3+(2*x+2+dim)%dim)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3-dim_2-dim+(2*x+2+dim)%dim+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        cnx = ((y+1)*dim3-dim_2-dim+(2*x+1+dim)%dim+dimSq3_4)%dimSq3_4;
                        nodes[id].addCnx(cnx);
                        ++id;
                    }
                }
            }
        }
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate rings to nodes
    int id=0;
    int ring;
    for(int y=0; y<dim/4; ++y){
        for(int yy=0; yy<4; ++yy){
            if(yy%4==0){
                for(int x=0; x<dim; ++x){
                    ring = y*dim3+(x+dim-1)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(int(ceil(x/2.0)))%dim_2;
                    nodes[id].addRing(ring);
                    ring = y*dim3+(x+0)%dim;
                    nodes[id].addRing(ring);
                    ring = (y*dim3-dim_2+int(floor(x/2.0))+dimSq3_4)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==1){
                for(int x=0; x<dim_2; ++x){
                    ring = y*dim3+dim+x;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+2*x;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(dim_2+x+1)%dim_2;
                    nodes[id].addRing(ring);
                    ring = y*dim3+2*x;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==2){
                for(int x=0; x<dim; ++x){
                    ring = y*dim3+dim+dim_2+(x+dim-1)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+dim+int(floor(x/2.0));
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+dim_2+(x+0)%dim;
                    nodes[id].addRing(ring);
                    ring = y*dim3+dim+(int(ceil(x/2.0)))%dim_2;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
            if(yy%4==3){
                for(int x=0; x<dim_2; ++x){
                    ring = (y+1)*dim3-dim_2+x;
                    nodes[id].addRing(ring);
                    ring = ((y+1)*dim3+(2*x+1+dim)%dim)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ring = (y+1)*dim3-dim_2+(dim_2+x+1)%dim_2;
                    nodes[id].addRing(ring);
                    ring = ((y+1)*dim3-dim_2-dim+(2*x+1+dim)%dim+dimSq3_4)%dimSq3_4;
                    nodes[id].addRing(ring);
                    ++id;
                }
            }
        }
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseHexLattice(int dim, bool restart) {
    //Initialise hexagonal lattice

    //Initialise nodes and rings
    latticeCnd = 3;
    int dimX = dim/2;
    int dimSq = dim*dimX;
    if(restart) {
        for(int i=0; i<nodes.n; ++i) nodes[i].resetRings();
    }
    else {
        nodes = VecR<Node>(0,dimSq);
        for(int i=0; i<dimSq; ++i) nodes.addValue(Node(latticeCnd));
    }
    Ring::autoId = 0;
    rings = VecR<Ring>(0,2*dimSq+(latticeCnd-nodeCnd)*2*dimSq);
    for(int i=0; i<dimSq/2; ++i) rings.addValue(Ring(dim));
    periodicity[0] = dimX*sqrt(3);
    periodicity[1] = 3*dim/4;

    //Make coordinates and allowed node connections if not previously
    if(!restart) {
        //Make coordinates
        VecF<double> crd(2);
        int id = 0;
        double dx = sqrt(3.0);
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dimX; ++x) {
                crd[0] += dx;
                nodes[id].setCrd(crd);
                ++id;
            }
            if (y % 4 == 0) {
                crd[0] = dx / 2;
                crd[1] += 0.5;
            } else if (y % 4 == 1) {
                crd[0] = dx / 2;
                crd[1] += 1.0;
            } else if (y % 4 == 2) {
                crd[0] = 0;
                crd[1] += 0.5;
            } else {
                crd[0] = 0;
                crd[1] += 1.0;
            }
        }

        //Make node connections
        id = 0;
        int cnx;
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dimX; ++x) {
                if (y % 4 == 0) {
                    cnx = (y + 1) * dimX + (x + dimX - 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y+1)*dimX+x;
                    nodes[id].addCnx(cnx);
                    cnx=(dimSq+(y-1)*dimX+x)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                else if (y % 4 == 1) {
                    cnx = (y - 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                    cnx = (y - 1) * dimX + (x + 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y + 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                }
                else if (y % 4 == 2) {
                    cnx = (y + 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                    cnx = (y + 1) * dimX + (x + 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y - 1) * dimX + x;
                    nodes[id].addCnx(cnx);
                }
                else {
                    cnx = (y - 1) * dimX + (x + dimX - 1) % dimX;
                    nodes[id].addCnx(cnx);
                    cnx = (y-1)*dimX+x;
                    nodes[id].addCnx(cnx);
                    cnx=((y+1)*dimX+x)%dimSq;
                    nodes[id].addCnx(cnx);
                }
                ++id;
            }
        }
    }

    //Associate rings to nodes
    int id=0,ring;
    int dimR=dimX;
    int dimRSq=dimR*dimR;
    int set4=0;
    for(int y=0; y<dim; ++y){
        for(int x=0; x<dimX; ++x){
            if(y%4==0){
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+(x+dimR-1)%dimR;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+x;
                nodes[id].addRing(ring);
            }
            else if(y%4==1){
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+(x+1)%dimR;
                nodes[id].addRing(ring);
                ring=((set4+dimR/2-1)%(dimR/2))*dimR*2+dimR+x;
                nodes[id].addRing(ring);
            }
            else if(y%4==2){
                ring=set4*dimR*2+dimR+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+(x+1)%dimR;
                nodes[id].addRing(ring);
            }
            else{
                ring=set4*dimR*2+dimR+(x+dimR-1)%dimR;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+dimR+x;
                nodes[id].addRing(ring);
                ring=set4*dimR*2+x;
                nodes[id].addRing(ring);
            }
            ++id;
        }
        if(y%4==3) ++set4;
    }

    //Make all node connections
    for(int i=0; i<nodes.n; ++i){
        nodes[i].maximiseCnxs();
    }

    //Associate nodes to rings
    for(int i=0; i<nodes.n; ++i){
        for(int j=0; j<nodes[i].rings.n; ++j){
            rings[nodes[i].rings[j]].addNode(i);
        }
    }
    //Order
    for(int i=0; i<rings.n; ++i){
        VecR<int> unordered=rings[i].nodes;
        VecR<int> ordered(0,unordered.n);
        ordered.addValue(unordered[0]);
        for(int j=1; j<unordered.n; ++j){
            VecR<int> shared=nodes[ordered[j-1]].getSharedCnxs(rings[i]);
            if(!vContains(ordered,shared[0])) ordered.addValue(shared[0]);
            else ordered.addValue(shared[1]);
        }
        rings[i].nodes=ordered;
    }
}


void Lattice::initialiseProLattice() {
    //Remove random connections from crystal lattice to make procrystalline

    //Write starting crystalline structure
    for(int i=0; i<nodes.n; ++i) nodes[i].maximiseCnxs();

    if(!pattern) for(int i=0; i<nodes.n; ++i) nodes[i].randomCnxs(nodeCnd,mtGen);
    else{
        if(randRL(mtGen)==0) for(int i=0; i<nodes.n; ++i) nodes[i].randomCnxs(nodeCnd,patternR,mtGen);
        else for(int i=0; i<nodes.n; ++i) nodes[i].randomCnxs(nodeCnd,patternL,mtGen);
    }

}


VecF<int> Lattice::generate(int maxIterations, double temperature) {
    //Monte Carlo search for fully coordinated lattice

    //Get initial energy as number of unsatisfied connections
    int energy=nodes.n*nodeCnd;
    for(int i=0; i<nodes.n; ++i){
        energy -= nodes[i].getReciprocalCnxs(nodes).n;
    }

    //Monte Carlo optimisation
    int iterations = 0;
    int cutoffIterations = 0;
    int cutoff = nodes.n*10;
    converged = false;
    uniform_int_distribution<int> randNode(0,nodes.n-1);

    for(;;) {
        int id=randNode(mtGen);
        int currCnd=nodes[id].getReciprocalCnxs(nodes).n;
        VecR<int> currCnxs=nodes[id].cnxs;
        if(!pattern) nodes[id].randomCnxs(nodeCnd,mtGen);
        else{
            if(randRL(mtGen)==0) nodes[id].randomCnxs(nodeCnd,patternR,mtGen);
            else nodes[id].randomCnxs(nodeCnd,patternL,mtGen);
        }
        int trialCnd=nodes[id].getReciprocalCnxs(nodes).n;
        int deltaCnd=trialCnd-currCnd;
        if(deltaCnd<0){
            if(rand01(mtGen)<exp(deltaCnd/temperature)){
                energy-=2*deltaCnd;
                cutoffIterations=0;
            }
            else{
                nodes[id].setCnxs(currCnxs);
                cutoffIterations+=1;
            }
        }
        else{
            energy-=2*deltaCnd;
            cutoffIterations=0;
        }
        if(energy==0){
            converged=true;
            break;
        }
        else if(cutoffIterations==cutoff) break;
        else if(iterations==maxIterations) break;
//        if(iterations%nodes.n*nodes.n==0){
//            cout<<energy<<endl;
//        }
        ++iterations;
//        if(iterations%1==0) writeNetwork("./vis",iterations);
    }

    //Set up status
    VecF<int> opt(2);
    opt[0]=converged;
    opt[1]=iterations;

    if(converged){
        if(calcRings) {
            int ringStatus = findRings();
            if (ringStatus == 1) opt[0] = 2;
            else if (ringStatus == 2) opt[0] = 3;
            if (opt[0] == 1) if (calcEnvs) findEnvironments();
        }
        if(calcChains){
            findChains();
        }
    }
//    writeNetwork("./output/mc",iterations);

    return opt;
}


int Lattice::generate(VecR<int> &pairA, VecR<int> &pairB) {
    //Generate from defined list of edges to break

    //Loop over edges and remove connections
    for(int i=0; i<pairA.n; ++i){
        int idA=pairA[i];
        int idB=pairB[i];
        nodes[idA].breakCnx(idB);
        nodes[idB].breakCnx(idA);
    }

    //Find rings
    findRings();

    //Find environments
    if(calcEnvs) findEnvironments();

    //Calculate energy
    int energy=nodes.n*nodeCnd;
    for(int i=0; i<nodes.n; ++i){
        energy -= nodes[i].getReciprocalCnxs(nodes).n;
    }

    return energy;
}


int Lattice::findRings() {

    //Merge rings
    for(int i=0; i<nodes.n; ++i){
        //Get broken connections
        VecR<int> vacantCnxs=nodes[i].getVacantCnxs();
        for(int j=0; j<vacantCnxs.n; ++j){
            if(i<vacantCnxs[j]){//prevent double counting
                 //Get ids of rings that nodes share
                 VecR<int> shared=nodes[i].getSharedRings(nodes[vacantCnxs[j]]);
                 VecR<int> rIds(0,shared.n);
                 for(int k=0; k<shared.n; ++k){
                     if(rings[shared[k]].checkEdge(i,vacantCnxs[j])) rIds.addValue(shared[k]);
                 }
                 if(rIds.n!=2){
//                     writeCrds("./split");
//                     writeNetwork("./split",i);
                     return 1; //sample split in two so fails
                 }
                 int rid0=rIds[0];
                 int rid1=rIds[1];
                 //Generate new ring by merging shared rings
                 Ring ring=rings[rid0].merge(rings[rid1],i,vacantCnxs[j]);
                 int rId2=ring.id;
                 //Deactivate old rings and change ids of associated nodes
                 rings[rid0].clear();
                 rings[rid1].clear();
                 for(int k=0; k<ring.nodes.n; ++k){
                     nodes[ring.nodes[k]].swapRing(rid0,rId2);
                     nodes[ring.nodes[k]].swapRing(rid1,rId2);
                 }
                 rings.addValue(ring);
            }
        }
    }

    //Find ring adjacencies
    int nid0,nid1;
    for(int i=0; i<nodes.n; ++i){
        nid0=i;
        for(int j=0; j<nodes[i].cnxs.n; ++j){
            nid1=nodes[i].cnxs[j];
            if(nid0<nid1){
                VecR<int> shared=nodes[nid0].getSharedRings(nodes[nid1]);
                VecR<int> rIds(0,shared.n);
                for(int k=0; k<shared.n; ++k){
                    if(rings[shared[k]].checkEdge(nid0,nid1)) rIds.addValue(shared[k]);
                }
                if(rIds.n==1){//self interaction
                    rings[rIds[0]].addCnx(rIds[0]);
                    rings[rIds[0]].addCnx(rIds[0]);
                }
                else if(rIds.n==2){
                    rings[rIds[0]].addCnx(rIds[1]);
                    rings[rIds[1]].addCnx(rIds[0]);
                }
                else{//should not occur
                    return 2;
                }
            }
        }
    }

    //Check number of adjacencies match number of nodes
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            if(rings[i].cnxs.n!=rings[i].nodes.n) return 2;
        }
    }

    //Calculate coordinates
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            VecF<double> x(rings[i].nodes.n),y(rings[i].nodes.n);
            for(int j=0; j<rings[i].nodes.n; ++j){
                x[j] = nodes[rings[i].nodes[j]].crd[0];
                y[j] = nodes[rings[i].nodes[j]].crd[1];
            }
            VecF<double> com=periodicCentreOfMass(x,y);
            //rings[i].setCrd(com);
        }
    }

    return 0;
}


int Lattice::findChains() {
    //Find chains of 2-coordinate nodes

    //Initialise chain vector
    chains=VecR<Chain>(0,nodes.n/2+1);

    //Find chains by following node connections until reach original node
    VecF<bool> nodeInc(nodes.n);
    nodeInc=false;
    for(int i=0; i<nodes.n; ++i){
        if(!nodeInc[i]){
            Chain chain(nodes.n);
            chain.addNode(i);
            int nId=i;
            int nIdPrev=nodes[nId].cnxs[0];
            for(;;){
                int nIdNext;
                if(nodes[nId].cnxs[0]==nIdPrev) nIdNext=nodes[nId].cnxs[1];
                else nIdNext=nodes[nId].cnxs[0];
                if(nIdNext==i) break;
                else{
                    chain.addNode(nIdNext);
                    nIdPrev=nId;
                    nId=nIdNext;
                }
            }
            for(int j=0; j<chain.nodes.n; ++j) nodeInc[chain.nodes[j]]=true;
            chains.addValue(chain);
        }
    }

    return 0;
}


void Lattice::chainAnalysis(VecR<int> &chainLengths) {
    //Analyse chains

    //Get chain lengths
    for(int i=0; i<chains.n; ++i){
        chainLengths.addValue(chains[i].nodes.n);
    }
}


VecF<double> Lattice::periodicCentreOfMass(VecF<double> &x, VecF<double> &y) {
    //Find centre of mass accounting for periodic boundary conditions

    VecF<double> thetaX,thetaY;
    thetaX = x*(2*M_PI/periodicity[0]);
    thetaY = y*(2*M_PI/periodicity[1]);
    VecF<double> alphaX,betaX,alphaY,betaY;
    alphaX = vCos(thetaX);
    betaX = vSin(thetaX);
    alphaY = vCos(thetaY);
    betaY = vSin(thetaY);
    double avAlphaX=vMean(alphaX);
    double avBetaX=vMean(betaX);
    double avAlphaY=vMean(alphaY);
    double avBetaY=vMean(betaY);
    double phiX = atan2(-avBetaX,-avAlphaX)+M_PI;
    double phiY = atan2(-avBetaY,-avAlphaY)+M_PI;
    VecF<double> com(2);
    com[0]=phiX*periodicity[0]/(2*M_PI);
    com[1]=phiY*periodicity[1]/(2*M_PI);

    return com;
}


void Lattice::findEnvironments() {
    //Assign environment type to each original ring in lattice


    if(latticeCode=="sq3"){
        //Square 3-coordinate lattice
        int dim=sqrt(nodes.n);
        //Find environments
        for(int i=0; i<nodes.n; ++i){
            //Active means all sides complete
            if(rings[i].active) rings[i].environment=0;
            else{
                //Otherwise determine missing sides
                VecF<int> ids(4);
                int row=floor(i/dim);
                ids[0]=row*dim+(i+dim-1)%dim;
                ids[1]=((row+1)*dim+(i+dim-1)%dim)%nodes.n;
                ids[2]=((row+1)*dim+i%dim)%nodes.n;
                ids[3]=row*dim+i%dim;
                VecF<int> vacancies(4);
                for(int j=0,k=1; j<4; ++j, ++k){
                    int id0=ids[j];
                    int id1=ids[k%4];
                    vacancies[j]=!vContains(nodes[id0].cnxs,id1);
                }
                //Single vacancy
                int nVac=vSum(vacancies);
                if(nVac==1){
                    for(int j=0; j<4; ++j){
                        if(vacancies[j]==1) rings[i].environment=j+1;
                    }
                }
                else if(nVac==2){
                    if(vacancies[0]==1 && vacancies[2]==1) rings[i].environment=5;
                    if(vacancies[1]==1 && vacancies[3]==1) rings[i].environment=6;
                }
            }
        }
    }
}

VecF<double> Lattice::networkAnalysis(VecF<int> &k, VecF<int> &pk, VecF< VecF<int> > &ejk) {
    //Calculate ring statistics and assortativity

    //Calculate ring statistics and adjacency matrix
    pk=VecF<int>(k.n);
    ejk=VecF< VecF<int> >(k.n);
    for(int i=0; i<k.n; ++i) ejk[i]=VecF<int>(k.n);
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            int ki=rings[i].getNumCnxs();
            for(int j=0; j<rings[i].cnxs.n; ++j){
                int kj=rings[rings[i].cnxs[j]].getNumCnxs();
                ++ejk[ki][kj];
            }
            ++pk[ki];
        }
    }

    VecF<double> res=calculateNetworkProperties(k,pk,ejk);
    return res;
}


VecF<double> Lattice::calculateNetworkProperties(VecF<int> &k, VecF<int> &pk, VecF<VecF<int> > &ejk) {
    //Calculate ring statistics, assortativity and Aboav-Weaire parameters

    //Normalisation
    double pNorm=vSum(pk),eNorm=0.0;
    for(int i=0; i<k.n; ++i) eNorm+=vSum(ejk[i]);

    //Moments and variance
    double k1,k2,k3;
    k1=vSum(k*pk)/pNorm;
    k2=vSum(k*k*pk)/pNorm;
    k3=vSum(k*k*k*pk)/pNorm;
    double mu_2 = k2-k1*k1;

    //Assortativity
    double assort=0.0;
    for(int i=0; i<k.n; ++i){
        for(int j=0; j<k.n; ++j){
            assort+=i*j*ejk[i][j];
//            if(i>=3 && i<=12){
//                if(j>=3 && j<=12) cout<<i<<" "<<j<<" "<<ejk[i][j]<<endl;
//            }
        }
    }
    assort/=eNorm;
    assort=k1*k1*assort-k2*k2;
    assort/=k1*k3-k2*k2;

    //Aboav-Weaire
    VecF<double> qk(k.n),mk(k.n);
    for(int i=0; i<k.n; ++i){
        qk[i]=vSum(ejk[i]);
        mk[i]=vSum(k*ejk[i]);
    }
    VecR<double> awX(0,k.n),awY(0,k.n);
    for(int i=0; i<k.n; ++i){
        if(qk[i]>0){
            awX.addValue(k1*(k[i]-k1));
            awY.addValue(k[i]*mk[i]/qk[i]);
        }
    }
    VecR<double> aw=vLinearRegression(awX,awY);

    //Get principle ring size
    int kMain;
    if(nodeCnd==3) kMain=6;
    else if(nodeCnd==4) kMain=4;
    else if(nodeCnd==5) kMain=3;

    //Return summary
    VecF<double> summary(7);
    summary[0]=pk[kMain]/pNorm;
    summary[1]=k1;
    summary[2]=mu_2;
    summary[3]=assort;
    summary[4]=1.0-aw[0];
    summary[5]=aw[1]-k1*k1;
    summary[6]=aw[2];

    return summary;
}


VecF<int> Lattice::getEnvironments() {
    //Get ring environments

    //Get environment for each ring
    VecF<int> envs(nodes.n);
    for (int i = 0; i < nodes.n; ++i){
        envs[i] = rings[i].environment;
    }

    return envs;
}


string Lattice::getEnvironmentCode() {
    //Get ring environment code

    //Get environment for each ring
    string envCode="";
    for (int i = 0; i < nodes.n; ++i){
        envCode += to_string(rings[i].environment);
    }

    return envCode;
}


void Lattice::rdfAnalysis(VecF<double> &latticeRDF, VecF<double> &dualRDF, double delta) {
    //Calculate rdfs for lattice points and dual

    //Lattice coordinates
    VecF<double> xLat(nodes.n),yLat(nodes.n);
    for(int i=0; i<nodes.n; ++i){
        xLat[i]=nodes[i].crd[0];
        yLat[i]=nodes[i].crd[1];
    }
    calculateRDF(xLat,yLat,latticeRDF,delta);

    //Dual coordinates
    VecF<double> xDual(Ring::totalActive),yDual(Ring::totalActive);
    int j=0;
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            xDual[j]=rings[i].crd[0];
            yDual[j]=rings[i].crd[1];
            ++j;
        }
    }
    calculateRDF(xDual,yDual,dualRDF,delta);
}


void Lattice::calculateRDF(VecF<double> &x, VecF<double> &y, VecF<double> &rdf, double delta) {
    //Calculate rdf

    //Get periodic information
    VecF<double> pbc=periodicity,rpbc(2);
    VecF<double> mic=periodicity/2.0;
    rpbc[0]=1.0/periodicity[0];
    rpbc[1]=1.0/periodicity[1];
    int maxBin;
    double maxDist,maxDistSq;
    if(pbc[0]<=pbc[1]) maxDist=pbc[0]/2.0;
    else maxDist=pbc[1]/2.0;
    maxDistSq=maxDist*maxDist;
    maxBin=floor(maxDist/delta)+1;

    //Calculate pairwise distances and bin
    double xI,yI,b;
    double dx,dy,dSq,d;
    for (int i=0; i<x.n-1; ++i) {
        xI = x[i];
        yI = y[i];
        for (int j=i+1; j<x.n; ++j) {
            dx = xI - x[j];
            dy = yI - y[j];
            dx -= pbc[0] * nearbyint(dx * rpbc[0]);
            dy -= pbc[1] * nearbyint(dy * rpbc[1]);
            dSq = dx * dx + dy * dy;
            if (dSq < maxDistSq) {
                d = sqrt(dSq);
                b = floor(d / delta);
                rdf[b] += 2;
            }
        }
    }
}


void Lattice::skAnalysis(VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax) {
    //Calculate structure factors for lattice points and dual

    //Lattice coordinates
    VecF<double> xLat(nodes.n),yLat(nodes.n);
    for(int i=0; i<nodes.n; ++i){
        xLat[i]=nodes[i].crd[0];
        yLat[i]=nodes[i].crd[1];
    }
    VecF<double> latticeSk;
//    calculateSk(xLat,yLat,sk,multiplicity,delta,nMax);

    //Dual coordinates
    VecF<double> xDual(Ring::totalActive),yDual(Ring::totalActive);
    int j=0;
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active){
            xDual[j]=rings[i].crd[0];
            yDual[j]=rings[i].crd[1];
            ++j;
        }
    }
    calculateSk(xDual,yDual,sk,multiplicity,delta,nMax);

}


void Lattice::calculateSk(VecF<double> &x, VecF<double> &y, VecF<double> &sk, VecF<int> &multiplicity, double delta, int nMax) {
    //Calculate structure factor

    //Variables
    double kx=2*M_PI/periodicity[0];
    double ky=2*M_PI/periodicity[1];

    //Set up real and imaginary components
    int numN=nMax*2+1;
    VecF< VecF<double> > real(numN),imag(numN),skk(numN);
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        real[i]=VecF<double>(numN);
        imag[i]=VecF<double>(numN);
        skk[i]=VecF<double>(numN);
    }
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            VecF<double> v=x*nx*kx+y*ny*ky;
            real[i][j]=vSum(vCos(v));
            imag[i][j]=vSum(vSin(v));
        }
    }

    //Calculate structure factor
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            skk[i][j]=(real[i][j]*real[i][j]+imag[i][j]*imag[i][j])/x.n;
        }
    }

    //Normalise
    double kSq,rDelta=1.0/delta;
    int b;
    for(int nx=-nMax,i=0; nx<=nMax; ++nx,++i){
        for(int ny=-nMax,j=0; ny<=nMax; ++ny,++j) {
            kSq=kx*nx*kx*nx+ky*ny*ky*ny;
            cout<<nx<<" "<<ny<<" "<<skk[i][j]<<endl;
            if(kSq>0) {
                b=floor(sqrt(kSq)*rDelta);
                sk[b]+=skk[i][j];
                ++multiplicity[b];
            }
        }
    }
}


int Lattice::getNumNodes() {
    //Get number of nodes

    return nodes.n;
}


int Lattice::getNumActiveRings() {
    //Get number of active rings

    return Ring::totalActive;
}


VecF<double> Lattice::getPeriodicity() {
    //Get periodicity information

    return periodicity;
}


void Lattice::getEdgeCombinations(VecR<int> &pairA, VecR<int> &pairB, int &n, int &r) {
    //Get unique pairs of edges

    pairA=VecR<int>(0,nodes.n*nodeCnd);
    pairB=VecR<int>(0,nodes.n*nodeCnd);
    int idA,idB;
    for(int i=0; i<nodes.n; ++i){
        idA=i;
        for(int j=0; j<nodes[i].cnxs.n; ++j){
            idB=nodes[i].cnxs[j];
            if(idA<idB){
                pairA.addValue(idA);
                pairB.addValue(idB);
            };
        }
    }

    n = pairA.n;
    r = nodes.n*(latticeCnd-nodeCnd)/2;
}


void Lattice::writeCrds(string prefix) {
    //Write periodicity, node and ring coordinates

    //Set up output files
    OutputFile crdFile(prefix+"_crds.dat");
    OutputFile rcrdFile(prefix+"_rcrds.dat");
    crdFile.initVariables(6,3,60,12);
    rcrdFile.initVariables(6,3,60,12);

    //Write coordination
    crdFile.write(latticeCnd,nodeCnd);

    //Write periodicity
    crdFile.writeRowVector(periodicity);
    rcrdFile.writeRowVector(periodicity);

    //Write node coordinates
    for(int i=0; i<nodes.n; ++i) crdFile.writeRowVector(nodes[i].crd);

    //Write ring coordinates
    for(int i=0; i<rings.n; ++i) rcrdFile.writeRowVector(rings[i].crd);
}


void Lattice::writeNetwork(string prefix, int index) {
    //Write node connections and rings

    //Set up output file
    OutputFile sampleFile(prefix+"_sample_"+to_string(index)+".dat");

    //Write node connections 
    sampleFile.write(nodes.n);
    for(int i=0; i<nodes.n; ++i) sampleFile.writeRowVector(nodes[i].cnxs);
    //Write ring ids
    sampleFile.write(Ring::totalActive);
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.write(i);
    }
    //Write ring connections
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].cnxs);
    }
    //Write rings
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].nodes);
    }
    //Write ring centres of mass
    for(int i=0; i<rings.n; ++i){
        if(rings[i].active) sampleFile.writeRowVector(rings[i].crd);
    }
    //Write chains
    sampleFile.write(chains.n);
    for(int i=0; i<chains.n; ++i){
        sampleFile.writeRowVector(chains[i].nodes);
    }
    //Write ring environments
    if(calcEnvs) for(int i=0; i<nodes.n; ++i) sampleFile.write(rings[i].environment);
}
