import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import matplotlib.pylab as pylab
import numpy as np
import sys

np.random.seed(0)

class Plot:
    """
    Plot procrystalline lattice.
    """


    def __init__(self,nodes=False,cnxs=False,rings=False,periodic=False,dual=False,envs=False):
        """
        Initialise with plot options.
        """

        # Get options
        self.nodes = nodes
        self.cnxs = cnxs
        self.rings = rings
        self.periodic = periodic
        self.dual = dual
        self.envs = envs

        # Set up empty figure
        params = {"figure.figsize": (6, 6)}
        pylab.rcParams.update(params)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_axis_off()


    def __call__(self,prefix,sample,**kwargs):
        """
        Plot procrystalline lattice.
        """

        # Images to generate based on periodicity
        if self.periodic:
            images = [-1,0,1]
        else:
            images = [0]

        # Unpack options
        self.lw = kwargs.get("lw",1.0)
        self.lc = kwargs.get("lc","k")
        self.ms = kwargs.get("ms",10)
        self.mc = kwargs.get("mc","k")
        save = kwargs.get("save",False)

        # Load data
        self.load_sample(prefix,sample)

        # Add images to plot
        for y in images:
            for x in images:
                self.plot_rings(x,y)
                self.plot_nodes(x,y)
                self.plot_cnxs(x,y)
                self.plot_dual(x,y)
                self.plot_envs(x,y)

        # Display
        if save:
            plt.savefig("config.png",dpi=300)
            plt.savefig("config.pdf")
        plt.show()


    def load_sample(self,prefix,sample):
        """
        Load coordinates and network information.
        """

        # Node coordination, periodicity and node crds
        data = np.genfromtxt(prefix+"_crds.dat")
        self.cnd = data[0,:].astype(int)
        self.pbc = data[1,:]
        self.mic = self.pbc/2
        self.node_crds = data[2:,:]
        if self.cnd[1] == 2:
            self.mean_ring_size = 0
        elif self.cnd[1] == 3:
            self.mean_ring_size = 6
        elif self.cnd[1] == 4:
            self.mean_ring_size = 4
        elif self.cnd[1] == 5:
            self.mean_ring_size = 3+1/3.0
        print(self.mean_ring_size)

        # Reciprocal crds
        data = np.genfromtxt(prefix+"_rcrds.dat")
        self.rec_crds = data[1:,:]

        # Node connections and rings
        self.node_cnxs = []
        self.node_rings = []
        self.node_chains = []
        self.dual_crds = []
        self.dual_cnxs = []
        self.dual_ids = {}
        with open("{}_sample_{}.dat".format(prefix,sample),"r") as f:
            n = int(f.readline())
            for i in range(n):
                cnxs = f.readline().split()
                for j in cnxs:
                    self.node_cnxs.append([i,j])
            self.node_cnxs = np.array(self.node_cnxs,dtype=int)
            n = int(f.readline())
            for i in range(n):
                j = int(f.readline())
                self.dual_ids[j] = i
            for i in range(n):
                cnxs = f.readline().split()
                for j in cnxs:
                    self.dual_cnxs.append([int(i),self.dual_ids[int(j)]])
            for i in range(n):
                ring = f.readline().split()
                self.node_rings.append(np.array(ring,dtype=int))
            for i in range(n):
                self.dual_crds.append(np.array([float(x) for x in f.readline().split()]))
            self.dual_crds = np.array(self.dual_crds)
            # n = int(f.readline())
            # for i in range(n):
            #     chain = np.array(f.readline().split(),dtype=int)
            #     self.node_chains.append(chain)
            # if self.envs:
            #     self.rec_envs = np.zeros(self.rec_crds.shape[0],dtype=int)
            #     for i in range(self.rec_envs.size):
            #         self.rec_envs[i] = int(f.readline())

        # Make connections accounting for periodicity
        self.node_cnx_crds = np.zeros((self.node_cnxs.shape[0],4))
        crd_i = np.zeros(2)
        crd_j = np.zeros(2)
        for i,p in enumerate(self.node_cnxs):
            crd_i[:] = self.node_crds[p[0]][:]
            crd_j[:] = self.node_crds[p[1]][:]
            x = crd_j[0]-crd_i[0]
            y = crd_j[1]-crd_i[1]
            if x>self.mic[0]: x-=self.pbc[0]
            elif x<-self.mic[0]: x+=self.pbc[0]
            if y>self.mic[1]: y-=self.pbc[1]
            elif y<-self.mic[1]: y+=self.pbc[1]
            self.node_cnx_crds[i,0] = crd_i[0]
            self.node_cnx_crds[i,1] = crd_i[0]+x/2
            self.node_cnx_crds[i,2] = crd_i[1]
            self.node_cnx_crds[i,3] = crd_i[1]+y/2

        # Make rings accounting for periodicity
        self.ring_crds = []
        self.max_ring_size = 0
        for ring in self.node_rings:
            crds = np.zeros((ring.size,2))
            for i,j in enumerate(ring):
                crds[i,:] = self.node_crds[j,:]
            for i in range(1,ring.size):
                x = crds[i,0] - crds[i-1,0]
                y = crds[i,1] - crds[i-1,1]
                if x>self.mic[0]: x -= self.pbc[0]
                elif x<-self.mic[0]: x += self.pbc[0]
                if y>self.mic[1]: y -= self.pbc[1]
                elif y<-self.mic[1]: y += self.pbc[1]
                crds[i,0] = crds[i-1,0] + x
                crds[i,1] = crds[i-1,1] + y
            x_com = np.average(crds[:,0])
            y_com = np.average(crds[:,1])
            if x_com>self.mic[0]: crds[:,0]-=self.pbc[0]
            elif x_com<-self.mic[0]: crds[:,0]+=self.pbc[0]
            if y_com>self.mic[1]: crds[:,1]-=self.pbc[1]
            elif y_com<-self.mic[1]: crds[:,1]+=self.pbc[1]
            self.ring_crds.append(crds)
            if ring.size > self.max_ring_size:
                self.max_ring_size = ring.size
        self.init_ring_colours(self.mean_ring_size,self.max_ring_size)

        # Make chain colouring
        if len(self.node_chains)>0:
            self.chains = True
        else:
            self.chains = False
        self.init_chain_colours()


    def plot_nodes(self,x_shift,y_shift):
        """
        Plot image nodes as scatter. 
        """

        if not self.nodes: return # Bounce if option not selected

        if not self.chains:
            self.ax.scatter(self.node_crds[:,0]+x_shift*self.pbc[0],self.node_crds[:,1]+y_shift*self.pbc[1],
                        marker="o",s=self.ms,c=self.mc,zorder=1)
        else:
            for i,chain in enumerate(self.node_chains):
                self.ax.scatter(self.node_crds[chain,0]+x_shift*self.pbc[0],self.node_crds[chain,1]+y_shift*self.pbc[1],
                            marker="o",s=self.ms,color=self.chain_colours[i],edgecolor='k',zorder=1)
        # for r in self.ring_crds:
        #     self.ax.scatter(r[:,0]+x_shift*self.pbc[0],r[:,1]+y_shift*self.pbc[1],
        #                 marker="o",s=self.ms,c=self.mc,zorder=1)

        # for i,c in enumerate(self.node_crds):
        #     self.ax.text(c[0],c[1],i,size=8)


    def plot_cnxs(self,x_shift,y_shift):
        """
        Plot image connections as lines. 
        """

        if not self.cnxs: return # Bounce if option not selected

        # Find under coordinated
        cc_count = {}
        for c in self.node_cnx_crds:
            cc = [c[1],c[3]]
            if cc[0]<0:  cc[0] += self.pbc[0]
            if cc[1]<0:  cc[1] += self.pbc[1]
            cc = tuple(cc)
            if cc in cc_count:
                cc_count[cc] += 1
            else:
                cc_count[cc] = 1
        complete = np.zeros(self.node_cnx_crds[:,0].size,dtype=bool)
        for i,c in enumerate(self.node_cnx_crds):
            cc = [c[1],c[3]]
            if cc[0]<0:  cc[0] += self.pbc[0]
            if cc[1]<0:  cc[1] += self.pbc[1]
            cc = tuple(cc)
            if cc_count[cc]==2:
                complete[i]=True

        self.node_cnx_crds[:,:2] += x_shift*self.pbc[0]
        self.node_cnx_crds[:,2:] += y_shift*self.pbc[1]
        for cnx_crd in self.node_cnx_crds[complete]:
            self.ax.plot(cnx_crd[:2],cnx_crd[2:],c=self.lc,lw=self.lw,zorder=-1)
        for cnx_crd in self.node_cnx_crds[~complete]:
            self.ax.plot(cnx_crd[:2],cnx_crd[2:],c='r',lw=self.lw,zorder=-1)
        self.node_cnx_crds[:,:2] -= x_shift*self.pbc[0]
        self.node_cnx_crds[:,2:] -= y_shift*self.pbc[1]


    def plot_rings(self,x_shift,y_shift):
        """
        Plot rings as polygons.
        """

        if not self.rings: return # Bounce if option not selected

        patches = []
        colours = []
        for ring in self.ring_crds:
            ring[:,0] += x_shift*self.pbc[0]
            ring[:,1] += y_shift*self.pbc[1]
            # xbox = (ring[:,0]>0)*(ring[:,0]<6) # square
            # ybox = (ring[:,1]>0)*(ring[:,1]<6) # square
            # xbox = (ring[:,0]>1.35)*(ring[:,0]<7.25) # snub
            # ybox = (ring[:,1]>0.97)*(ring[:,1]<7.2) # snub
            # xbox = (ring[:,0]>1.35)*(ring[:,0]<6.3) # isosnub
            # ybox = (ring[:,1]>0.97)*(ring[:,1]<4.8) # isosnub
            # xbox = (ring[:,0]>1.35)*(ring[:,0]<6.71) # tri
            # ybox = (ring[:,1]>0.97)*(ring[:,1]<5.5) # tri
            # xbox = (ring[:,0]>0.84)*(ring[:,0]<8.4) # trihex
            # ybox = (ring[:,1]>0.56)*(ring[:,1]<8.1) # trihex
            xbox=np.ones_like(ring[:,0])
            ybox=np.ones_like(ring[:,1])
            if np.all(xbox*ybox==1):
                patches.append(Polygon(np.array(ring), True))
                colours.append(self.ring_colours[ring[:,0].size])
                # self.ax.scatter(ring[:,0],ring[:,1],c=self.mc,s=self.ms)
            ring[:,0]-=x_shift*self.pbc[0]
            ring[:,1]-=y_shift*self.pbc[1]
        self.ax.add_collection(PatchCollection(patches,facecolor=colours,linewidths=self.lw,edgecolor="k",zorder=0))


    def plot_dual(self,x_shift,y_shift):
        """
        Plot image ring com as scatter.
        """

        if not self.dual: return # Bounce if option not selected

        self.ax.scatter(self.dual_crds[:,0]+x_shift*self.pbc[0],self.dual_crds[:,1]+y_shift*self.pbc[1],
                        marker="s",s=self.ms,zorder=1,color='grey')
        # for cnx in self.dual_cnxs:
        #     print(cnx)
        #     c0 = self.dual_crds[cnx[0],:]
        #     c1 = self.dual_crds[cnx[1],:]
        #     v = c1-c0
        #     if v[0]>self.mic[0]: v[0]-=self.pbc[0]
        #     elif v[0]<-self.mic[0]: v[0]+=self.pbc[0]
        #     if v[1]>self.mic[1]: v[1]-=self.pbc[1]
        #     elif v[1]<-self.mic[1]: v[1]+=self.pbc[1]
        #     self.ax.plot([c0[0],c0[0]+v[0]],[c0[1],c0[1]+v[1]],color='grey',lw=self.lw)

        # for i,c in enumerate(self.dual_crds):
        #     self.ax.text(c[0],c[1],i,size=8)


    def plot_envs(self,x_shift,y_shift):
        """
        Plot image environments as scatter.
        """

        if not self.envs: return # Bounce if option not selected

        cmap = cm.get_cmap('Set1')
        colours = []
        for e in self.rec_envs:
            colours.append(cmap(e))
        self.ax.scatter(self.rec_crds[:,0]+x_shift*self.pbc[0],self.rec_crds[:,1]+y_shift*self.pbc[1],
                        marker="o",s=self.ms,facecolors=colours,edgecolors='k',zorder=1)


    def init_ring_colours(self,av_ring_size=6,max_ring_size=10):
        """
        Initialise colouring for rings.
        """
        av_ring_size=6
        map_lower = cm.get_cmap('Blues_r', 128)
        map_upper = cm.get_cmap('Reds', 128)
        map_mean=cm.get_cmap("Greys")
        map_lower=ListedColormap(map_lower(np.arange(20,100)))
        map_upper=ListedColormap(map_upper(np.arange(20,100)))

        norm_lower=Normalize(vmin=av_ring_size-3,vmax=av_ring_size)
        norm_upper=Normalize(vmin=av_ring_size,vmax=av_ring_size+6)
        colour_mean=map_mean(50)
        self.ring_colours=[]
        for i in range(max_ring_size+1):
            if i < 3:
                self.ring_colours.append("white")
            elif np.abs(i-av_ring_size)<1e-6:
                self.ring_colours.append(colour_mean)
            elif i<av_ring_size:
                self.ring_colours.append(map_lower(norm_lower(i)))
            else:
                self.ring_colours.append(map_upper(norm_upper(i)))
            # if i%2==0:
            #     self.ring_colours[-1] = 'whitesmoke'
            # else:
                # self.ring_colours[-1] = 'gold'
              # self.ring_colours[-1] = 'whitesmoke'


    def init_chain_colours(self):
        """
        Initialising colouring for chains
        """

        n_chains = len(self.node_chains)
        cmap=cm.get_cmap("rainbow")
        self.chain_colours=cmap(np.linspace(0,1,n_chains))
        # cmap=ListedColormap(cmap(np.arange(0,1,0.001)))
        # lengths=np.array([c.size for c in self.node_chains])
        # norm=Normalize(vmin=0,vmax=np.max(lengths))
        # rand = np.random.uniform(0,1,len(self.node_chains))
        # rand = np.ones(len(self.node_chains))
        # self.chain_colours = [cmap(norm(l)) for l in lengths]


if __name__ == "__main__":

    prefix = sys.argv[1]
    sample = int(sys.argv[2])

    if len(sys.argv) <= 3:
        nodes = True
        cnxs = False
        rings = True
        periodic = False
        dual = False
        envs = False
        save = False
    else:
        flags = sys.argv[3]
        nodes = 'n' in flags
        cnxs = 'c' in flags
        rings = 'r' in flags
        periodic = 'p' in flags
        dual = 'd' in flags
        envs = 'e' in flags
        save = 's' in flags

    plot=Plot(nodes=nodes,cnxs=cnxs,rings=rings,periodic=periodic,dual=dual,envs=envs)
    plot(prefix,sample,ms=0,lw=0.75,save=save)
