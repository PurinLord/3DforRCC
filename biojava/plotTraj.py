
import sys
from matplotlib import pyplot as plt
from sklearn.manifold import Isomap
from scipy.spatial.distance import euclidean
import matplotlib.cm as cm
import numpy as np

def plotTrajectory(dfile):
    fin = open(dfile)

    fin.readline()
    fin.readline()
    fin.readline()
    Vsteps = []
    Vsolutions = []
    Vtarget = fin.readline().strip().split()
    Vtarget = map(float,Vtarget)
    Vsteps.append(Vtarget)
    solFlag = False
    for l in fin:
        l = l.strip().split()
        if len(l) == 2: solFlag = True; continue 
        if len(l) != 26: continue
        if solFlag:
				    Vsolutions.append(map(float,l))
				    solFlag = False
        l = map(float,l)
        Vsteps.append(l)

    longitud = len(Vsteps)
    #print  Vsolutions 
    Vsteps.extend(Vsolutions)
    distances = [euclidean(a,Vsteps[0]) for a in Vsteps[1:]]
    print len(distances)
    colors = cm.cool(np.linspace(0, 1, len(distances)))

    _map = plt.get_cmap("winter")
    distcolors = _map(distances)


    dimred = Isomap(n_components=2)
    Vsteps = dimred.fit_transform(Vsteps)



    #objective vector
    plt.scatter(Vsteps[0,0],Vsteps[0,1],color='red',s=30,marker=(5,1))
    #Optimization steps
    plt.scatter(Vsteps[1:,0],Vsteps[1:,1],color=colors,alpha=0.5)
		#
    plt.scatter(Vsteps[longitud:,0],Vsteps[longitud:,1],color='yellow',s=25)

    #plt.show()
    plt.savefig('gabrile_purin/plot.png', dpi=100)

def main():
    fn = sys.argv[1]
    plotTrajectory(dfile=fn)

if __name__ == '__main__':
    main()
