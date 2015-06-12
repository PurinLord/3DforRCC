
import sys
from matplotlib import pyplot as plt
from sklearn.manifold import Isomap
from scipy.spatial.distance import euclidean

def plotTrajectory(dfile):
    fin = open(dfile)

    Vsteps = []
    Vtarget = fin.readline().strip().split()
    Vtarget = map(float,Vtarget)
    Vsteps.append(Vtarget)
    for l in fin:
        l = l.strip().split()
        if len(l) != 26: continue
        l = map(float,l)
        Vsteps.append(l)


    distances = [euclidean(a,Vsteps[0]) for a in Vsteps[1:]]
    print len(distances)

    _map = plt.get_cmap("winter")
    distcolors = _map(distances)


    dimred = Isomap(n_components=2)
    Vsteps = dimred.fit_transform(Vsteps)



    #objective vector
    plt.scatter(Vsteps[0,0],Vsteps[0,1],color='red',s=30,marker=(5,1))
    #Optimization steps
    plt.scatter(Vsteps[1:,0],Vsteps[1:,1],color=distcolors,alpha=0.5)

    plt.show()

def main():
    fn = sys.argv[1]
    plotTrajectory(dfile=fn)

if __name__ == '__main__':
    main()
