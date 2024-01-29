import numpy as np
from .Traj3D import * 
import time
from random import uniform

def evaluation(traj):
        vec1 = traj.getTraj()[0]
        vec2 = traj.getTraj()[-1]
        diff = vec1 - vec2
        return diff.length
        #xyz = np.array(Traj3D.getTraj(traj))
        #x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        #return np.sqrt((x[0]-x[-1])**2 + (y[0]-y[-1])**2 + (z[0]-z[-1])**2)

def exponentielle(x, y):
    return math.exp(x / y)

def recuit_simule(seq):
        tps_init = time.clock()
        temps = 0
        s = RotTable()
        e = evaluation(s,seq)
        temperature = 1 

        while(temps < 100 and temperature > 0.01):
                nombre_aleatoire = uniform()
                for sn in voisin(s):
                        en = evaluation(s,seq)
                        if(en<e or nombre_aleatoire < exponentielle(en-e,temperature)):
                                s = sn
                                e = en
                temperature = 0.99 * temperature 
                tpsi = time.clock()
                temps = tpsi-tps_init     
        print(temperature)
        return s