import numpy as np
from .Traj3D import * 
import time
from random import uniform

def evaluation(traj): #On calcule la distance euclidienne entre le premier point et le dernier point, on cherchera donc à minimiser cette valeur 
        vec1 = traj.getTraj()[0]
        vec2 = traj.getTraj()[-1]
        diff = vec1 - vec2
        return diff.length

def test_evaluation():
    rot_table = RotTable()
    traj = Traj3D()
    
    # Read file
    lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
    # Formatting
    seq = ''.join(lineList[1:])
    
    traj.compute(seq, rot_table)
    print(evaluation(traj))
    assert(evaluation(traj)==2)

def exponentielle(x, y): #On caclule l'exponentielle de x sur y 
    return math.exp(x / y)

assert(exponentielle(0,5)==1)

def recuit_simule(trajectoire,table,sequence): 
        tps_init = time.clock() # On regarde le temps actuel en seconde
        temps = 0 #On initialise une variable temps qui permet de savoir quand notre algorithme se terminera 
        s = table #Il s'agit de la variable qui contiendra la table qu'on conserve
        e = evaluation(trajectoire) #L'énergie d'une trajectoire 
        temperature = 10 #Temperature initial

        while(temps < 100 and temperature > 0.1):
                nombre_aleatoire = uniform() #On choisi un nombre uniformement dans [0,1]
                for sn in voisin(s): #Pour tous les voisins de la table s, on calcule sa trajectoire ainsi que son énergie
                        trajectoire.compute(sequence,sn) 
                        en = evaluation(trajectoire)
                        if(en<e or nombre_aleatoire < exponentielle(en-e,temperature)):
                                s = sn
                                e = en
                temperature = 0.99 * temperature #On change la temperature
                tpsi = time.clock() 
                temps = tpsi-tps_init #On actualise le temps    
        print(temperature)
        trajectoire = trajectoire.compute(s)
        return trajectoire

