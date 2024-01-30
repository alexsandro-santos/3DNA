import numpy as np
from .Traj3D import * 
import time
from random import *
from copy import deepcopy

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

def exponentielle(x, y):
        return math.exp(x / y)

def voisins(table:RotTable, p: float):          # On a une table et on veut renvoyer une liste de 32 tables tel que l'on a modifié chaque coefficient proportionnellement à son écart type
        new_table_list = []
        table = table.getTable
        
        # lcouple contient la liste des couple possibles suivi de leur complémentaire : "AA" suivi de "TT". /!\ On inclut pas les couples qui sont leur propres complémentaires (AT, GC, GC, TA)
        lcouple = ["AA", "TT", "AC", "GT", "AG", "CT", "CA", "TG", "CC", "GG", "GA", "TC"]
        # lcouple2 contient les couples qui sont leur propres complémentaires
        lcouple2 = ["AT", "GC", "GC", "TA"]
        for i in range(0,lcouple,2):
                # A chaque couple (et son complémentaire) on fait +- son écart type pour le twist et aussi pour le wedge
                couple = lcouple[i]
                compl = lcouple[i+1]
                
                # Modif le twist 
                
                new_plus_Twist = deepcopy(table)
                new_plus_Twist.setTwist(couple, table[couple][0] + p*table[couple][3])  # table["AA"][3] écart type pour le twist de AA et # table["AA"][0] pour la valeur
                new_plus_Twist.setTwist(compl, table[compl][0] + p*table[compl][3])  # Implique changement du complémentaire  
                new_table_list.append(new_plus_Twist)

                new_moins_Twist = deepcopy(table)
                new_moins_Twist.setTwist(couple, table[couple][0] - p*table[couple][3]) 
                new_moins_Twist.setTwist(compl, table[compl][0] - p*table[compl][3])  
                new_table_list.append(new_moins_Twist)
                
                # Modif le wedge
                
                new_plus_Wedge = deepcopy(table)
                new_plus_Wedge.setWedge(couple, table[couple][1] + p*table[couple][4]) 
                new_plus_Wedge.setWedge(compl, table[compl][1] + p*table[compl][4])  
                new_table_list.append(new_plus_Wedge)
                
                new_moins_Wedge = deepcopy(table)
                new_moins_Wedge.setWedge(couple, table[couple][1] - p*table[couple][4]) 
                new_moins_Wedge.setWedge(compl, table[compl][1] - p*table[compl][4])  
                new_table_list.append(new_moins_Wedge)
        
        for i in range(len(lcouple2)):
                couple = lcouple2[i]
                
                new_plus_Twist = deepcopy(table)
                new_plus_Twist.setTwist(couple, table[couple][0] + p*table[couple][3])  
                new_table_list.append(new_plus_Twist)

                new_moins_Twist = deepcopy(table)
                new_moins_Twist.setTwist(couple, table[couple][0] - p*table[couple][3]) 
                new_table_list.append(new_moins_Twist)
                
                # Modif le wedge
                
                new_plus_Wedge = deepcopy(table)
                new_plus_Wedge.setWedge(couple, table[couple][1] + p*table[couple][4]) 
                new_table_list.append(new_plus_Wedge)
                
                new_moins_Wedge = deepcopy(table)
                new_moins_Wedge.setWedge(couple, table[couple][1] - p*table[couple][4]) 
                new_table_list.append(new_moins_Wedge)
        return new_table_list


        
        
def recuit_simule(seq):
        tps_init = time.clock()
        temps = 0
        s = RotTable()
        e = evaluation(s,seq)
        temperature = 1 

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

