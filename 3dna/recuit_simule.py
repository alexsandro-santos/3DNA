import numpy as np
from .Traj3D import * 
import time
from random import *
from copy import deepcopy

def evaluation(traj): #On calcule la distance euclidienne entre le premier point et le dernier point, on cherchera donc à minimiser cette valeur 
        xyz = np.array(traj.getTraj())
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        return np.sqrt((x[0]-x[-1])**2 + (y[0]-y[-1])**2 + (z[0]-z[-1])**2)


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


def voisins(table:RotTable, p: float):        
        """
        Entrée : table de rotation, p : floattant paramètre
        Sortie : Liste de 40 tables de rotation, où l'on a modifié chaque coefficient (et son complémentaire s'il existe). On fait une table ou on ajoute p*ecart_type et une ou on fait - p*ecart_type

        """
        new_table_list = []
        table_dic = table.getTable()            # On prends table_dic pour pouvoir acceder aux écarts-types
        
        # lcouple contient la liste des couple possibles suivi de leur complémentaire : "AA" suivi de "TT". /!\ On inclut pas les couples qui sont leur propres complémentaires (AT, GC, GC, TA)
        lcouple = ["AA", "TT", "AC", "GT", "AG", "CT", "CA", "TG", "CC", "GG", "GA", "TC"]
        # lcouple2 contient les couples qui sont leur propres complémentaires
        lcouple2 = ["AT", "GC", "GC", "TA"]
        for i in range(0,len(lcouple),2):
                # A chaque couple (et son complémentaire) on fait +- son écart type pour le twist et aussi pour le wedge
                couple = lcouple[i]
                compl = lcouple[i+1]
                
                # Modif le twist 
                
                new_plus_Twist = deepcopy(table)
                new_plus_Twist.setTwist(couple, table.getTwist(couple) + p*table_dic[couple][3])  # on récupere les données grace aux méthodes get.  table["AA"][3] écart type pour le twist de AA 
                new_plus_Twist.setTwist(compl, table.getTwist(couple) + p*table_dic[compl][3])  # Implique changement du complémentaire  
                new_table_list.append(new_plus_Twist)   

                new_moins_Twist = deepcopy(table)
                new_moins_Twist.setTwist(couple, table.getWedge(couple) - p*table_dic[couple][3])      
                new_moins_Twist.setTwist(compl, table.getWedge(compl) - p*table_dic[compl][3])  
                new_table_list.append(new_moins_Twist)
                
                # Modif le wedge
                
                new_plus_Wedge = deepcopy(table)
                new_plus_Wedge.setWedge(couple, table.getTwist(couple) + p*table_dic[couple][4]) 
                new_plus_Wedge.setWedge(compl, table.getTwist(compl) + p*table_dic[compl][4])  
                new_table_list.append(new_plus_Wedge)
                
                new_moins_Wedge = deepcopy(table)
                new_moins_Wedge.setWedge(couple, table.getWedge(couple) - p*table_dic[couple][4]) 
                new_moins_Wedge.setWedge(compl, table.getWedge(compl) - p*table_dic[compl][4])  
                new_table_list.append(new_moins_Wedge)
        
        for i in range(len(lcouple2)):                  # Cas des dinucléotides qui sont leur propres complémentaires
                couple = lcouple2[i]
                
                new_plus_Twist = deepcopy(table)
                new_plus_Twist.setTwist(couple, table.getTwist(couple) + p*table_dic[couple][3])  
                new_table_list.append(new_plus_Twist)

                new_moins_Twist = deepcopy(table)
                new_moins_Twist.setTwist(couple, table.getTwist(couple) - p*table_dic[couple][3]) 
                new_table_list.append(new_moins_Twist)
                
                # Modif le wedge
                
                new_plus_Wedge = deepcopy(table)
                new_plus_Wedge.setWedge(couple, table.getWedge(couple) + p*table_dic[couple][4]) 
                new_table_list.append(new_plus_Wedge)
                
                new_moins_Wedge = deepcopy(table)
                new_moins_Wedge.setWedge(couple, table.getWedge(couple) - p*table_dic[couple][4]) 
                new_table_list.append(new_moins_Wedge)
        return new_table_list
               
def recuit_simule(seq,trajectoire):
        tps_init = time.process_time()
        temps = 0
        table = RotTable()
        eval = evaluation(trajectoire)
        temperature = 100
        coeff = 1
        while(temps < 100 and temperature > 0.1):
                nombre_aleatoire = uniform(0,1) # On choisit un nombre uniformement dans [0,1]
                for table_n in voisins(table, coeff*nombre_aleatoire): # Pour tous les voisins de la table s, on calcule sa trajectoire ainsi que son énergie
                        trajectoire.compute(seq, table_n) 
                        eval_n = evaluation(trajectoire)
                        if(eval_n < eval or nombre_aleatoire < math.exp((eval-eval_n)/temperature)):
                                table = table_n
                                eval = eval_n
                temperature = 0.95*temperature # On change la temperature
                tpsi = time.process_time()
                temps = tpsi-tps_init # On actualise le temps
                coeff = 0.95*coeff    
        print(temperature)
        trajectoire = trajectoire.compute(seq,table)
        print(table.getTable())



