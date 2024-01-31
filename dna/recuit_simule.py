import numpy as np
from .Traj3D import * 
import time
from random import *
from copy import deepcopy

def voisins(table:RotTable, p: float):        
        """
        Entrée : table de rotation, p : floattant paramètre
        Sortie : Liste de 40 tables de rotation, où l'on a modifié chaque coefficient (et son complémentaire s'il existe). On fait une table ou on ajoute p*ecart_type et une ou on fait - p*ecart_type

        """
        new_table_list = []
        table_dic = table.getTable()            # On tranforme la table en un dictionnaire table_dic pour pouvoir acceder aux écarts-types
        
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
               
               
               
def recuit_simule(seq, trajectoire, temp_init, temps_max):
        tps_init = time.process_time()          # Permet de calculer le temps passé 
        temps = 0
        table = RotTable()
        eval = trajectoire.getLength()
        coeff = 10
        temperature = temp_init
        while(temps < temps_max and temperature > 0.1):
                nombre_aleatoire = uniform(0,1) # On choisit un nombre uniformement dans [0,1]
                for table_n in voisins(table, coeff*nombre_aleatoire): # Pour tous les voisins de la table s, on calcule sa trajectoire ainsi que son énergie
                        trajectoire.compute(seq, table_n) 
                        eval_n = trajectoire.getLength()
                        if (eval_n < eval or nombre_aleatoire < math.exp((eval-eval_n)/temperature)):            
                                # print("on s'est amélioré de", eval-eval_n)
                                table = table_n
                                eval = eval_n
                temperature = 0.99*temperature # On change la temperature
                temps = time.process_time()-tps_init # On actualise le temps
                coeff = 0.90*coeff    
        print("temperature de fin", temperature)
        print("temps d'execution", temps) 
        print(coeff)
        trajectoire = trajectoire.compute(seq,table)
        print(table.getTable())



