import unittest
from dna.Traj3D import *
from dna.recuit_simule import *
from dna.RotTable import * 

# Fichier de test pour le recuit_simule

class TestRecuitSimule(unittest.TestCase):
    
    # Test des méthodes 
    
    def test_getLength(self):
        traj = Traj3D()
        rot_table = RotTable()
        # Read file
        lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
        # Formatting
        seq = ''.join(lineList[1:])
        traj.compute(seq,rot_table)
        result = traj.getLength()
        self.assertAlmostEqual(result, 5404.77, delta=1)
        
    
    def test_getDerivatives(self):
        traj = Traj3D()
        rot_table = RotTable()
        # Read file
        lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
        # Formatting
        seq = ''.join(lineList[1:])
        traj.compute(seq,rot_table)
        result = traj.getDerivatives()
        self.assertAlmostEqual(result, 0.086, delta=0.01)
        
    # Test de voisins 
    
    def test_voisins(self):
        rot_table = RotTable()
        result = voisins(rot_table,1)
        self.assertEqual(len(result), 40)
        new_table = result[0]
        self.assertEqual(new_table.getTwist("AA"), 35.68)
        self.assertEqual(new_table.getWedge("AA"), 7.2)
        self.assertEqual(new_table.getTwist("TT"), 35.68)
        
    # Test de recuit_simule (permet de vérifier l'éxécution et que on a bien un distance qui diminiue) 
    
    def test_recuit_simule(self):
        traj = Traj3D()
        rot_table = RotTable()
        # Read file
        lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
        # Formatting
        seq = ''.join(lineList[1:])
        traj.compute(seq,rot_table)
        recuit_simule(seq,traj,1,20)
        result = traj.getLength()
        self.assertLessEqual(result, 100)           # Compliqué de tester quelque chose dans lequel il y a de l'aléatoire, mais j'ai mis 100 je pense que c'est une valeur raisonnable dans le sens où on part de 5000 et qu'on veut un test rapide 
         
        
    
if __name__ == '__main__':
    unittest.main()