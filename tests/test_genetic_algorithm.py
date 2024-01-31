import unittest
from dna.genetic_algorithm import *
from dna.RotTable import RotTable

class TestGeneticAlgorithm(unittest.TestCase):
    def test_populate(self):
        rot_table = RotTable()
        population_size1 = 10
        gen_alg1 = GeneticAlgorithm(population_size1, rot_table)
        self.assertEqual(len(gen_alg1.population), population_size1)
        for table in gen_alg1.population:
            self.assertIsNotNone(table)
            for row in table.getTable().values():
                self.assertIsNotNone(row)
        
        population_size2 = 50
        gen_alg2 = GeneticAlgorithm(population_size2, rot_table)
        self.assertEqual(len(gen_alg2.population), population_size2)
        for table in gen_alg2.population:
            self.assertIsNotNone(table)
            for row in table.getTable().values():
                self.assertIsNotNone(row)
