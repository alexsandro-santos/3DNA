from random import getrandbits,gauss,uniform
from .RotTable import RotTable
from .Traj3D import Traj3D

class GeneticAlgorithm:
    def __init__(self, population_size,og_table) -> None:
        self.__population_size = population_size
        self._population = []
        self.og_table=og_table
        self._populate()

    @property
    def population(self):
        return self._population
    
    @population.setter
    def population(self, new_population):
        self._population = new_population

    def _populate(self):
        self.population += [uniform_mutation(self.population) for _ in range(self.__population_size)]
        self.population.append(self.og_table)

def gaussian_mutation(table):
    new_table = table.copy()

    for dinucleotide in table:
        for i in range(2):
            new_table[dinucleotide][i] += gauss(0, dinucleotide[i+3])
    return new_table

def uniform_mutation(table):
    #uniform: (-2sigma,+2sigma)
    new_table = table.copy()
    for dinucleotide in table:
        for i in range(2):
            lb=-2*dinucleotide[i+3]
            ub=2*dinucleotide[i+3]
            new_table[dinucleotide][i] += uniform(lb,ub)
    return new_table