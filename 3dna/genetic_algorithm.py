from random import randint, gauss, uniform, choices
from copy import deepcopy
from .RotTable import RotTable
from .Traj3D import Traj3D

class GeneticAlgorithm:
    def __init__(self, population_size: int, og_table: RotTable) -> None:
        self.__population_size = population_size
        self._population = []
        self.og_table = og_table
        self._populate()

    @property
    def population(self):
        return self._population
    
    @population.setter
    def population(self, new_population):
        self._population = new_population

    def _populate(self):
        self.population += [uniform_randomize(self.population) for _ in range(self.__population_size)]
        self.population.append(self.og_table)
    
    def evaluate(self,seq):
        self.scores = []
        for table in self.population:
            traj = Traj3D()
            traj.compute(seq,table)
            self.scores += [traj.getLength()]


    def crossover(self):
        new_population = []
        for i in range(self.__population_size):
            parent1, parent2 = choices(self.population,k=2)
            child1, child2 = self._crossover(parent1,parent2) #TODO: Implement crossover
            new_population.append(child1,child2)
        

##############################################################################################################

def symmetrizeTable(incomplete_table: RotTable):
    symmetry = {"AA":"TT","AC":"GT","AG":"CT","CA":"TG","CC":"GG","GA":"TC"}
    table = deepcopy(incomplete_table)
    for base_pair in symmetry:
        table.setTwist(symmetry[base_pair], table.getTwist(base_pair))
        table.setWedge(symmetry[base_pair], table.getWedge(base_pair))
    
    return table


def gaussian_randomize(table: RotTable) -> RotTable:
    new_table = deepcopy(table)
    means = {
        "AA": [35.62 , 7.2 , -154],
        "AC": [34.4  , 1.1 ,  143],
        "AG": [27.7  , 8.4 ,    2],
        "AT": [31.5  , 2.6 ,    0],
        "CA": [34.5  , 3.5 ,  -64],
        "CC": [33.67 , 2.1 ,  -57],
        "CG": [29.8  , 6.7 ,    0],
        "CT": [27.7  , 8.4 ,   -2],
        "GA": [36.9  , 5.3 ,  120],
        "GC": [40    , 5   ,  180],
        "GG": [33.67 , 2.1 ,   57],
        "GT": [34.4  , 1.1 , -143],
        "TA": [36    , 0.9 ,    0],
        "TC": [36.9  , 5.3 , -120],
        "TG": [34.5  , 3.5 ,   64],
        "TT": [35.62 , 7.2 ,  154]
    }

    non_symmetric_table = table.getNonSymmetric()
    for dinucleotide, rotations in non_symmetric_table.items():
        new_table.setTwist(dinucleotide, gauss(means[dinucleotide][0], rotations[3]))
        new_table.setWedge(dinucleotide, gauss(means[dinucleotide][1], rotations[4]))
        
    return symmetrizeTable(new_table)


def uniform_randomize(table: RotTable) -> RotTable:
    #uniform: (-2sigma,+2sigma)
    new_table = deepcopy(table)
    means = {
        "AA": [35.62 , 7.2 , -154],
        "AC": [34.4  , 1.1 ,  143],
        "AG": [27.7  , 8.4 ,    2],
        "AT": [31.5  , 2.6 ,    0],
        "CA": [34.5  , 3.5 ,  -64],
        "CC": [33.67 , 2.1 ,  -57],
        "CG": [29.8  , 6.7 ,    0],
        "CT": [27.7  , 8.4 ,   -2],
        "GA": [36.9  , 5.3 ,  120],
        "GC": [40    , 5   ,  180],
        "GG": [33.67 , 2.1 ,   57],
        "GT": [34.4  , 1.1 , -143],
        "TA": [36    , 0.9 ,    0],
        "TC": [36.9  , 5.3 , -120],
        "TG": [34.5  , 3.5 ,   64],
        "TT": [35.62 , 7.2 ,  154]
    }

    non_symmetric_table = table.getNonSymmetric()
    for dinucleotide, rotations in non_symmetric_table.items():
        lb = means[dinucleotide][0] - 2*rotations[3]
        ub = means[dinucleotide][0] + 2*rotations[3]
        new_table.setTwist(dinucleotide, uniform(lb,ub))

        lb = means[dinucleotide][1] - 2*rotations[4]
        ub = means[dinucleotide][1] + 2*rotations[4]
        new_table.setWedge(dinucleotide, uniform(lb,ub))

    return symmetrizeTable(new_table)


def read_file(path):
    lineList = [line.rstrip('\n') for line in open(path)]
    seq = ''.join(lineList[1:])

    return seq
