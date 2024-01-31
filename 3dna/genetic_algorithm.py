from random import randint, gauss, uniform, choice, sample, random
from copy import deepcopy
from .RotTable import RotTable
from .Traj3D import Traj3D

class GeneticAlgorithm:
    def __init__(self, population_size: int, og_table: RotTable, mutation_prob: float) -> None:
        self.__population_size = population_size
        self._population = []
        self.og_table = og_table
        self.__mutation_prob = mutation_prob
        self.scores = []
        self._populate()

    @property
    def population(self):
        return self._population
    
    @population.setter
    def population(self, new_population):
        self._population = new_population

    def _populate(self):
        self.population.append(self.og_table)
        self.population.extend([uniform_randomize(self.og_table) for _ in range(self.__population_size - 1)])
    
    def evaluate(self,seq,traj):
        self.scores = []
        for table in self.population:
            traj.compute(seq,table)
            self.scores += [traj.getLength()]
        return self.scores

    def crossover(self):
        new_population = []
        for i in range(0, self.__population_size, 2):
            parent1 = self.population[i]
            parent2 = self.population[i+1]
            child1, child2 = simple_crossover(parent1, parent2)
            new_population.extend([child1,child2])
        
        self.population = new_population
        self.population_size = len(new_population)

    def mutation(self):
        new_population = []
        for table in self.population:
            if random() <= self.__mutation_prob:
                mutated_table = mutate(table)
                new_population.append(mutated_table)
            else:
                new_population.append(table)

        self.population = new_population
    
    def selection(self,seq,traj):
        populations = self.population
        self.evaluate(seq,traj)
        score = self.scores
        population_score = [(populations[i],score[i]) for i in range(len(self.population))]
        new_population = []
        while(len(population_score)>1):
            two_random_elements = sample(population_score, 2)
            for element in two_random_elements:
                population_score.remove(element)
            t1,s1 = two_random_elements[0]
            t2,s2 = two_random_elements[1]
            if (s1>s2):
                new_population.append(t2)
                #print(f"les élements aléatoires : {deux_elements_aleatoires}, le selectionné : {(t2,s2)}")
            else:
                new_population.append(t1)
                #print(f"les élements aléatoires : {deux_elements_aleatoires}, le selectionné : {(t1,s1)}")
        for element in population_score: #Maybe change this, if we have an odd number of population,we keep the one who didn't fight 
            a,_ = element
            new_population.append(a)
        self.population = new_population
        self.population_size = len(new_population)
    
    def run(self,seq,traj):
        for i in range(self.population_size//2-1):
            self.selection(seq,traj)
            self.crossover()
            self.mutation()

    def get_results(self, seq, traj)->(RotTable, float): #returns the best table and its score
        self.evaluate(seq,traj)
        max_score = max(self.scores)
        max_index = self.scores.index(max_score)
        return self.population[max_index], max_score

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

def simple_crossover(parent1: RotTable, parent2: RotTable):
        cross_point = randint(1,9)
        non_symmetric_elements = ["AA","AC","AG","CA","CC","GA","AT","GC","CG","TA"]
        child1 = deepcopy(parent1)
        child2 = deepcopy(parent2)

        for i in range(cross_point):
            child1.setTwist(non_symmetric_elements[i], parent2.getTwist(non_symmetric_elements[i]))
            child1.setWedge(non_symmetric_elements[i], parent2.getWedge(non_symmetric_elements[i]))
            
            child2.setTwist(non_symmetric_elements[i], parent1.getTwist(non_symmetric_elements[i]))
            child2.setWedge(non_symmetric_elements[i], parent1.getWedge(non_symmetric_elements[i]))
        
        return child1, child2

def mutate(table: RotTable) -> RotTable:
    mutated_table = deepcopy(table)
    non_symmetric_table = table.getNonSymmetric()
    dinucleotide = choice(list(non_symmetric_table.keys()))

    if randint(0,1):
        twist = mutated_table.getTwist(dinucleotide)
        mutated_table.setTwist(dinucleotide, gauss(twist, non_symmetric_table[dinucleotide][3]/10))
    else:
        wedge = mutated_table.getWedge(dinucleotide)
        mutated_table.setWedge(dinucleotide, gauss(wedge, non_symmetric_table[dinucleotide][4]/10))

    return symmetrizeTable(mutated_table)