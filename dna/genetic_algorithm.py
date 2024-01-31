from random import randint, gauss, uniform, choice,choices, sample, random, seed
from copy import deepcopy
from .RotTable import RotTable
from .Traj3D import Traj3D

class GeneticAlgorithm:
    def __init__(self, population_size: int, og_table: RotTable, mutation_prob: float, seq: str, traj: Traj3D, seed = None) -> None:
        self.__population_size = population_size
        self._population = []
        self.scores = []
        self.og_table = og_table
        self.__mutation_prob = mutation_prob
        self.seq=seq
        self.traj=traj
        self.seed=seed
        self.scores = []
        self._populate()
        self.evaluate()

    @property
    # access to our population (list of rot_table)
    def population(self):
        return self._population
    
    @population.setter
    # allows modifying the population
    def population(self, new_population):
        self._population = new_population

    def _populate(self):
    # initialize a population with a uniform distribution around the og_table
        self.population.append(self.og_table)
        self.population.extend([uniform_randomize(self.og_table) for _ in range(self.__population_size - 1)])
    
    def evaluate(self): #passed the self.score to init
    # evaluation fonction : return a list where scores[i] is the distance between the last and first point using the rot_table i
        new_scores = []
        for table in self.population:
            self.traj.compute(self.seq,table)
            new_scores += [self.traj.getLength()]
        self.scores = new_scores

    def crossover(self):
    # shuffle the population : takes 2 individuals to create 2 children => doubling the total population
        new_population = []
        weights=[]
        max_score=max(self.scores)
        for score in self.scores:
            weights.append(max_score-score)
        for _ in range(0, self.__population_size, 2):
            parent1,parent2 = choices(self.population, k=2, weights=weights)
            while (parent1 == parent2):
                parent1,parent2 = choices(self.population, k=2, weights=weights)
            child1, child2 = double_crossover(parent1, parent2)
            new_population.extend([child1,child2])
        
        self.population += new_population
        self.__population_size = len(new_population)
        # print(f"new population size : {self.__population_size}")

    def mutation(self):
        # apply mutations with a low probability : allows us to not stuck in a local minimum 
        if self.seed is not None:
            seed(self.seed)
        new_population = []
        for table in self.population:
            if random() <= self.__mutation_prob:
                mutated_table = mutate(table)
                new_population.append(mutated_table)
            else:
                new_population.append(table)

        self.population = new_population
    
    def selection(self):
    #
        populations = self.population
        self.evaluate()
        score = self.scores
        population_score = [(populations[i],score[i]) for i in range(len(self.population))]
        new_population = []
        new_score = []
        while(len(population_score)>1):
            two_random_elements = sample(population_score, 2)
            for element in two_random_elements:
                population_score.remove(element)
            t1,s1 = two_random_elements[0]
            t2,s2 = two_random_elements[1]
            if (s1>s2):
                new_population.append(t2)
                new_score.append(s2)
                #print(f"les élements aléatoires : {deux_elements_aleatoires}, le selectionné : {(t2,s2)}")
            else:
                new_population.append(t1)
                new_score.append(s1)
                #print(f"les élements aléatoires : {deux_elements_aleatoires}, le selectionné : {(t1,s1)}")
        for element in population_score: #Maybe change this, if we have an odd number of population,we keep the one who didn't fight 
            a,s = element
            new_population.append(a)
            new_score.append(s)
        self.population = new_population
        self.__population_size = len(new_population)
        self.scores = new_score
    
    def run(self):
       i = 0
       while i<100:
        #    print(f"population size : {self.__population_size}")
        #    print(f"best score : {min(self.scores)}")
           self.selection()
           self.crossover()
           self.mutation()
           i+=1


    def get_results(self)->(RotTable, float): #returns the best table and its score
        self.evaluate()
        best_score = min(self.scores)
        min_index = self.scores.index(best_score)
        return self.population[min_index], best_score

##############################################################################################################

def symmetrizeTable(incomplete_table: RotTable):
    symmetry = {"AA":"TT","AC":"GT","AG":"CT","CA":"TG","CC":"GG","GA":"TC"}
    table = deepcopy(incomplete_table)
    for base_pair in symmetry:
        table.setTwist(symmetry[base_pair], table.getTwist(base_pair))
        table.setWedge(symmetry[base_pair], table.getWedge(base_pair))
        table.setDirection(symmetry[base_pair], -table.getDirection(base_pair))
    
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
        
        return symmetrizeTable(child1), symmetrizeTable(child2) #FIXED SYMMETRIZE

def double_crossover(parent1: RotTable, parent2: RotTable):
    cross_point1 = randint(1,9)
    cross_point2 = randint(1,9)
    non_symmetric_elements = ["AA","AC","AG","CA","CC","GA","AT","GC","CG","TA"]
    child1 = deepcopy(parent1)
    child2 = deepcopy(parent2)
    if cross_point1 == cross_point2:
        return simple_crossover(parent1,parent2)
    
    elif cross_point1 > cross_point2:
        cross_point1, cross_point2 = cross_point2, cross_point1

    for i in range(cross_point1):
        child1.setTwist(non_symmetric_elements[i], parent2.getTwist(non_symmetric_elements[i]))
        child1.setWedge(non_symmetric_elements[i], parent2.getWedge(non_symmetric_elements[i]))
        
        child2.setTwist(non_symmetric_elements[i], parent1.getTwist(non_symmetric_elements[i]))
        child2.setWedge(non_symmetric_elements[i], parent1.getWedge(non_symmetric_elements[i]))

    for j in range(cross_point2,10):
        child1.setTwist(non_symmetric_elements[j], parent2.getTwist(non_symmetric_elements[j]))
        child1.setWedge(non_symmetric_elements[j], parent2.getWedge(non_symmetric_elements[j]))
        
        child2.setTwist(non_symmetric_elements[j], parent1.getTwist(non_symmetric_elements[j]))
        child2.setWedge(non_symmetric_elements[j], parent1.getWedge(non_symmetric_elements[j]))

    return symmetrizeTable(child1), symmetrizeTable(child2)

def mutate(table: RotTable) -> RotTable:
    mutated_table = deepcopy(table)
    non_symmetric_table = table.getNonSymmetric()
    dinucleotide = choice(list(non_symmetric_table.keys()))

    if randint(0,1):
        twist = mutated_table.getTwist(dinucleotide)
        mutated_table.setTwist(dinucleotide, gauss(twist, non_symmetric_table[dinucleotide][3]))
    else:
        wedge = mutated_table.getWedge(dinucleotide)
        mutated_table.setWedge(dinucleotide, gauss(wedge, non_symmetric_table[dinucleotide][4]))

    return symmetrizeTable(mutated_table)