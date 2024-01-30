from random import getrandbits

class GeneticAlgorithm:
    def __init__(self, population_size) -> None:
        self.__population_size = population_size
        self._population = []
        self._populate()

    @property
    def population(self):
        return self._population
    
    @population.setter
    def population(self, new_population):
        self._population = new_population

    def _populate(self):
        self.population += [getrandbits(32) for _ in range(self.__population_size)]

alg = GeneticAlgorithm(10)
print(*[f'{gene:032b}' for gene in alg.population],sep='\n')
