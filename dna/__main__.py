from .RotTable import RotTable
from .Traj3D import Traj3D
from .genetic_algorithm import *
from .recuit_simule import *
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("algo_choice", help="key to the desired algorithm to be executed (1 - Genetic algorithm, 2 - Simulated annealing)")
parser.add_argument("filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()

def main():

    ####### Parameters for Genetic Algorithm #########
    population_size = 40 # Size of the population
    mutation_prob = 0.3 # Mutation probabality
    ##################################################

    ###### Parameters for Simulated Annealing ########
    init_temp = 1 # Initial temperature
    max_time = 60 # Maximum execution time
    coeff = 3    # Coefficient (see report)
    ##################################################

    rot_table = RotTable()
    traj = Traj3D()

    # Read file
    lineList = [line.rstrip('\n') for line in open(args.filename)]
    # Formatting
    seq = ''.join(lineList[1:])
    traj.compute(seq, rot_table)

    choice = args.algo_choice
    try:
        match int(choice):
            case 1:
                algo = GeneticAlgorithm(population_size, rot_table, mutation_prob, seq, traj)
                algo.run()
                table, score=algo.get_results()

                print(f"Best score: {score}")
                print(f"Best table: {table.getTable()}")
                algo.write_results("./dna/results.json")

                table=RotTable("./dna/results.json")
                traj.compute(seq,table)
                traj.draw()
            case 2:
                recuit_simule(seq, traj, init_temp, max_time, coeff)
                print(traj.getLength2())
                print(traj.getAngle())
                traj.draw()
            case _:
                print('Invalid choice! Please choose 1 for Genetic algorithm and 2 for Simulated annealing')

    except ValueError:
        print("Invalid command. You should use a command in the format 'python -m dna <number> <input_file>'")


if __name__ == "__main__" :
    main()
