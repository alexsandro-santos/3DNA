from .RotTable import RotTable
from .Traj3D import Traj3D
from .recuit_simule import *
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="input filename of DNA sequence")
parser.parse_args()
args = parser.parse_args()

def dist(x,y): #On calcule la distance euclidienne entre le premier point et le dernier point, on cherchera donc à minimiser cette valeur 
        return np.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)


def main():
    #Initialisation de la table, et de la trajectoire
    rot_table = RotTable()
    traj = Traj3D()
    # Read file
    lineList = [line.rstrip('\n') for line in open(args.filename)]
    # Formatting
    seq = ''.join(lineList[1:])
    traj.compute(seq,rot_table)
    recuit_simule(seq,traj,1,60)
    print(traj.getLength())
    print(traj.getDerivatives())
    traj.draw()
    traj.write(args.filename+".png")


if __name__ == "__main__" :
    main()
