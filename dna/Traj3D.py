import math

#For drawing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from .RotTable import RotTable


class Traj3D:
    """Represents a 3D trajectory"""

    # Vertical translation (elevation) between two di-nucleotides
    __MATRIX_T = np.array([[1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 3.38/2.0],
                       [0, 0, 0, 1]])

    def __init__(self):
        self.__Traj3D = {}
        self.fig = plt.figure()
        self.ax = plt.axes(projection='3d')

    def getTraj(self) -> dict:
        return self.__Traj3D

    def compute(self, dna_seq: str, rot_table: RotTable):

        # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = np.eye(4)  # Identity matrix

        # On enregistre la position du premier nucléotide
        self.__Traj3D = [np.array([0.0, 0.0, 0.0, 1.0])]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                theta = math.radians(rot_table.getTwist(dinucleotide)/2)
                # Create rotation matrix of theta on Z axis
                matrices_Rz[dinucleotide] = \
                    np.array([[math.cos(theta), -math.sin(theta), 0, 0],
                              [math.sin(theta), math.cos(theta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])

                alpha = math.radians((rot_table.getWedge(dinucleotide)))
                beta = math.radians((rot_table.getDirection(dinucleotide)-90))
                # Rotate of -beta on Z axis
                # Rotate of -alpha on X axis
                # Rotate of beta on Z axis
                matrices_Q[dinucleotide] = \
                    np.array([[math.cos(-beta), -math.sin(-beta), 0, 0],
                              [math.sin(-beta), math.cos(-beta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]]) \
                    @ np.array([[1, 0, 0, 0],
                                 [0, math.cos(-alpha), -math.sin(-alpha), 0],
                                 [0, math.sin(-alpha), math.cos(-alpha), 0],
                                 [0, 0, 0, 1]]) \
                    @ np.array([[math.cos(beta), -math.sin(beta), 0, 0],
                              [math.sin(beta), math.cos(beta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix = \
                total_matrix \
                @ self.__MATRIX_T \
                @ matrices_Rz[dinucleotide] \
                @ matrices_Q[dinucleotide] \
                @ matrices_Rz[dinucleotide] \
                @ self.__MATRIX_T

            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
            self.__Traj3D.append(total_matrix @ self.__Traj3D[0])

    def draw(self):
        xyz = np.array(self.__Traj3D)
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        self.ax.plot(x,y,z)
        self.ax.scatter([x[0]], [y[0]], [z[0]], color='green', marker='o', s=100)

        # Add a star marker for the end of the trajectory
        self.ax.scatter([x[-1]], [y[-1]], [z[-1]], color='red', marker='o', s=100)

        plt.show()

    def write(self, filename: str):
        self.fig.savefig(filename)

    def getLength(self, dna_seq, rot_table) -> float:
    # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = np.eye(4)  # Identity matrix

        # On enregistre la position du premier nucléotide
        self.__Traj3D = [np.array([0.0, 0.0, 0.0, 1.0])]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                theta = math.radians(rot_table.getTwist(dinucleotide)/2)
                # Create rotation matrix of theta on Z axis
                matrices_Rz[dinucleotide] = \
                    np.array([[math.cos(theta), -math.sin(theta), 0, 0],
                              [math.sin(theta), math.cos(theta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])

                alpha = math.radians((rot_table.getWedge(dinucleotide)))
                beta = math.radians((rot_table.getDirection(dinucleotide)-90))
                # Rotate of -beta on Z axis
                # Rotate of -alpha on X axis
                # Rotate of beta on Z axis
                matrices_Q[dinucleotide] = \
                    np.array([[math.cos(-beta), -math.sin(-beta), 0, 0],
                              [math.sin(-beta), math.cos(-beta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]]) \
                    @ np.array([[1, 0, 0, 0],
                                 [0, math.cos(-alpha), -math.sin(-alpha), 0],
                                 [0, math.sin(-alpha), math.cos(-alpha), 0],
                                 [0, 0, 0, 1]]) \
                    @ np.array([[math.cos(beta), -math.sin(beta), 0, 0],
                              [math.sin(beta), math.cos(beta), 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix = \
                total_matrix \
                @ self.__MATRIX_T \
                @ matrices_Rz[dinucleotide] \
                @ matrices_Q[dinucleotide] \
                @ matrices_Rz[dinucleotide] \
                @ self.__MATRIX_T

            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
        return np.linalg.norm(total_matrix @ self.__Traj3D[0] - self.__Traj3D[0])

    
    def getAngle(self) -> float:
        end, end_1, start = self.__Traj3D[-1], self.__Traj3D[-2], self.__Traj3D[0]
        vecteur1 = [end_1[0]-end[0], end_1[1]-end[1], end_1[2]-end[2]]
        vecteur2 = [start[0]-end[0], start[1]-end[1], start[2]-end[2]]
        produit_scalaire = np.dot(vecteur1, vecteur2)
        return produit_scalaire/(np.linalg.norm(vecteur1)*np.linalg.norm(vecteur2))
        # On a calculé cos(theta) où theta est l'angle formé par l'avant dernier, le dernier et le premier point de la trajectoire.
    
    def getEval(self,dna_seq, rot_table) -> float:
        # On ajoute une condition : getAngle() < -1/2 singnifie que nous points sont "relativement alignés" 
        # if self.getAngle() > -0.90:
        #     return 10000
        # else:
        return self.getLength(dna_seq, rot_table)