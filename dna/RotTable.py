import json
from os import path as os_path

here = os_path.abspath(os_path.dirname(__file__))

class RotTable:
    """Represents a rotation table"""

    # 3 first values: 3 angle values
    # 3 last values: SD values
    
    #def __init__(self, filename: str = None):
        #if filename is None:
            #filename = os_path.join(here, 'table.json')
        #self.rot_table = json.load(open(filename))


    def __init__(self, filename: str = None): #I modified the way we init a rottable to close the file after we reading it
        if filename is None:
            filename = os_path.join(here, 'table.json')
        with open(filename) as file:
            self.rot_table = json.load(file)


    ###################
    # WRITING METHODS #
    ###################
    def setTwist(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][0] = value

    def setWedge(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][1] = value

    def setDirection(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][2] = value

    def setSDTwist(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][3] = value

    def setSDWedge(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][4] = value

    def setSDDirection(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][5] = value

    ###################
    # READING METHODS #
    ###################
    def getTwist(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][0]

    def getWedge(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][1]

    def getDirection(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][2]
    
    def getSDTwist(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][3]

    def getSDWedge(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][4]

    def getSDDirection(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][5]
    
    def getTable(self) -> dict:
        return self.rot_table
    
    def getNonSymmetric(self) -> dict:
        NonSymmetric_elements = ["AA","AC","AG","CA","CC","GA","AT","GC","CG","TA"]
        table = self.getTable()
        return {elem:table[elem] for elem in NonSymmetric_elements}

    def __eq__(self, object2) -> bool:
        return self.getTable() == object2.getTable()
    
    def __str__(self):
        return str(self.getTable())

    ###################
    # TO JSON METHODS #
    ###################

    def toJSON(self, filename: str):
        with open(filename, 'w') as outfile:
            json.dump(self.getTable(), outfile, indent=4)