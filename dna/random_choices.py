import random

def choose_2(lst:list, weights:list):
    choice1, choice2 =random.choices(lst,k=2,weights=weights)
    while choice1==choice2:
        choice1, choice2=random.choices(lst,k=2,weights=weights)
    return choice1, choice2

