from Subsystem import Subsystem
import numpy as np
import random
from itertools import combinations_with_replacement

class System: 
    def __init__(self, N, E_i, f):
        self.N = N  # Number of dice
        self.E_i = E_i  # Dice face values
        self.f = f      #Number of rolls 
        self.subsystems = []  # List to store subsystems
        self.unique_sums = set()  # Set to track unique sums

        # Generate all possible combinations of {E_i}
        for comb in combinations_with_replacement(self.E_i, self.N):
            group_sum = np.sum(comb)        #find the energy of the partition

            if group_sum not in self.unique_sums:       # make subsystem if it doesnt exist
                self.unique_sums.add(group_sum)
                sub = Subsystem(energy=group_sum)  # Create a subsystem with the energy
                self.subsystems.append(sub)         # append it to subsystems

        self.populate()
        
    
    def populate(self):
        for _ in range(self.f):  # Perform 'f' rolls
            
            # Simulate rolling N dice (random values from E_i)
            roll = random.choices(self.E_i, k=self.N)

            #roll = [unfair_die_roll() for _ in range(10)] # Uncomment to use unfair die roll
            
            group_sum = np.sum(roll)  # Calculate the sum (energy) of the group

            # Find the corresponding subsystem by energy
            for sub in self.subsystems:
                if sub.energy == group_sum:
                    sub.add_values(roll)  # Add the partition (roll) to the subsystem
                    break
        
        for i in self.subsystems:
            i.calculate(self.E_i)


def unfair_die_roll():
    # Define outcomes and their (unfair) probabilities
    outcomes = [1, 2, 3, 4, 5, 6]
    probabilities = [0.05, 0.1, 0.1, 0.15, 0.1, 0.5]  # Sum must be 1

    return random.choices(outcomes, weights=probabilities, k=1)[0]