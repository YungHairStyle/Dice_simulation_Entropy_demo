import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

class Subsystem:
    def __init__(self, energy):
        self.values = np.array([])  # Initially empty array of rolls
        self.counter = 0            # Counts how many times subsystem appears in system
        self.energy = energy  # Initialize with the sum of the partition
        self.p = np.array([])  # List to store probabilities
        self.p_i = np.array([]) # Statistical probability model
        self.S = 0              # Entropy
        self.beta = 0           # Beta 
        self.C = 0              # Heat capacity
        self.N = 0

    def calculate_p(self, E_i):
        """
        Calculate probabilities based on the subsystem's values.
        """
        counts = Counter(self.values)
        counts_list = [counts.get(face, 0) for face in E_i]  # Get counts for each face
        total_rolls = sum(counts_list)  # For normalization

        # Normalize frequencies (avoid division by zero)
        norm = [f / total_rolls for f in counts_list] if total_rolls > 0 else [0] * len(E_i)
        
        # Convert to NumPy array
        self.p = np.array(norm)

    def calculate_E(self, E_i):
        """
        Calculate the energy of the subsystem using its probabilities and E_i.
        """

        self.E = int(self.N *sum([self.p[i]*E_i[i] for i in range(len(E_i))]))

    def calculate_beta(self, E_i):
        '''
        Calculate beta for the subsystem.
        '''
        temp = []

        for i in range(len(E_i)):
            p_i = self.p[i]
            E_i_val = E_i[i]
            
            for j in range(i + 1, len(E_i)):  # Avoid comparing same indices
                p_j = self.p[j]
                E_j_val = E_i[j]
                
                if p_i > 0 and p_j > 0 and (E_j_val - E_i_val) != 0:  # Avoid log(0) & divide by zero
                    nom = np.log(p_i / p_j)
                    denom = E_j_val - E_i_val
                    temp.append(nom / denom)
        
        if len(temp) > 0:
            self.beta = np.average(temp)  # Define beta as average over all pairs
        else:
            self.beta = 0  # Default if no valid pairs found

    def calculate_S(self):
        """
        Calculate the entropy of the subsystem.
        """

        if len(self.p) > 0:
            # Filter out zero probabilities to avoid -inf
            self.p = np.clip(self.p, 1e-10, None)  # Replace zero values with a small number to avoid log(0)

            # Shannon entropy calculation
            self.S = -self.N*np.sum(self.p * np.log(self.p))

    def calculate_C(self, E_i):
        """
        Calculate the thermodynamics heat capacity of the subsystem.
        """
        self.C=np.sum(self.p*(np.array(E_i)-self.E/10)**2)*(self.beta**2)
    
    def calculate_p_i(self, E_i):
        """
        Calculate Boltzmann probabilities of the subsystem.
        """
        # Calculate the partition function Z
        Z = np.sum(np.exp(-self.beta * np.array(E_i)))  # Ensure E_i is a numpy array for efficient vectorized operations

        # Initialize self.p_i if it isn't already initialized
        self.p_i = np.exp(-self.beta * np.array(E_i)) / Z  # Vectorized approach for Boltzmann probabilities
        
    def add_values(self, group):
        """
        Add new values to the subsystem.
        """
        self.N = len(group)
        self.values = np.concatenate((self.values, group))  # Add new values
        self.counter += 1  # Increment the counter for each new group added

    def calculate(self, E_i):
        """
        Call all the calculation methods for the subsystem.
        
        Args:
        - E_i: List or array of energies corresponding to the states.
        """
        self.calculate_p(E_i)  # Calculate probabilities
        self.calculate_E(E_i)  # Calculate energy
        self.calculate_beta(E_i)  # Calculate beta (inverse temperature)
        self.calculate_S()  # Calculate entropy
        self.calculate_C(E_i) # Calculate heat cap
        self.calculate_p_i(E_i) # Calculate probabilities from beta values

