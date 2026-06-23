import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os 

# Set default parameters for plots
tickfont = 16
titlefont = 20
labelfont = 16
legendfont = 16
size = (10,6)
scatter_color = "blue"
format = "pdf"
#format = "svg"

def fair_dice_check(system, E_i, save=False, save_dir="plots", filename="fair_dice_check."+format):
    """
    Generate a bar plot of the frequencies of E_i in all subsystem values to check the fairness of the dice.
    """

    # Aggregate all values from all subsystems
    all_values = np.concatenate([sub.values for sub in system.subsystems])

    # Count frequencies of each dice face
    counts = Counter(all_values)
    frequencies = [counts.get(face, 0) for face in E_i]
    frequencies = np.array(frequencies)
    frequencies=frequencies/sum(frequencies)

    # Plot
    plt.figure(figsize=size)  # Consistent figure size
    plt.bar(E_i, frequencies, edgecolor='black')
    plt.xlabel(f"$V_n$", fontsize=labelfont)
    plt.ylabel("Probability", fontsize=labelfont)
    #plt.title(f"Frequencies of Dice Faces", fontsize=titlefont)
    plt.xticks(E_i, fontsize=tickfont)
    plt.yticks(fontsize=tickfont)
    plt.grid(axis='y', alpha=0.7)

    plt.tight_layout()

    # Save plot as SVG if requested
    if save:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, format=format, dpi=300)
        

    # plt.show()

def plot_energy_frequency(system, save=False, save_dir="plots", filename="energy_frequency."+format, max_ticks=10):
    """
    Generate a bar plot of the frequency of energy values in all subsystems.
    """
    # Collect all energy values from all subsystems
    energies = [sub.energy for sub in system.subsystems]
    counts = [sub.counter for sub in system.subsystems]
    
    # Plot
    plt.figure(figsize=size)  # Consistent figure size
    plt.bar(energies, counts, edgecolor='black')
    plt.xlabel(f"$\\Sigma$", fontsize=labelfont)
    plt.ylabel("Frequency", fontsize=labelfont)
    #plt.title(f"Frequency of Energy Values Across Subsystems", fontsize=titlefont)
    # Reduce number of x-ticks
    if len(energies) > max_ticks:
        step = len(energies) // max_ticks  # Calculate step size for reduced ticks
        energies_to_show = energies[::step]  # Only show every `step`-th energy value
        plt.xticks(energies_to_show, fontsize=12, rotation=45)
    else:
        plt.xticks(energies, fontsize=12, rotation=45)

    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()

    # Save plot as SVG if requested
    if save:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, format=format, dpi=300)
        

    # plt.show()

def plot_subsystem_details(system, subsystems_indices, save=False, save_dir="plots", filename="subsystem_details."+format):
    """
    Generate a figure with 4 subplots corresponding to 4 selected subsystems, each showing a bar plot of self.p vs E_i.
    """
    fig, axes = plt.subplots(2, 2, figsize=size)

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    for idx, ax in zip(subsystems_indices, axes):
        sub = system.subsystems[idx]  # Get the subsystem by index
        
        # Plot p vs E_i
        ax.bar(system.E_i, sub.p, label="Data probability")
        
        # Add text with the details (S, beta, E, C)
        text = (f"S = {sub.S:.2f}\n"
                f"$\\beta$ = {sub.beta:.2f}\n"
                f"E = {sub.E}\n"
                f"C = {sub.C:.2f}")
        ax.text(0.1, 0.63, text, transform=ax.transAxes, fontsize=10,
        ha='center', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='black', boxstyle='round,pad=0.5'))

        
        
        # Set axis labels and title
        ax.set_xlabel("Energy", fontsize=labelfont)
        ax.set_ylabel("Probability", fontsize=labelfont)
        ax.set_title(f"E={sub.energy}", fontsize=titlefont)

        # Set the x-ticks and y-ticks with the proper font size
        ax.tick_params(axis='x', which='major', labelsize=tickfont)
        ax.set_ylim(0, 1)

    # Adjust layout
    plt.tight_layout()

    # Save plot as SVG if requested
    if save:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, format=format, dpi=300)

    # plt.show()

    # Create a figure with 4 subplots (2 rows and 2 columns)
    fig, axes = plt.subplots(2, 2, figsize=size)

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    filename= "boltzmann_"+filename
    
    for idx, ax in zip(subsystems_indices, axes):
        sub = system.subsystems[idx]  # Get the subsystem by index
        
        # Plot p vs E_i
        ax.bar(system.E_i, sub.p, label="Data probability")

        ax.scatter(system.E_i,sub.p_i, color ="red" , label="Boltzmann probability")
        #ax.legend(loc='upper right', fontsize=legendfont, frameon=True)
        
        
        # Add text with the details (S, beta, E, C)
        # text = (f"S = {sub.S:.2f}\n"
        #         f"$\\beta$ = {sub.beta:.2f}\n"
        #         f"E = {sub.E}\n"
        #         f"C = {sub.C:.2f}")
        text = (f"$\\beta$ = {sub.beta:.2f}")
        ax.text(0.1, 0.87, text, transform=ax.transAxes, fontsize=10,
        ha='center', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='black', boxstyle='round,pad=0.5'))

        
        
        # Set axis labels and title
        ax.set_xlabel(f"$V_n$", fontsize=labelfont)
        ax.set_ylabel("Probability", fontsize=labelfont)
        ax.set_title(f"$\\Sigma$={sub.energy}", fontsize=titlefont)

        # Set the x-ticks and y-ticks with the proper font size
        ax.tick_params(axis='x', which='major', labelsize=tickfont)
        ax.set_ylim(0, 1)

    # Adjust layout
    plt.tight_layout()

    # Save plot as SVG if requested
    if save:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, format=format, dpi=300)

    # plt.show()

def plot_quantities(system, save=False, save_dir="plots/", filename=["S."+format,"C."+format,"beta."+format]):
    E = [sub.energy for sub in system.subsystems if sub.S != 0]  # Extract energies from subsystems
    S = [sub.S for sub in system.subsystems if sub.S != 0]       # Entropy values
    C = [sub.C for sub in system.subsystems]       # Heat capacity values
    beta = [sub.beta for sub in system.subsystems if sub.beta != 0] # Beta values

    # Plot S vs E
    fig_S = plt.figure(figsize=(10, 6))
    plt.plot(E, S, color="black", linewidth = 3)
    plt.xlabel(f'$\\Sigma_R$', fontsize=labelfont)
    plt.ylabel(f'$S_R$', fontsize=labelfont)
    #plt.title('S vs E', fontsize=titlefont)
    plt.grid()
    if save:
        fig_S.savefig(f"{save_dir+'/'+filename[0]}", format=format)
    # plt.show()
    E = [sub.energy for sub in system.subsystems]
    # # Plot C vs E
    # fig_C = plt.figure(figsize=(10, 6))
    # plt.scatter(E, C, color=scatter_color)
    # plt.xlabel('Energy (E)', fontsize=labelfont)
    # plt.ylabel('Heat Capacity (C)', fontsize=labelfont)
    # #plt.title('C vs E', fontsize=titlefont)
    # plt.grid()
    # if save:
    #     fig_C.savefig(f"{save_dir+'/'+filename[1]}", format=format)
    # plt.show()
    # E = [sub.energy for sub in system.subsystems if sub.beta != 0]
    # # Plot beta vs E with best fit line and confidence interval
    # fig_beta = plt.figure(figsize=(10, 6))
    # plt.scatter(E, beta, color=scatter_color, label='Beta vs E')
    # plt.xlabel('Energy (E)', fontsize=labelfont)
    # plt.ylabel('Beta (β)', fontsize=labelfont)
    # plt.ylim(-1,1)
    # #plt.title('Beta vs E', fontsize=titlefont)
    # plt.ylim(-1,1)
    # plt.grid()
    # if save:
    #     fig_beta.savefig(f"{save_dir+'/'+filename[2]}", format=format)
    # plt.show()
    
def plot_all(system1, system2, system3, save=False, save_dir="plots", filename="Entropy_comparison."+format):
    systems = [system1,system2,system3]
    colors = ["blue","red","green"]
    plt.figure(figsize=size)
    for i in range(len(systems)):
        energies = [sub.energy for sub in systems[i].subsystems if sub.S != 0]
        entropies = [sub.S for sub in systems[i].subsystems if sub.S != 0 ]
        plt.scatter(energies,entropies, label=f"{systems[i].subsystems[0].energy} Dice", color=colors[i])
    plt.xlabel("Energy", fontsize=labelfont)
    plt.ylabel("Entropy",fontsize=labelfont)
    #plt.title("Entropy Comparison for three systems", fontsize=titlefont)
    plt.legend(loc='upper right', fontsize=legendfont, frameon=True)
    plt.grid()
    plt.tight_layout
    if save:
        save_path=os.path.join(save_dir,filename)
        plt.savefig(save_path, format=format, dpi=300)
    # plt.show()

def plot_comparison(system1, system2, system3, subsystems_indices ,save=False, save_dir="plots", filename=["Entropy_comparison_subsystems."+format,"beta_comparison_subsystems."+format,"beta_approach_eq."+format]):

    energy_vals = []
    entropy_vals = []
    color_vals = []  # Q values for coloring
    b_vals = []       # Equilibrium betas
    beta_A = []
    beta_B = []

    energies = [sub.energy for sub in system3.subsystems if sub.S!=0]
    entropies = [sub.S for sub in system3.subsystems if sub.S!=0]
    betas = [sub.beta for sub in system3.subsystems if sub.beta!=0]

    processed_pairs = set()  # Set to store unique (energy_A, energy_B) pairs

    for sub1 in system1.subsystems:
        for sub2 in system2.subsystems:
            energy_A, energy_B = sub1.energy, sub2.energy

            # Ensure uniqueness by storing sorted energy pairs
            energy_pair = tuple(sorted((energy_A, energy_B)))

            if energy_pair in processed_pairs:
                continue  # Skip already processed pairs

            processed_pairs.add(energy_pair)  # Mark this pair as processed

            energy_sum = energy_A + energy_B

            if any(np.isclose(energy_sum, subsystems_indices[i], atol=1e-6) for i in range(len(subsystems_indices))) and sub1.beta != 0 and sub2.beta!=0:
                # Assign variables
                C_A, C_B = sub1.C, sub2.C

                # Calculate the equilibrium beta
                b = ((C_A / sub1.beta) + (C_B / sub2.beta)) / (C_A/(sub1.beta**2) + C_B/(sub2.beta**2))

                # Find heat term for color map
                Q = int(np.abs(energy_A - (energy_A + energy_B) / 2))

                # Append values to respective arrays
                energy_vals.append(energy_sum)          # Sum of energies
                entropy_vals.append(sub1.S + sub2.S)    # Sum of entropies
                color_vals.append(Q)                    # Heat value for color map
                b_vals.append(b)                        # Combined heat capacities
                beta_A.append(sub1.beta)                # Initial values for beta
                beta_B.append(sub2.beta)
    
    if len(subsystems_indices)>1:
        fig1 = plt.figure(figsize=size)

        plt.plot(energies,entropies,color="black", linewidth = 3)

        # Create scatter plot with color mapping
        scatter = plt.scatter(energy_vals, entropy_vals, c=color_vals, cmap="viridis", label = f"$S_R^0+S_B^0$")

        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label(f'$\\Sigma_Q$', fontsize=labelfont)

        # Labels and title
        plt.xlabel(f"$\\Sigma_{{R+B}}$", fontsize=labelfont)
        plt.ylabel(f"$S_{{R+B}}$",fontsize=labelfont)
        #plt.title("Entropy Comparison for Systems in and Out of Equilibrium", fontsize=titlefont)
        plt.legend(fontsize=legendfont)
        plt.grid()
        plt.tight_layout
        if save:
            save_path=os.path.join(save_dir,filename[0])
            fig1.savefig(save_path, format=format, dpi=300)
        #plt.show()

        fig2= plt.figure(figsize=size)
        energies = [sub.energy for sub in system3.subsystems if sub.beta!=0]
        plt.plot(energies, betas, color="black", label=f"$\\beta^0$", linewidth = 3)

        scatter = plt.scatter(energy_vals, beta_A, c=color_vals, cmap="viridis", label=f"$\\beta^A$")
        plt.scatter(energy_vals, beta_B, c=color_vals, cmap="viridis", label=f"$\\beta^B$")

        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label(f'$\\Sigma_Q$', fontsize=labelfont)
        
        # Labels and title
        plt.xlabel(f"$\\Sigma_{{R+B}}$", fontsize=labelfont)
        plt.ylabel(f"$\\beta$",fontsize=labelfont)
        #plt.title(f"$\\beta$ Comparison for Systems in and Out of Equilibrium", fontsize=titlefont)
        plt.grid()
        # plt.legend(fontsize=legendfont) 
        plt.tight_layout
        if save:
            save_path=os.path.join(save_dir,filename[1])
            fig2.savefig(save_path, format=format, dpi=300)
        #plt.show()

    else:
        x_vals = np.arange(0,len(beta_A),1)

        # Create the plot
        plt.figure(figsize=size)
        b_vals[0]=b_vals[1]

        # Plot the first half of the data
        plt.plot(x_vals, b_vals, c="black", linewidth = 3 ,label=f"calculated $\\beta$")
        plt.axhline(y=b_vals[-1], color='r', linestyle='--', label=f"$\\beta^0$")
        scatter = plt.scatter(x_vals, beta_A, c=color_vals, cmap="viridis", label=f"$\\beta^0_A$")
        plt.scatter(x_vals, beta_B, c=color_vals, cmap="viridis", label=f"$\\beta^0_B$")
        
        # Add colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Heat transfer', fontsize=labelfont)
        
        top_text = (f"$\\beta^0_R$")
        bottom_text = (f"$\\beta^0_B$")
        plt.text(0.85, 0.91, top_text, transform=plt.gca().transAxes, fontsize=labelfont+4, ha='center', va='bottom')
        plt.text(0.85, 0.02, bottom_text, transform=plt.gca().transAxes, fontsize=labelfont+4, ha='center', va='bottom')


        # Labels and title
        plt.xlabel("Steps to equilibrium", fontsize=labelfont)
        plt.ylabel(f"$\\beta$", fontsize=labelfont)
        #plt.title(f"E = {subsystems_indices[0]} $\\beta$ Steps to Equilibrium", fontsize=titlefont)
        #plt.legend(fontsize=legendfont)
        plt.grid()
        plt.tight_layout()

        # Save if needed
        if save:
            save_path = os.path.join(save_dir, f"{subsystems_indices[0]}"+filename[2])
            plt.savefig(save_path, format=format, dpi=300)

        # Show plot
        # plt.show()


