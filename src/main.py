import numpy as np 
from System import System
import plot

exp = 1  # exponent for energy values, can adjust to show process for non-linear energy systems
E_i = np.array([1**(exp), 2**(exp), 3**(exp), 4**(exp), 5**(exp), 6**(exp)])  # Dice face values
N1 = 10     # Number of dice for system 1
N2 = 10     # Number of dice for system 2
N3 = N1+N2  # Number of dice for system 3
f = 100000  # Number of rolls, can adjust
system1 = System(N1, E_i, f) # Create system 1 with N1 dice and energy values E_i and f rolls
save = False     # Save images (svg format)
more_systems = False     # Compares different systems
dir = "plot"


#plot.fair_dice_check(system1, E_i, save=save, save_dir=dir, filename="fair_dice_check.svg")   # uncomment to plot fair dice check
#plot.plot_energy_frequency(system1, save=save, save_dir=dir, filename="energy_frequency.svg") # uncomment to plot energy frequency
#plot.plot_quantities(system1, save=save, save_dir=dir, filename=["S.svg","C.svg","beta.svg"]) # uncomment to plot entropy, heat capacity, and beta
#plot.plot_subsystem_details(system1, [25,30,40,50], save=save, save_dir=dir, filename="subsystem_details.svg") # uncomment to plot subsystem E_i vs p_i


if more_systems:
    system2 = System(N2, E_i, f) # Create system 2 with N2 dice and energy values E_i and f rolls
    system3 = System(N3,E_i, f) # Create system 3 with N3 dice and energy values E_i and f rolls
    plot.plot_all(system1, system2, system3, save=save, save_dir=dir, filename="Entropy_comparison.svg") # uncomment to plot entropy comparison of all systems
    for i in [50,60,70,80,90]:     # energy values for beta approach plots
        # uncomment to plot beta approach for different energy values
        plot.plot_comparison(system1, system2, system3, [i] ,save=save, save_dir=dir, filename=["Entropy_comparison_subsystems.svg","beta_comparison_subsystems.svg","beta_approach_eq.svg"])
    ## Uncomment to plot beta approach for equilibrium
    plot.plot_comparison(system1, system2, system3, [50,60,70,80,90] ,save=save, save_dir=dir, filename=["Entropy_comparison_subsystems.svg","beta_comparison_subsystems.svg","beta_approach_eq.svg"]) #for equilibrium plots

