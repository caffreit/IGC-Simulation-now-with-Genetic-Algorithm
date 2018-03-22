# IGC Simulation now with Genetic Algorithm

#IGC Simulation
Computes and fits IGC profile for a material described by 3 Gaussians.
Each Gaussian distribution decribes a site type; Basal Plane Defect, Basal Plane, Edge.
The mean surface energy and standard deviation of each gaussian can be changed.
IGC Surface Energy profile is computed for the site energies given for a range of Edge sites % (equivalent to changing the flake length).
The profiles are then fitted to a stretched exponential. Parameters are plotted in a variety of combinations, figures are then saved.

#Genetic Fit
Using the code to simulate a surface energy profile, one is produced from a known surface energy distribution. That profile is fitted to a stretched exponential. The fitting parameters are entered into this code as the target.

The code then uses a genetic algorithm to find the original surface energy distribution. It guesses some distributions, simulates the profiles, fits them, compares those fitting parameters to the target parameters and mates the fittest parents. It loops and hopefully finds the a surface energy distribution close to the original distribution.


UPDATE:
The Genetic Algorithm quickly selects the ideal parameters that give outputs that match the real data. This is useful because the parameter space is very large. The fittest individuals are selected for reproduction, with the probability of passing genetic material onto the next generation dependent on the individuals fitness. There is also a chance that a gene will mutate between generations. Overnight, it was able to find parameters better than I had found by using knowlegde of the system.
