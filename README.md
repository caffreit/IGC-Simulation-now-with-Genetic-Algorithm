# Triple-Gaussian now with Genetic Algorithm
Computes and fits IGC profile for a material described by 3 Gaussians.
Each Gaussain distribution decribes a site type; Basal Plane Defect, Basal Plane, Edge.
The positions and widths of each gaussian can be changed.
IGC Surface Energy profile is computed for the site energies given for a range of Edge sites % (equivalent to changing the flake length).
The profiles are then fitted to a logistic function. Parameters are plotted in a variety of combinations, figures are then saved.

UPDATE:
The Genetic Algorithm quickly selects the ideal parameters that give outputs that match the real data. This is useful because the parameter space is very large. The fittest individuals are selected for reproduction, with the probability of passing genetic material onto the next generation dependent on the individuals fitness. There is also a chance that a gene will mutate between generations. Overnight, it was able to find parameters better than I had found by using knowlegde of the system.
