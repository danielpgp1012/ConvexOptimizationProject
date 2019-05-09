# ConvexOptimizationProject
Convex optimization of Biological systems. 
We are trying to model the growth of E-Coli by doing a flux balance analysis (FBA). 
We model the bacteria as a set of chemical reactions where nutrients are consumed to produce biomass and other biproducts. 
What we want to solve is the flux of components "x", which contain metabolites and biomass by applying different objective functions and constraints. 

One of the constraints is that the mass balance of the metabolites must be conserved: A*x = b where b is the 0 vector, A is a mxn matrix where
m=number of metabolites or components, n= number of reactions. Each column of A is a chemical reaction and each row is a component. 

