# PhD-Py-Modules
My Modules used for Numerics for my PhD


Hamiltonians:

	oned_h(oned_h(structure ,N , inpt=None, oupt=None,  A=(1/50), B=(1/93), Zin=(1/50), Zout=(1/50))
		This function requires inputs of structure (list) and N (number of connections/cables), it creates a Hamiltonian for a 1D structures, and can include optional arguments for input and output sites, 
		it currently can only take 2 values for the impedances, A and B, and by default these are set to the two impedances of RG58 and RG62 cables respectively. It returns two values, the adjancy matrix,
		and the scaled Hamiltonian.
	
	bdg_h(i, N, u, t, d)
		This functions requires a positonal argument 'i', length argument, N, and hopping terms, u, t, d, and then generates a Bogliobov de-Gennes Hamiltonian, it then returns this Hamiltonian with an applied scaling.
