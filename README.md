AMPLITUDE-FREQUENCY ESTIMATION


PREAMBLE:

This package contains a new computational strategy to estimate synaptic conductances by taking into account the peak amplitude and the peak frequency of the membrane potential of the target neuron. Apart from estimating the synaptic conductances, the method also discerns between the excitatory and the inhibitory contributions. 

Next, we describe all functions contained in the package.



(1) MAIN FUNCTION:

(1.1) EstimationProcedure(A,T,gE,gI,V, t)
	
- description:

  Routine to estimate the excitatory and inhibitory conductances from a given membrane potential trace.

- input:
  
  A=nxm matrix containing the v(t) amplitude corresponding to the (gE,gI) pair
  
  T=nxm matrix containing the v(t) period corresponding to the (gE,gI) pair
  
  gE=1xm vector containing discretized values of excitatory conductance gE
   
  gI=1xn vector containing discretized values of inhibitory conductance gI
  
  V=1xp vector of the membrane potential from which we want to estimate gE and gI
  
  t=1xp vector corresponding to the time vector of V
  
- output:
  
  EstimatedgE:1xp vector of the estimated gE conductance obtained
  
  EstimatedgI:1xp vector of the estimated gI conductance obtained


(2) AUXILIAR FUNCTIONS:

(2.1) TableSearch(gE,gI,T,t,A,a)
	
- description:

    Routine to find the closest pair (gE,gI) that can be found in tables A and T providing an amplitude value equal to 'a' and a period value equal to 't'. 

    If more than one pair (gE,gI) is found, the code returns the first (gE,gI) found.

- input:

    A:nxm matrix containing the v(t) amplitude corresponding to the (gE,gI) pair 
  
    T:nxm matrix containing the v(t) period corresponding to the (gE,gI) pair
  
    gE:1xm vector containing discretized values of excitatory conductance gE 
  
    gI:1xn vector containing discretized values of inhibitory conductance gI
  
    t:current interspike period for which we want to find the correspondig pair (gE,gI) in tables A and T
  
    a:current spike amplitude for which we want to find the correspondig pair (gE,gI) in tables A and T

- output:

    outgE:closest value of gE found in the tables providing an amplitude value equal to 'a' and a period value equal to 't'
  
    outgI:closest value of gI found in the tables providing an amplitude value equal to 'a' and a period value equal to 't'

    numCases:number of pairs (gE,gI) providing the same pair (a,t)

	
(2.2) InterpNewton(gE,gI,T,t,A,a,gE0,gI0,nmax,tol)
	
- description:

    Numerical Newton's method to solve system (A-a,T-t)=(0,0) where A and T depends on the unknowns (gE,gI)

- input:

    A:nxm matrix containing the v(t) amplitude corresponding to the (gE,gI) pair 

    T:nxm matrix containing the v(t) period corresponding to the (gE,gI) pair

    gE:1xm vector containing discretized values of excitatory conductance gE 

    gI:1xn vector containing discretized values of inhibitory conductance gI

    t:current interspike period for which we want to find the corresponding pair (gE,gI) in tables A and T

    a:current spike amplitude for which we want to find the corresponding pair (gE,gI) in tables A and T

    gE0:initial gE passed as a parameter to the Newton's method

    gI0:initial gI passed as a parameter to the Newton's method

    nmax:maximum number of iterations that can be performed when applying the Newton's method

    tol:Newton's tolerance

- output:

    outgE:component to the gE component when solving the system (A-a,T-t)=(0,0)

    outgI:component to the gI component when solving the system (A-a,T-t)=(0,0)

	
(2.3) rk45(RHS, t0, x0, tf, N ,gEv, gIv, param):
	
- description:

    Numerical analysis algorithm for the numerical solution of ordinary differential equations describing the dynamics of the membrane potential and the gatting variables of a neuron given a specific base model

- input:

    RHS:base model describing spiking neuronal activity 

    t0:initial time considered for the simulation

    x0:initial conditions of the membrane potential and the gating variables of the system to be solved

    tf:final time considered for the simulation

    N:length of the time vector considered for the simulation

    gEv:vector containing discretized values of excitatory conductance gE

    gIv:vector containing discretized values of inhibitory conductance gI
  
    param:biophysical parameters of the base model
  
- output:
  
    wi:system solution vector

    ti:system solution time vector



(3) OTHER FUNCTIONS:

(3.1) StellateModelOriginal(t,x,g,param)
	
- description:

    Base model used to describe the spiking neuronal activity 

- input:

    t: time variable

    x: vector of initial conditions of each variable in the system

    g: vector containing the excitatory and inhibitory synaptic conductances g=(gE,gI)

    param: biophysical parameters of the base model

- output:
  
    dx: vector with the equations of the dynamical system describing the neuronal activity


(3.2) Results:

- description:

    Function to obtain the matrices of v(t) amplitude, A(gE,gI), and period, T(gE,gI), used to estimate the conductances gE, gI. It uses the Stellate Model but it can be changed for any other model by using a different model function.
    
- output:

    A:nxm matrix containing the v(t) amplitude corresponding to the (gE,gI) pair 

    T:nxm matrix containing the v(t) period corresponding to the (gE,gI) pair

    gE:1xm vector containing discretized values of excitatory conductance gE used to compute A and T

    gI:1xm vector containing discretized values of inhibitory conductance gI used to compute A and T

  
(3.3) Example:

- description: 
 
  	Program to illustrate the results of executing the main function in different scenarious of activity. Four options are presented to simulate the excitatory and inhibitory conductance traces, gE and gI respectively, and to see the result of the estimation procedure in each case.

    The user can choose from the following options:

  	0:constant conductances
  
  	1:simple frequency conductances
  
  	2:double frequency conductances
  
  	3:in silico conductances

