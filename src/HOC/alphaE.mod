COMMENT
ENDCOMMENT
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS ALPHAE
	POINTER pre
	RANGE g, h, f, deltat
	NONSPECIFIC_CURRENT i
	GLOBAL etd,etr,tau0,tau1
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(umhoms) = (micromho ms)
	(mM) = (milli/liter)
}

PARAMETER {
	Erev	= 0	(mV)		  : reversal potential
	Prethresh = 0 			  : voltage level nec for release
	tau0 = 0.0988142      (ms): open channels at end of release
	tau1 = 8.33333        (ms): open channels at start of release
}



ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	f 		(umhoms)		: jump size conductance
	g 		(umho)		: smooth conductance
	h 		(umho)		: coarse conductance
    deltat      (ms)
	pre 				: pointer to presynaptic variable
    etr
    etd
    c
}

INITIAL {
	g = 0
    h = 0
    etr = exp(-deltat/tau0)
	etd = exp(-deltat/tau1)
	c = tau0/(tau1-tau0) * (etd-etr)
    :printf("Initial dt: %.15e, f: %.15e\n", deltat,f)
    :printf("    tau0: %.15e, tau1: %.15e\n", tau0, tau1)
}

BREAKPOINT {
	SOLVE release
	i = g*(v - Erev)
}

PROCEDURE release() {
    :printf("    release at v:%.15e, t:%.3f\n", v, t)
    :printf("    with g:%.15e, h:%.15e \n", g, h)
	if (pre > Prethresh) {
        h = h + f/tau0
        :printf("h jumped:%.15e to %.15e \n", f/tau0, h)
	}

	g = g * etd + c * h
	h = h * etr
    :printf("    after release g:%.15e, h:%.15e \n", g, h)

	VERBATIM
	return 0;
	ENDVERBATIM
}

:   FUNCTION exptable(x) { 
:   	TABLE  FROM -10 TO 10 WITH 2000
:   
:   	if ((x > -10) && (x < 10)) {
:   
:   		exptable = exp(x)
:   	} else {
:   		exptable = 0.
:   	}
:   }
