COMMENT
ENDCOMMENT
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS ALPHAI 
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
	Erev	= -80	(mV)		: reversal potential
	Prethresh = 0 			    : voltage level nec for release
	tau0 = 0.98814229249  (ms)    : open channels at end of release
	tau1 = 50.0           (ms)    : open channels at start of release
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
}

BREAKPOINT {
	SOLVE release
	i = g*(v - Erev)
}

PROCEDURE release() {
	if (pre > Prethresh) {
        h = h + f/tau0
	}

	g = g * etd + c * h
	h = h * etr

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
