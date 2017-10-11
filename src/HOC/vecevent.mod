:  Vector stream of events

NEURON {
	THREADSAFE
	ARTIFICIAL_CELL VecStim
	POINTER ptr
    RANGE y, dt
}

ASSIGNED {
    y
    dt    (ms)
	index
	etime (ms)
	ptr
}


INITIAL {
    y = 0
	index = 0
	element()
	if (index > 0) {
		net_send(etime - t, 1)
	}
}

NET_RECEIVE (w) {
	if (flag == 1) {
        y = 2
		net_event(t)
		element()
		if (index > 0) {
			net_send(etime - t, 1)
		}
        net_send(dt,2)
	}
    if (flag == 2) {
        y = 0
    }
}

DESTRUCTOR {
VERBATIM
	void* vv = (void*)(_p_ptr);  
        if (vv) {
		hoc_obj_unref(*vector_pobj(vv));
	}
ENDVERBATIM
}

PROCEDURE element() {
VERBATIM	
  { void* vv; int i, size; double* px;
	i = (int)index;
	if (i >= 0) {
		vv = (void*)(_p_ptr);
		if (vv) {
			size = vector_capacity(vv);
			px = vector_vec(vv);
			if (i < size) {
				etime = px[i];
				index += 1.;
			}else{
				index = -1.;
			}
		}else{
			index = -1.;
		}
	}
  }
ENDVERBATIM
}

PROCEDURE play() {
VERBATIM
	void** pv;
	void* ptmp = NULL;
	if (ifarg(1)) {
		ptmp = vector_arg(1);
		hoc_obj_ref(*vector_pobj(ptmp));
	}
	pv = (void**)(&_p_ptr);
	if (*pv) {
		hoc_obj_unref(*vector_pobj(*pv));
	}
	*pv = ptmp;
ENDVERBATIM
}
