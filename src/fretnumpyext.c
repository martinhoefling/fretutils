#include <Python.h>
#include <numpy/arrayobject.h>
#include "fretnumpyext.h"

/*#define USE_MKL*/

#ifdef USE_MKL
#include <stdio.h>
#include <mkl_vsl.h>
#define BUFSIZE 100
#endif

#ifdef USE_MT
#include "SFMT.h"
#endif

static PyObject *tryGetCPhoton(PyObject *self, PyObject *args)
{
	PyArrayObject *pvar;
    long startpos,rngseed;
    double pconst,*pvarp,rnd;
    npy_intp *pvarlen;
    int retval,i,j;

    if (!PyArg_ParseTuple(args,"dOll",&pconst,&pvar,&startpos,&rngseed)) return NULL;

    /*Checking if array is C-Contiguous, double and then initializing pointers*/

    if (!PyArray_ISCONTIGUOUS(pvar) || (sizeof(double) != PyArray_ITEMSIZE(pvar))) return NULL;
    pvarp=(double*)PyArray_DATA(pvar);
    pvarlen = PyArray_DIMS(pvar);

    retval=-1;

    /*Choose between RNG*/

#ifndef USE_MKL
#ifndef USE_MT
    srand(rngseed);
    for (i=startpos;i<pvarlen[0];i++){
    	rnd = (double)random() / RAND_MAX * 1.0f;

    	if (rnd < pconst) {
    		retval=0;
    		break;
    	}
    	if (rnd < pvarp[i]){
    		retval=1;
    		break;
    	}
    }
#endif
#endif

#ifdef USE_MT
    init_gen_rand(rngseed);
    for (i=startpos;i<pvarlen[0];i++){
    	rnd = genrand_res53();

    	if (rnd < pconst) {
    		retval=0;
    		break;
    	}
    	if (rnd < pvarp[i]){
    		retval=1;
    		break;
    	}
    }
#endif


#ifdef USE_MKL
    float r[BUFSIZE]; /* buffer for random numbers */
    VSLStreamStatePtr stream;
    vslNewStream( &stream, VSL_BRNG_MT19937, rngseed );

    for ( i=0; i<pvarlen[0]; i+=BUFSIZE )
    {
       vsRngUniform( VSL_METHOD_DUNIFORM_STD, stream, BUFSIZE, r, 0.0, 1.0 );
       for ( j=i; j<BUFSIZE; j++ )
       {
          rnd=r[j];
   		  if (rnd < pconst) {
			retval=0;
			break;
		  }
		  if (rnd < pvarp[i]){
			retval=1;
			break;
		  }
       }
    }

    vslDeleteStream( &stream );
#endif

    return Py_BuildValue("il",retval,i);
}

static PyMethodDef fretmethods[] = {
	{"tryGetCPhoton", tryGetCPhoton,METH_VARARGS,"Trys generation of a photon in plain C-Code"},
	{NULL,NULL,0,NULL}
	};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC initfretnumpyext(void)
{
    (void) Py_InitModule("fretnumpyext",fretmethods);
    import_array();
#ifndef USE_MT
#ifndef USE_MKL
    printf("-> FRET numpy extension initialized using rand()\n");
#endif
#endif
#ifdef USE_MT
    printf("-> FRET numpy extension initialized using SF Mersenne Twister\n");
#endif
#ifdef USE_MKL
    printf("-> FRET numpy extension initialized using Mersenne Twister from MKL\n");
#endif

}
