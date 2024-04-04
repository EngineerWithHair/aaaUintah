#include <aba_for_c.h>  //to use FOR_NAME, which is equaivalent to umat_
#include "FastMatrix.cc"
#include "GFMS_full.cc"
#include "Matrix3x3.cc"
#include "PortableMieGruneisenEOSTemperature.cc"
#include "PortableTongeRamesh.cc"
#include "UMATTR.cc"
#include "Vector3.cc"
//#define PTR_NUM_MAT_PARAMS 76
//#include <iostream> 
//#include <vector> 


extern "C" void FOR_NAME(vumat) (int  *nblock, int *ndir, int *nshr, 
                             int *nstatev, int *nfieldv, int *nprops, 
                             int *lanneal, double *stepTime, double *totalTime, 
                             double *dt, double *cmname, void *coordMp, double *charLength,
                             double props[], void *density, 
                             double *strainInc, void *relSpinInc,
                             double *tempOld, void *stretchOld, double *defgradOld, 
                             void *fieldOld, double *stressOld, double *stateOld, 
                             double *enerInternOld, double *enerInelasOld,
                             double *tempNew, void *stretchNew, double *defgradNew, void *fieldNew,
                             //variables can be modified
                             double *stressNew, double *stateNew, double *enerInternNew, double *enerInelasNew)
{
	   
    PTR_vumat_stressUpdate(nblock, ndir, nshr, 
                            nstatev, nfieldv, nprops, 
                            lanneal, stepTime, totalTime, 
                            dt, cmname, coordMp, charLength,
                            props, density, 
                            strainInc, relSpinInc,
                            tempOld, stretchOld, defgradOld, 
                            fieldOld, stressOld, stateOld, 
                            enerInternOld, enerInelasOld,
                            tempNew, stretchNew, defgradNew, fieldNew,
                            //variables can be modified
                            stressNew, stateNew, enerInternNew, enerInelasNew);

}
  
