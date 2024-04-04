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


extern "C" void FOR_NAME(umat) (double STRESS[6], double STATEV[], double DDSDDE[6][6],
                             double *SSE, double *SPD, double *SCD, double *RPL,
                             double DDSDDT[6], double DRPLDE[6], double *DRPLDT,
                             const double STRAN[6], const double DSTRAN[6], const double TIME[2],
                             const double *DTIME, double *TEMP, double *DTEMP,
                             const double *PREDEF, const double *DPRED, const double *CMNAME,
                             const int *NDI, const int *NSHR, const int *NTENS, const int *NSTATV,
                             const double PROPS[PTR_NUM_MAT_PARAMS], const int *NPROPS, const double COORDS[3],
                             const double DROT[3][3], double *PNEWDT, const double *CELENT,
                             const double DFGRD0[3][3], const double DFGRD1[3][3],
                             const int *NOEL, const int *NPT, const int *LAYER,
                             const int *KSPT, const int *KSTEP, const int *KINC)
{
	
   //std::cout<<"PTR_NUM_MAT_PARAMS="<<PTR_NUM_MAT_PARAMS<<std::endl; 
   //std::cout<<"KSTEP="<<*KSTEP<<" KINC="<<*KINC<<std::endl; 
   //std::cout<<"t1="<<TIME[0]<<" t2="<<TIME[1]<<" dt="<<*DTIME<<std::endl; 
   //if (*KINC==2)   
   //std::cin.get();
   if (*KINC==1||*KINC==0){
      *SCD=PROPS[73];//The tempture is initialized as reference temp
      PTR_umat_getInitialValuesNoProps(*NSTATV, STATEV);
   }      
   std::ofstream amorphizationIni;
   
   double ptemp = *SCD;
   PTR_umat_stressUpdate(STRESS, STATEV, DDSDDE,
                             SSE, SPD, SCD, RPL,
                             DDSDDT, DRPLDE, DRPLDT,
                             STRAN, DSTRAN, TIME,
                             DTIME, &ptemp, DTEMP,
                             PREDEF, DPRED, CMNAME,
                             NDI, NSHR, NTENS, NSTATV,
                             PROPS, NPROPS, COORDS,
                             DROT, PNEWDT, CELENT,
                             DFGRD0, DFGRD1,
                             NOEL, NPT,LAYER,
                             KSPT, KSTEP, KINC, amorphizationIni); 
    *SCD = ptemp;

    //std::cout<<"Iel="<<STATEV[0]-1<<std::endl; 
    /*double E=460e9;
    double v=0.17;
    double mu = E/2/(1+v);
    double lambda = v*E/(1+v)/(1-2*v);
    for (int i=0; i<6; ++i){
      for (int j=0; j<6; ++j){
        DDSDDE[i][j] = 0.0;
      }
    }

    DDSDDE[0][0] = lambda + 2.0*mu;
    DDSDDE[0][1] = lambda;
    DDSDDE[0][2] = lambda;
    DDSDDE[1][0] = lambda;
    DDSDDE[1][1] = lambda + 2.0*mu;
    DDSDDE[1][2] = lambda;
    DDSDDE[2][0] = lambda;
    DDSDDE[2][1] = lambda;
    DDSDDE[2][2] = lambda + 2.0*mu;
    DDSDDE[3][3] = mu;
    DDSDDE[4][4] = mu;
    DDSDDE[5][5] = mu;  

    for (int i=0; i<6; ++i){
      for (int j=0; j<6; ++j){
        STRESS[i]=STRESS[i]+DDSDDE[i][j]*DSTRAN[j];
      }
    }*/
}
  
