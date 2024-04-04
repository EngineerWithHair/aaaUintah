/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/Materials/ConstitutiveModel/CompGUSS.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Patch.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Math/Matrix3.h>
//#include <Core/Math/TangentModulusTensor.h> //from the HGOC model (was added for stiffness)
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Math/MinMax.h>
#include <Core/Malloc/Allocator.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Uintah;

// Material Constants are G_inf (Quasi-static shear modulus), Jm (Limiting chain extensibility parameter),
// k11 (Linear rate-sensitivity control parameter), k21 (Non-linear rate-sensitivity control parameter), c21 (Rate-sensitivity index),
// and kappa (Initial bulk modulus).  
// G_inf = 2*C1 (neo-Hookean parameter).

CompGUSS::CompGUSS(ProblemSpecP& ps,MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  //______________________material properties
  ps->require("G_inf",d_initialData.G_inf);
  ps->require("Jm",d_initialData.Jm);
  ps->require("k11",d_initialData.k11);
  ps->require("k21",d_initialData.k21);
  ps->require("c21",d_initialData.c21);
  ps->require("kappa",d_initialData.kappa); // This is the "bulk" paramater of the Comp_NeoHookean model
}

/*
CompGUSS::CompGUSS(const CompGUSS* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.G_inf = cm->d_initialData.G_inf;
  d_initialData.Jm = cm->d_initialData.Jm;
  d_initialData.k11 = cm->d_initialData.k11;
  d_initialData.k21 = cm->d_initialData.k21;
  d_initialData.c21 = cm->d_initialData.c21;
  d_initialData.kappa = cm->d_initialData.kappa;
}
*/

CompGUSS::~CompGUSS()
  // _______________________DESTRUCTOR
{
}


void CompGUSS::outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type","CompGUSS");
  }

  cm_ps->appendElement("G_inf", d_initialData.G_inf);
  cm_ps->appendElement("Jm", d_initialData.Jm);
  cm_ps->appendElement("k11", d_initialData.k11);
  cm_ps->appendElement("k21", d_initialData.k21);
  cm_ps->appendElement("c21", d_initialData.c21);
  cm_ps->appendElement("kappa", d_initialData.kappa);
}


CompGUSS* CompGUSS::clone()
{
  return scinew CompGUSS(*this);
}

void CompGUSS::initializeCMData(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  computeStableTimestep(patch, matl, new_dw);
}

void CompGUSS::addParticleState(std::vector<const VarLabel*>& from,
                                     std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
}

void CompGUSS::computeStableTimestep(const Patch* patch,
                                          const MPMMaterial* matl,
                                          DataWarehouse* new_dw)
{
  // Time step depends on cell spacing (dx), particle velocity (pvelocity), and
  // material dilatational wave speed (sqrt(bulk modulus + (4*shear modulus/3))/density)
  // computed at every particle (see the CFL condition). Minimum of all dTs from every patch is used.

  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pmass, pvolume;
  constParticleVariable<Vector> pvelocity;

  new_dw->get(pmass,     lb->pMassLabel, pset);
  new_dw->get(pvolume,   lb->pVolumeLabel, pset);
  new_dw->get(pvelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12,1.e-12,1.e-12);

  // __Compute wave speed at each particle, store the maximum
  
  double G_inf = d_initialData.G_inf;
  double kappa = d_initialData.kappa;

  for(ParticleSubset::iterator iter = pset->begin();iter != pset->end();iter++){
    particleIndex idx = *iter;

    c_dil = sqrt((kappa + 4.*G_inf/3.)*pvolume[idx]/pmass[idx]);

    WaveSpeed=Vector(Max(c_dil+fabs(pvelocity[idx].x()),WaveSpeed.x()),
                     Max(c_dil+fabs(pvelocity[idx].y()),WaveSpeed.y()),
                     Max(c_dil+fabs(pvelocity[idx].z()),WaveSpeed.z()));
  }
  WaveSpeed = dx/WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}


void CompGUSS::computeStressTensor(const PatchSubset* patches,
                                        const MPMMaterial* matl,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
//COMPUTES THE STRESS ON ALL THE PARTICLES IN A GIVEN PATCH FOR A GIVEN MATERIAL
//CALLED ONCE PER TIME STEP CONTAINS A COPY OF computeStableTimestep
{
  for(int pp=0;pp<patches->size();pp++){
    const Patch* patch = patches->get(pp);

    double J_new,J,p;
    Matrix3 Identity; Identity.Identity();
    Matrix3 F_bar_new, F_bar;
    Matrix3 FRate_new, D_new;
    Matrix3 B_bar_new, C_bar_new, L_bar_new, D_bar_new, BDB_bar_new, BBDBB_bar_new;
    Matrix3 FRate_bar_new, CRate_bar_new, CRate_bar_new_sq;
    Matrix3 dev_B_bar_new, dev_BDB_bar_new, dev_BBDBB_bar_new;
    Matrix3 pressure, deviatoric_elasticstress, deviatoric_viscousstress;
    double invar1_bar,invar2_bar,invar7_bar,invar12_bar,invar12_bar_sq;
    double alpha_he1, alpha_ve1, alpha_ve2;
    double U,W_he,W_ve,se=0.;
    double c_dil=0.0;
    Vector WaveSpeed(1.e-12,1.e-12,1.e-12);

    Vector dx = patch->dCell();
    int dwi = matl->getDWIndex();

    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    constParticleVariable<Matrix3> deformationGradient_new;
    constParticleVariable<Matrix3> deformationGradient;
    constParticleVariable<double> pmass,pvolume;
    constParticleVariable<Vector> pvelocity; 
    ParticleVariable<Matrix3> pstress;
    ParticleVariable<double> pdTdt,p_q;
    constParticleVariable<Matrix3> velGrad;

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    
    old_dw->get(pmass,               lb->pMassLabel,               pset);
    old_dw->get(pvelocity,           lb->pVelocityLabel,           pset);
    old_dw->get(deformationGradient, lb->pDeformationMeasureLabel, pset);

    new_dw->get(pvolume,          lb->pVolumeLabel_preReloc,              pset);
    new_dw->get(deformationGradient_new,
                                  lb->pDeformationMeasureLabel_preReloc,  pset);
    new_dw->get(velGrad,          lb->pVelGradLabel_preReloc,             pset);
    
    new_dw->allocateAndPut(pstress,          lb->pStressLabel_preReloc,   pset);
    new_dw->allocateAndPut(pdTdt,           lb->pdTdtLabel,      pset);
    //new_dw->allocateAndPut(pdTdt,           lb->pdTdtLabel_preReloc,      pset);
    new_dw->allocateAndPut(p_q,             lb->p_qLabel_preReloc,        pset);

    //_____________________________________________material parameters
    double G_inf = d_initialData.G_inf;
    double Jm = d_initialData.Jm;
    double k11 = d_initialData.k11;
    double k21 = d_initialData.k21;
    double c21 = d_initialData.c21;
    double kappa = d_initialData.kappa;

    double rho_orig = matl->getInitialDensity();

    for(ParticleSubset::iterator iter = pset->begin();
        iter != pset->end(); iter++){
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Compute the following deviatoric (_bar) tensors: deformation gradient (F), its rate,
      // left/right Cauchy-Green deformation tensor (B/C), rate of change of C, velocity gradient tensor (L),
      // and the left stretch tensor (D).
      J_new = deformationGradient_new[idx].Determinant();
      //J = deformationGradient[idx].Determinant();

      F_bar_new = deformationGradient_new[idx]* pow(J_new,-(1./3.));
      //F_bar = deformationGradient[idx]* pow(J,-(1./3.));
      B_bar_new = deformationGradient_new[idx]*deformationGradient_new[idx].Transpose()*pow(J_new,-(2./3.));
      C_bar_new = deformationGradient_new[idx].Transpose()*deformationGradient_new[idx]*pow(J_new,-(2./3.));

      FRate_new = velGrad[idx] * deformationGradient_new[idx];
      D_new = (velGrad[idx] + velGrad[idx].Transpose())*.5;
      FRate_bar_new = (FRate_new - ((1./3.)*D_new.Trace()*deformationGradient_new[idx])) * pow(J_new,-(1./3.));
      //FRate_bar_new = (F_bar_new - F_bar)/delT;
      CRate_bar_new = (FRate_bar_new.Transpose() * F_bar_new) + (F_bar_new.Transpose() * FRate_bar_new);
      L_bar_new = velGrad[idx] - Identity*((1.0/3.0)*D_new.Trace());
      D_bar_new = D_new;
      //L_bar_new = FRate_bar_new * F_bar_new.Inverse();
      //D_bar_new = 0.5*(L_bar_new + L_bar_new.Transpose());

      BDB_bar_new = B_bar_new * (D_bar_new * B_bar_new);
      BBDBB_bar_new = (B_bar_new * BDB_bar_new) + (BDB_bar_new * B_bar_new);

      // Compute the invariants
      invar1_bar = C_bar_new.Trace();
      invar2_bar = 0.5*((invar1_bar*invar1_bar) - (C_bar_new*C_bar_new).Trace()); // Keeping this to allow scaling to MR type models later.
      CRate_bar_new_sq = CRate_bar_new*CRate_bar_new;
      invar7_bar = 0.5*(CRate_bar_new_sq.Trace()); // This is twice the J_2 invariant of the USS model. 
      invar12_bar = (C_bar_new*CRate_bar_new_sq).Trace(); // This is same as the J_5 invariant of the USS model.
      invar12_bar_sq = fabs(invar12_bar * invar12_bar);

      // Compute the response coefficients
      alpha_he1 = G_inf/(1 - (invar1_bar - 3.0)/Jm); // This is 2dW/dI1tilde. There will be alpha2 as well for the MR type models.
      alpha_ve1 = 4.0*k11*sqrt(fabs(invar1_bar - 3.0));
      if (invar12_bar >= -1.e-5 && invar12_bar <= 1.e-5) {
        alpha_ve2 = 0;
      } else {
      alpha_ve2 = 2.0*k21*sqrt(fabs(invar2_bar - 3.0)) * (pow(invar12_bar_sq,(c21/2.0)) / invar12_bar);
      }

      // Compute dilatational stress (pressure)
      p = 0.5*kappa*(J_new*J_new-1)/J_new;
            if (p >= -1.e-5 && p <= 1.e-5)
              p = 0.;
            pressure = Identity*p;

      // Compute deviatoric elastic stress
      dev_B_bar_new = B_bar_new - Identity*((1.0/3.0)*B_bar_new.Trace());
      deviatoric_elasticstress = (alpha_he1*dev_B_bar_new)/J_new;

      // Compute deviatoric viscous stress
      dev_BDB_bar_new = BDB_bar_new - Identity*((1.0/3.0)*BDB_bar_new.Trace());
      dev_BBDBB_bar_new = BBDBB_bar_new - Identity*((1.0/3.0)*BBDBB_bar_new.Trace());
      deviatoric_viscousstress = ((alpha_ve1 * dev_BDB_bar_new) + (alpha_ve2 * dev_BBDBB_bar_new))*(2./J_new);

      // Compute total Cauchy stress as the sum of dilatational (pressure), deviatoric (isochoric) elastic, and deviatoric viscous stresses.
      pstress[idx] = pressure + deviatoric_elasticstress + deviatoric_viscousstress;      

      // shear = 2.*c1+c2+I4tilde*(4.*d2WdI4tilde2*lambda_tilde*lambda_tilde
      //                                -2.*dWdI4tilde*lambda_tilde); // I have no clue about source of this expression // Removed it - KU 
      
      // Compute deformed volume and local wave speed
      double rho_cur = rho_orig/J_new;
      //c_dil = sqrt((Bulk+1./3.*shear)/rho_cur); // I have no clue about source of this expression. 
      c_dil = sqrt((kappa+4./3.*G_inf)*pvolume[idx]/pmass[idx]); // Assuming fixed c_dil regardless of deformation and its rate.

      //________________________________end stress


      // Compute the strain energy for all the particles
      U = 0.5*kappa*((J_new*J_new-1)/2-log(J_new));
      W_he = -0.5*G_inf*Jm*log(1.0 - (invar1_bar-3.0)/Jm);
      W_ve = (2*k11*sqrt(invar1_bar - 3.0)*invar7_bar) + ((k21/c21)*sqrt(invar2_bar - 3.0) * pow(invar12_bar,c21));
      double e = (U + W_he + W_ve)*pvolume[idx]/J_new;
      se += e;
      
      Vector pvelocity_idx = pvelocity[idx];

      WaveSpeed=Vector(Max(c_dil+fabs(pvelocity_idx.x()),WaveSpeed.x()),
                       Max(c_dil+fabs(pvelocity_idx.y()),WaveSpeed.y()),
                       Max(c_dil+fabs(pvelocity_idx.z()),WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z())/3.0;
        double c_bulk = sqrt(kappa/rho_cur);
        Matrix3 D=(velGrad[idx] + velGrad[idx].Transpose())*0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    }  // end loop over particles

    WaveSpeed = dx/WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se),      lb->StrainEnergyLabel);
    }
  }
}


void CompGUSS::carryForward(const PatchSubset* patches,
                                 const MPMMaterial* matl,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
  //_____________________________used with RigidMPM
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models 
    // when using RigidMPM.visco_one_cell_expl_parallel.ups
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model 
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());
    
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.),     lb->StrainEnergyLabel);
    }
  }
}


void CompGUSS::addComputesAndRequires(Task* task,
                                           const MPMMaterial* matl,
                                           const PatchSet* patches) const
  //___________TELLS THE SCHEDULER WHAT DATA
  //___________NEEDS TO BE AVAILABLE AT THE TIME computeStressTensor IS CALLED
{
  // Add the computes and requires that are common to all explicit 
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.

  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);

  // Other constitutive model and input dependent computes and requires
  Ghost::GhostType  gnone = Ghost::None;
}

void CompGUSS::addComputesAndRequires(Task* ,
                                           const MPMMaterial* ,
                                           const PatchSet* ,
                                           const bool ) const
  //_________________________________________here this one's empty
{
}

// The "CM" versions use the pressure-volume relationship of the CNH model
double CompGUSS::computeRhoMicroCM(double /*pressure*/, 
                                        const double /*p_ref*/,
                                        const MPMMaterial* /*matl*/,
                                        double temperature,
                                        double rho_guess)
{
/*  double rho_orig = matl->getInitialDensity();
  double kappa = d_initialData.kappa;
  
  double p_gauge = pressure - p_ref;
  double rho_cur;

  if(d_useModifiedEOS && p_gauge < 0.0) {
    double A = p_ref;           // MODIFIED EOS
    double n = p_ref/kappa;
    rho_cur = rho_orig*pow(pressure/A,n);
  } else {                      // STANDARD EOS
    rho_cur = rho_orig*(p_gauge/kappa + sqrt((p_gauge/kappa)*(p_gauge/kappa) +1));
  }*/
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR CompGUSS"
       << endl;

  double rho_cur=0.;

  return rho_cur;

}

void CompGUSS::computePressEOSCM(const double /*rho_cur*/,double& /*pressure*/, 
                                      const double /*p_ref*/,
                                      double& /*dp_drho*/, double& /*tmp*/,
                                      const MPMMaterial* /*matl*/,
                                      double temperature)
{
/*  double kappa = d_initialData.kappa;
  double rho_orig = matl->getInitialDensity();

  if(d_useModifiedEOS && rho_cur < rho_orig){
    double A = p_ref;           // MODIFIED EOS
    double n = kappa/p_ref;
    pressure = A*pow(rho_cur/rho_orig,n);
    dp_drho  = (kappa/rho_orig)*pow(rho_cur/rho_orig,n-1);
    tmp      = dp_drho;         // speed of sound squared
  } else {                      // STANDARD EOS            
    double p_g = .5*kappa*(rho_cur/rho_orig - rho_orig/rho_cur);
    pressure   = p_ref + p_g;
    dp_drho    = .5*kappa*(rho_orig/(rho_cur*rho_cur) + 1./rho_orig);
    tmp        = kappa/rho_cur;  // speed of sound squared
  }*/
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR CompGent"
       << endl;

}

double CompGUSS::getCompressibility()
{
  //return 1.0/d_initialData.kappa;
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR CompGent"
       << endl;
  return 1.0;
}


namespace Uintah {
  
#if 0
  static MPI_Datatype makeMPI_CMData()
  {
    ASSERTEQ(sizeof(CompGUSS::StateData), sizeof(double)*0);
    MPI_Datatype mpitype;
    MPI_Type_vector(1, 0, 0, MPI_DOUBLE, &mpitype);
    MPI_Type_commit(&mpitype);
    return mpitype;
  }
  
  const TypeDescription* fun_getTypeDescription(CompGUSS::StateData*)
  {
    static TypeDescription* td = 0;
    if(!td){
      td = scinew TypeDescription(TypeDescription::Other,
                                  "CompGUSS::StateData", true, &makeMPI_CMData);
    }
    return td;
  }
#endif
} // End namespace Uintah

