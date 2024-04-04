/*
 * The MIT License
 *
 * Copyright (c) 1997-2017 The University of Utah
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

#ifndef __COMPGENT_CONSTITUTIVE_MODEL_H__
#define __COMPGENT_CONSTITUTIVE_MODEL_H__


#include <CCA/Components/MPM/Materials/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Ports/DataWarehouseP.h>

#include <cmath>

namespace Uintah {
  class MPMLabel;
  class MPMFlags;

  /**************************************
CLASS
   CompGent
   
   Short description...

GENERAL INFORMATION

   CompGent.h

   Author: Kshitiz Upadhyay
   Hopkins Extreme Materials Institute
   The Johns Hopkins University
  

KEYWORDS
   Comp_Gent

DESCRIPTION
   Long description...
  
WARNING
  
  ****************************************/

  class CompGent : public ConstitutiveModel {
    // Create datatype for storing model parameters
  public:
    struct CMData {
      double C1; // Twice the initial shear modulus
      double Jm; // Limiting chain extensibility parameter
      double k0; // Initial bulk modulus (same as the "bulk" paramater of the Comp_NeoHookean model)
    };
  private:
    CMData d_initialData;
         
    // Prevent copying of this class
    // copy constructor
    //CompMooneyRivlin(const CompMooneyRivlin &cm);
    CompGent& operator=(const CompGent &cm);

  public:
    // constructor
    CompGent(ProblemSpecP& ps, MPMFlags* flag);
    //    CompMooneyRivlin(const CompMooneyRivlin* cm);
         
    // destructor 
    virtual ~CompGent();

    virtual void outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag = true);

    // clone

    CompGent* clone();
         
    // compute stable timestep for this patch
    virtual void computeStableTimestep(const Patch* patch,
                                       const MPMMaterial* matl,
                                       DataWarehouse* new_dw);
         
    // compute stress at each particle in the patch
    virtual void computeStressTensor(const PatchSubset* patches,
                                     const MPMMaterial* matl,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);

    // carry forward CM data for RigidMPM
    virtual void carryForward(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

    virtual double computeRhoMicroCM(double pressure,
                                     const double p_ref,
                                     const MPMMaterial* matl, 
                                     double temperature,
                                     double rho_guess);

    virtual void computePressEOSCM(double rho_m, double& press_eos,
                                   double p_ref,
                                   double& dp_drho, double& ss_new,
                                   const MPMMaterial* matl, 
                                   double temperature);
         
    virtual double getCompressibility();

    // initialize  each particle's constitutive model data
    virtual void initializeCMData(const Patch* patch,
                                  const MPMMaterial* matl,
                                  DataWarehouse* new_dw);

    virtual void addComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet* patches) const;

    virtual void addComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet* patches,
                                        const bool recursion) const;

    virtual void addParticleState(std::vector<const VarLabel*>& from,
                                  std::vector<const VarLabel*>& to);
  };
} // End namespace Uintah


#endif  // __COMPGENT_CONSTITUTIVE_MODEL_H__ 

