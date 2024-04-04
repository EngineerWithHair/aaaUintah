/*
 * This project constitutes a work of the United States Government and is not
 * subject to domestic copyright protection under 17 USC Â§ 105.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef PTR_DEFS_H
#define PTR_DEFS_H

#define NDEBUG

#define PTR_NUM_MAT_PARAMS 89  // Originally it was 76, Mehmet added 12 more
#define PTR_NUM_SCALAR_ISV 68 //originally it is 14 and Qinglei adds 38 params for amorphization 9Fa+9F+18orientation+6kesai+6gama+shear+damage+ifamorphi and Mehmet added two state variables
#define PTR_NUM_TENS_ISV 1
#define PTR_FLAWHIST_OFFSET (PTR_NUM_SCALAR_ISV + 6*PTR_NUM_TENS_ISV)
#define PTR_NUM_FLAW_DIST_PARAM 12
#define PTR_NUM_FLAW_HIST_PARAM 3
#define PTR_BULKMOD_IDX 6
#define PTR_SHEARMOD_IDX 7
#define PTR_DESITY_IDX 8
#define PTR_LOCALIZED_IDX 9
#define PTR_NUM_FLAW_BIN_IDX 13
const char * const PTR_MatParamNames[PTR_NUM_MAT_PARAMS] = {
  "useDamage",
  "usePlasticity",
  "useGranularPlasticity",
  "useOldStress",
  "artificialViscosity",
  "artificialViscousHeating",
  "BulkModulus",
  "ShearModulus",
  "rho_orig",
  "FlowStress",
  "HardeningModulus",
  "InitialPlasticStrain",
  "J2RelaxationTime",
  "NumCrackFamilies",
  "MeanFlawSize",
  "FlawDensity",
  "StdFlawSize",
  "FlawDistType",
  "MinFlawSize",
  "MaxFlawSize",
  "FlawDistExponent",
  "RandomizeFlawDist",
  "RandomSeed",
  "RandomizeMethod",
  "BinBias",
  "KIc",
  "FlawFriction",
  "FlawAngle",
  "FlawGrowthExponent",
  "FlawGrowthAlpha",
  "CriticalDamage",
  "MaxDamage",
  "MicroMechPlaneStrain",
  "MaxDamageInc",               /* Remove? */
  "UseDamageTimestep",          /* Remove */
  "dt_increaseFactor",          /* Remove */
  "IncInitDamage",
  "DoFlawInteraction",
  "GPTimeConst",
  "JLoc",
  "GPGranularSlope",
  "GPCohesion",
  "GPYieldSurfType",            /* Merge with 46? */
  "GPPc",
  "GPJref",
  "GPPref",
  "GPSurfaceType",              /* 0-Two surface, 1- Single surface */
  "AbsToll",
  "RelToll",
  "MaxIter",
  "MaxLevels",
  "GFMSm0",
  "GFMSm1",
  "GFMSm2",
  "GFMSp0",
  "GFMSp1",
  "GFMSp2",
  "GFMSp3",
  "GFMSp4",
  "GFMSa1",
  "GFMSa2",
  "GFMSa3",
  "GFMSBeta",
  "GFMSPsi",
  "GFMSJ3Type",
  "ArtVisc1",
  "ArtVisc2",
  "MGC0",
  "MGGamma0",
  "MGS1",
  "MGS2",
  "MGS3",
  "MGCv",
  "MGTheta_0",
  "JMin",                       /* Remove */
  "MGNPts",                      /* Remove */
  "GPB_Ec",
  "GPB_kappa",
  "GPB_Mcs",
  "GPB_gammaB",
  "GPB_gamma",
  "GPB_theta",
  "GPB_eta",
  "GPB_Nvp",
  "GPB_por_u",   
  "GPB_por_l",        
  "GPB_porind_u",
  "GPB_porind_l",
  "GPB_por0"
  };

const char * const PTR_FlawDistParamNames[PTR_NUM_FLAW_DIST_PARAM] = {
  "NumCrackFamilies",
  "MeanFlawSize",
  "FlawDensity",
  "StdFlawSize",
  "FlawDistType",
  "MinFlawSize",
  "MaxFlawSize",
  "FlawDistExponent",
  "RandomizeFlawDist",
  "RandomSeed",
  "RandomizeMethod",
  "BinBias"
};

/* The vector of history variables consists of two parts, a fixed length portion
   representing the variable names in PTR_StateVarNames followed by a variable
   length portion that contains the flaw distribution and wing crack growth
   information. This variable length portion has the ordering shown in
   PTR_FlawDistStateBaseNames*/
const char * const PTR_StateVarNames[PTR_FLAWHIST_OFFSET] = {
  "Iel",
  "damage",
  "JGP",
  "GP_strain",
  "GP_energy",
  "plasticStrain",
  "plasticEnergy",
  "artViscPres",
  "volHeating",
  "localized",
  "epsVGP",
  "gamGP",
  "epsVGP_qs",
  "gamGP_qs",
  "sig_qs_11",
  "sig_qs_22",
  "sig_qs_33",
  "sig_qs_23",
  "sig_qs_13",
  "sig_qs_12",
  "Fa11",
  "Fa12",
  "Fa13",
  "Fa21",
  "Fa22",
  "Fa23",
  "Fa31",
  "Fa32",
  "Fa33",
  "F11",
  "F12",
  "F13",
  "F21",
  "F22",
  "F23",
  "F31",
  "F32",
  "F33",
  "n1x",
  "n1y",
  "n1z",
  "n2x",
  "n2y",
  "n2z",
  "n3x",
  "n3y",
  "n3z",
  "n4x",
  "n4y",
  "n4z",
  "n5x",
  "n5y",
  "n5z",
  "n6x",
  "n6y",
  "n6z",
  "kesai_1",
  "kesai_2",
  "kesai_3",
  "kesai_4",
  "kesai_5",
  "kesai_6",
  "gama_band_1",
  "gama_band_2",
  "gama_band_3",
  "gama_band_4",
  "gama_band_5",
  "gama_band_6",  
  "gama_shear",
  "damage_shear",
  "amorphized",
  "temp",
  "GP_porosity",
  "GP_B"
};

const char * const PTR_FlawDistStateBaseNames[PTR_NUM_FLAW_HIST_PARAM] = {
  "flawNumber",
  "starterFlawSize",
  "wingLength"
};

#endif // end ifdef PTR_DEFS_H
