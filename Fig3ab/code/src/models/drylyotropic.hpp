#ifndef MODELS_DRYLYOTROPIC_HPP_
#define MODELS_DRYLYOTROPIC_HPP_

#include "models.hpp"

class DryLyotropic : public Model
{
protected:

  //Defining fields of the model
  ScalarField QQxx, QNxx, QQyx, QNyx;
  ScalarField phi, phn, phi_tmp;
  ScalarField ux, uy, ux_phi, uy_phi; 
  ScalarField fusion, fusn, fusion_tmp,maxstress, fusionSum;

  //Derivatives, etc
  ScalarField HHxx, HHyx, MU, Px, Py, L1v, L2v, L3v;
  ScalarField FFx, FFy, outS, gradS;
  ScalarField dxQQxx, dyQQxx, dxQQyx, dyQQyx;
  ScalarField del2QQxx, del2QQyx;
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;

  /** Model parameters */
  double GammaQ, xi, friction, LL=0, K1, K2, K3, CC, zeta, alpha, omega, thrStress, epsilon=1.0, delta=0;
  double GammaPhi, Aphi, Kphi, Vphi;
  //Initial Configuration options 
  double level, conc=1.0, angle_deg, angle, noise, radius;
  std::string init_config;

  double fC=1, aF=0;
  double tmix = 1000000, ton=10000;
  /** Settings and checks */  
  unsigned npc = 1; 
  int n_preinit = 1000; bool preinit_flag = false;  
  bool backflow_on = true;
  bool conserve_phi = false; bool surftension_on = true;

  double ptot = 0;  
  double totalphi=0., countphi=0.;
  double gam=0;

  int inits=0;

  double ds=0.01;




  /** Update fields using predictor-corrector method
   *
   * Because of the way the predictor-corrector is implemented this function
   * can be called many times in a row in order to import the numerical
   * stability of the algorithm. Only the first call needs to have the parameter
   * set to true.
   * */
  virtual void UpdateNematicFields(bool);
  virtual void UpdateFluidFields(bool);
  virtual void UpdateFusionFields(bool);

  /** Compute chemical potential, stress and derivatives */
  virtual void UpdateNematicQuantities();
  virtual void UpdateFluidQuantities();
  /** UpdateFields() implementation */
  void UpdateNematicFieldsAtNode(unsigned, bool);
  void UpdateFluidFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdateNematicQuantitiesAtNode(unsigned);
  void UpdateFluidQuantitiesAtNode(unsigned);
  void UpdateFusionFieldsAtNode(unsigned, bool);
  void CalcMaxStress();
  void CalcTotalFusion();

  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields2();

public:
  DryLyotropic() = default;
  DryLyotropic(unsigned, unsigned, unsigned);
  DryLyotropic(unsigned, unsigned, unsigned, GridType);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example DryLyotropicFreeBoundary.
   * */
  virtual void ConfigureAtNode(unsigned);
  virtual void ConfigureAtNode1(unsigned, std::vector<float>, std::vector<float>);
  virtual std::vector<float> ParseInitConds(std::string);
  

  // functions from base class Model
  virtual void Initialize();
  virtual void Step();
  virtual void Configure();
  virtual void RuntimeChecks();
  virtual option_list GetOptions();

  /** Serialization of parameters (do not change) */
  template<class Archive>
  void serialize_params(Archive& ar)
  {
    ar & auto_name(level)
       & auto_name(conc)
       & auto_name(angle)
       & auto_name(noise)
       & auto_name(totalphi)
       & auto_name(GammaPhi)
       & auto_name(GammaQ)
       & auto_name(xi)
       & auto_name(zeta)
       & auto_name(friction)
       & auto_name(K1)
       & auto_name(K2)
       & auto_name(K3)
       & auto_name(Kphi)
       & auto_name(init_config)
       & auto_name(Aphi)
       & auto_name(Vphi)
       & auto_name(alpha)
       & auto_name(CC)
       & auto_name(fC)
       & auto_name(aF)
       & auto_name(thrStress)
       & auto_name(omega)
       & auto_name(epsilon);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(QQxx)
       & auto_name(QQyx)
       //& auto_name(phi)
       & auto_name(outS)
       & auto_name(fusionSum)
       & auto_name(maxstress)
       & auto_name(L3v)
       & auto_name(fusion);
  }
};

#endif//MODELS_LYOTROPIC_HPP_
