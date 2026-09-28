#ifndef MODELS_LYOADVECTION_HPP_
#define MODELS_LYOADVECTION_HPP_

#include "models.hpp"
#include <vector>


class Lyoadvection : public Model
{
protected:
  /** Lattice Boltzmann distribution
   *
   * Written as ff[k][v], where k is node and v is the direction.
   * */ 
  LBField ff, fn, ff_tmp, fn_tmp;

  //Defining fields of the model
  ScalarField QQxx, QNxx, QQyx, QNyx;
  ScalarField phi, phn, phi_tmp;
  ScalarField ux, uy, ux_phi, uy_phi; 
  ScalarField n;

  //Derivatives, etc
  ScalarField HHxx, HHyx, MU;
  ScalarField FFx, FFy;
  ScalarField dxQQxx, dyQQxx, dxQQyx, dyQQyx;
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;

  /** Model parameters */
  double rho = 40.;
  double GammaP, GammaQ, xi, tauNem, tauIso, friction1, friction2, LL1, LL2, L0, KK, AA, CC, C0, zeta1,zeta2;

  //Initial Configuration options 
  double level, conc=1.0, angle_deg, angle, noise, radius;
  std::string init_config;

  
  /** Settings and checks */  
  unsigned npc = 1; 
  int n_preinit = 1000; bool preinit_flag = false;  
  double tmix = 1000000;
  bool backflow_on = true; bool isGuo = true;
  bool conserve_phi = false;

  double ftot = 0; double ptot = 0;  
  double totalphi=0., countphi=0.;

  unsigned initson=0;



  /** Update fields using predictor-corrector method
   *
   * Because of the way the predictor-corrector is implemented this function
   * can be called many times in a row in order to import the numerical
   * stability of the algorithm. Only the first call needs to have the parameter
   * set to true.
   * */
  virtual void UpdateNematicFields(bool);
  virtual void UpdateFluidFields(bool);
  /** Compute chemical potential, stress and derivatives */
  virtual void UpdateNematicQuantities();
  virtual void UpdateFluidQuantities();
  /** UpdateFields() implementation */
  void UpdateNematicFieldsAtNode(unsigned, bool);
  void UpdateFluidFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdateNematicQuantitiesAtNode(unsigned);
  void UpdateFluidQuantitiesAtNode(unsigned);

  /** Boundary Conditions for the flow */
  virtual void BoundaryConditionsLB();
  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields2();
  /** Move the LB particles */
  void Move();

public:
  Lyoadvection() = default;
  Lyoadvection(unsigned, unsigned, unsigned);
  Lyoadvection(unsigned, unsigned, unsigned, GridType);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example LyoadvectionFreeBoundary.
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
       & auto_name(rho)
       & auto_name(GammaP)
       & auto_name(GammaQ)
       & auto_name(xi)
       & auto_name(zeta1)
       & auto_name(zeta2)
       & auto_name(tauNem)
       & auto_name(tauIso)
       & auto_name(friction1)
       & auto_name(friction2)
       & auto_name(LL1)
       & auto_name(LL2)
       & auto_name(L0)
       & auto_name(KK)
       & auto_name(init_config)
       & auto_name(AA)
       & auto_name(CC)
       & auto_name(initson);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(QQxx)
       & auto_name(QQyx)
       & auto_name(phi);
  }
};

#endif//MODELS_LYOTROPIC_HPP_
