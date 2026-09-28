#ifndef MODELS_INCOMPDRYPOLARLYO_HPP_
#define MODELS_INCOMPDRYPOLARLYO_HPP_

#include "models.hpp"

class IncompDryPolarLyo : public Model
{
protected:
  //Defining fields of the model
  ScalarField Px, PNx, Py, PNy;
  ScalarField phi, phn, phi_tmp;
  ScalarField pres, pres_tmp;
  ScalarField ux, uy, ux_phi, uy_phi; 

  //Derivatives, etc
  ScalarField Hx, Hy, MU;
  ScalarField FFx, FFy;
  ScalarField dxPx, dyPx, dxPy, dyPy;
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;

  /** Model parameters */
  double gamma, xi, friction, Kn, Kp, CC, zeta, zetap, Vp, alpha, beta;
  double GammaPhi, Kphi, Aphi, Vphi, nu;

  //Initial Configuration options 
  double level, conc=1.0, angle_deg, angle, noise, radius, init_order;
  std::string init_config;

  
  /** Settings and checks */  
  unsigned npc = 1; 
  int n_preinit = 1000; bool preinit_flag = false;  
  bool backflow_on = true; bool flow_alignment = true; bool isGuo = true;
  bool conserve_phi = false; bool surftension_on = true;

  double ftot = 0; double ptot = 0;  
  double totalphi=0., countphi=0.;

  double pres_tmp_sum=0.;
  double pres_sum=0.;
  double accuracy_inc=0.001/LX/LY;  

  /** Switch on Polar fluctuations? (is automatically switched if x_kBT!=0)*/
  bool P_fluct = false;
  /** Strength of polar fluctuations */
  double P_kBT = 0;
  /** Type of fluctuations: 1 - angle fluctuation, 2 - P fluctuation */
  int P_fluc_mode = 1;




  /** Update fields using predictor-corrector method
   *
   * Because of the way the predictor-corrector is implemented this function
   * can be called many times in a row in order to import the numerical
   * stability of the algorithm. Only the first call needs to have the parameter
   * set to true.
   * */
  virtual void UpdatePolarFields(bool);
  virtual void UpdateFluidFields(bool);
  /** Compute chemical potential, stress and derivatives */
  virtual void UpdatePolarQuantities();
  virtual void UpdateFluidQuantities();
  /** calculate pressure and apply incompressibility */
  virtual void Incompressibility();
  /** pressure() implementation */
  void IncompressibilityAtNode(unsigned);
  void UpdatePressureAtNode(unsigned);
  /** UpdateFields() implementation */
  /** UpdateFields() implementation */
  void UpdatePolarFieldsAtNode(unsigned, bool);
  void UpdateFluidFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdatePolarQuantitiesAtNode(unsigned);
  void UpdateFluidQuantitiesAtNode(unsigned);

  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields2();

public:
  IncompDryPolarLyo() = default;
  IncompDryPolarLyo(unsigned, unsigned, unsigned);
  IncompDryPolarLyo(unsigned, unsigned, unsigned, GridType);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example LyotropicFreeBoundary.
   * */
  virtual void ConfigureAtNode(unsigned);

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
       & auto_name(gamma)
       & auto_name(xi)
       & auto_name(beta)
       & auto_name(zeta)
       & auto_name(zetap)
       & auto_name(Vphi)
       & auto_name(Vp)
       & auto_name(nu)
       & auto_name(alpha)
       & auto_name(friction)
       & auto_name(Kn)
       & auto_name(Kp)
       & auto_name(Kphi)
       & auto_name(init_config)
       & auto_name(Aphi)
       & auto_name(CC);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(Px)
       & auto_name(Py)
       & auto_name(phi)
       & auto_name(ux)
       & auto_name(uy);
  }
};

#endif//MODELS_LYOTROPIC_HPP_
