#ifndef MODELS_LYOTROPIC_HPP_
#define MODELS_LYOTROPIC_HPP_

#include "models.hpp"

class Lyotropic : public Model
{
protected:
  /** Lattice Boltzmann distribution
   *
   * Written as ff[k][v], where k is node and v is the direction.
   * */
  LBField ff, fn, ff_tmp, fn_tmp;
  /** Q-Tensor */
  ScalarField QQxx, QNxx, QQyx, QNyx;
  /** Binary order */
  ScalarField phi, phn, phi_tmp;

  /** Velocity */
  ScalarField ux, uy, ux_phi, uy_phi;
  /** Density */
  ScalarField n, fric2save, quad2save;
  /** Molecular field */
  ScalarField HHxx, HHyx, MU;
  /** Derivatives */
  ScalarField dxQQxx, dyQQxx, dxQQyx, dyQQyx;
  /** Stress tensor */
  ScalarField sigmaXX, sigmaYY, sigmaYX, sigmaXY;

  /** Nematic region level, concentration, director inclination, inital noise,
   * intial radius of the circle */
  double level, conc, angle_deg, angle, noise, radius;
  /** Total binary phase value*/
  double totalphi=0., countphi=0.;
  /** Fluid density */
  double rho = 40.;
  /** Fluid parameters */
  double GammaP, GammaQ, xi, tauNem, tauIso, friction, LL, KK, AA, CC, zeta;
  double L0=0;
  /** Intial configuration */
  std::string init_config;
  /** Number of correction steps in the predictor/corrector scheme */
  unsigned npc = 1;
  unsigned typpe = 1;
  double frateF = 0;
  double frateQ = 0;
  /** Sum of f (for checking purposes) */
  double ftot = 0;
  /** Total phi (for checking purposes) */
  double ptot = 0;
  /** Flag indicating if we need to conserve the phi field */
  bool conserve_phi = false;
  /** Anisotropic friction */
  double epsilon=0;
  /** Quadrupole force */
  double zeta2 = 0;
  double quad;
  int inits=0;
  double tfric=0, tfric2=0;
  double sigma=0, zeta2a=0;

  double amp=0, freq=0, base=0;

  /** Update fields using predictor-corrector method
   *
   * Because of the way the predictor-corrector is implemented this function
   * can be called many times in a row in order to import the numerical
   * stability of the algorithm. Only the first call needs to have the parameter
   * set to true.
   * */
  virtual void UpdateFields(bool);
  /** Compute chemical potential, stress and derivatives */
  virtual void UpdateQuantities();
  /** UpdateFields() implementation */
  void UpdateFieldsAtNode(unsigned, bool);
  /** UpdateQuantities() implementation */
  void UpdateQuantitiesAtNode(unsigned);

  /** Boundary Conditions for the flow */
  virtual void BoundaryConditionsLB();
  /** Boundary Conditions for the fields */
  virtual void BoundaryConditionsFields();
  /** Boundary Conditions for the secondary fields */
  virtual void BoundaryConditionsFields2();
  /** Move the LB particles */
  void Move();

public:
  Lyotropic() = default;
  Lyotropic(unsigned, unsigned, unsigned);
  Lyotropic(unsigned, unsigned, unsigned, GridType);

  /** Configure a single node
   *
   * This allows to change the way the arrays are configured in derived
   * classes, see for example LyotropicFreeBoundary.
   * */
  virtual void ConfigureAtNode1(unsigned, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>);
  virtual void ConfigureAtNode(unsigned);
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
       & auto_name(zeta)
       & auto_name(tauNem)
       & auto_name(tauIso)
       & auto_name(friction)
       & auto_name(LL)
       & auto_name(KK)
       & auto_name(init_config)
       & auto_name(AA)
       & auto_name(CC)
       & auto_name(epsilon)
       & auto_name(quad)
       & auto_name(zeta2)
       & auto_name(sigma)
       & auto_name(frateF)
       & auto_name(frateQ);
  }

  /** Serialization of the current frame (time snapshot) */
  template<class Archive>
  void serialize_frame(Archive& ar)
  {
    ar & auto_name(ff)
       & auto_name(QQxx)
       & auto_name(QQyx)
       //& auto_name(ux)
       //& auto_name(uy)
       & auto_name(phi)
       //& auto_name(quad2save)
       & auto_name(fric2save);
  }
};

#endif//MODELS_LYOTROPIC_HPP_
