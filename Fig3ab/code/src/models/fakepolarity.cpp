#include "header.hpp"
#include "models/fakepolarity.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

FakePolarity::FakePolarity(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}
FakePolarity::FakePolarity(unsigned LX, unsigned LY, unsigned BC, GridType Type)
  : Model(LX, LY, BC, Type)
{}

void FakePolarity::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  ff.SetSize(LX, LY, Type);
  fn.SetSize(LX, LY, Type);
  ff_tmp.SetSize(LX, LY, Type);
  fn_tmp.SetSize(LX, LY, Type);
  QQxx.SetSize(LX, LY, Type);
  QQyx.SetSize(LX, LY, Type);
  QNxx.SetSize(LX, LY, Type);
  QNyx.SetSize(LX, LY, Type);
  phi.SetSize(LX, LY, Type);
  phi_tmp.SetSize(LX, LY, Type);
  phn.SetSize(LX, LY, Type);
  n.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  ux_phi.SetSize(LX, LY, Type);
  uy_phi.SetSize(LX, LY, Type);
  HHxx.SetSize(LX, LY, Type);
  HHyx.SetSize(LX, LY, Type);
  MU.SetSize(LX, LY, Type);
  dxQQxx.SetSize(LX, LY, Type);
  dyQQxx.SetSize(LX, LY, Type);
  dxQQyx.SetSize(LX, LY, Type);
  dyQQyx.SetSize(LX, LY, Type);
  sigmaXX.SetSize(LX, LY, Type);
  sigmaYY.SetSize(LX, LY, Type);
  sigmaYX.SetSize(LX, LY, Type);
  sigmaXY.SetSize(LX, LY, Type);
  Px.SetSize(LX, LY, Type);
  Py.SetSize(LX, LY, Type);

  FFx.SetSize(LX, LY, Type);
  FFy.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");
}

void FakePolarity::ConfigureAtNode(unsigned k)
{
  double nematicOrder = 1.0;
  double theta;
  double xtemp,ytemp;
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);
  xtemp=x;
  ytemp=y;

 //theta   = atan2(ytemp-LY/2,xtemp-LX/2)+angle + noise*M_PI*(random_real() - .5);
  theta   = angle + noise*M_PI*(random_real() - .5);
	
  QQxx[k] = nematicOrder*(cos(2*theta));
  QQyx[k] = nematicOrder*(sin(2*theta));
  Px[k] = cos(theta);
  Py[k] = sin(theta);
  phi[k]  = conc + noise*(random_real() - .5);
  totalphi += phi[k];
  // equilibrium dist
  ux[k] = uy[k] = ux_phi[k] = uy_phi[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
  ptot += phi[k];
}
void FakePolarity::Configure()
{
  for(unsigned k=0; k<DomainSize; ++k)
    ConfigureAtNode(k);

  //and do preinitialization
  cout << "Preinitialization started. ... ";
  preinit_flag = true;
  for (int i = 0; i< n_preinit; i++){
    Step();
  }
  preinit_flag = false;
  cout << "Preinitialization done. ... ";

}

void FakePolarity::UpdateNematicQuantitiesAtNode(unsigned k)
{
  const auto& d = get_neighbours(k);  

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k]; 
   
  const double del2Qxx  = laplacian(QQxx,  d, sD);
  const double dxQxx    = derivX   (QQxx,  d, sB);
  const double dyQxx    = derivY   (QQxx,  d, sB);
  const double del2Qyx  = laplacian(QQyx,  d, sD);
  const double dxQyx    = derivX   (QQyx,  d, sB);
  const double dyQyx    = derivY   (QQyx,  d, sB);

  const double p        = phi[k];
  const double dxPhi    = derivX   (phi, d, sB);
  const double dyPhi    = derivY   (phi, d, sB);
  const double del2p    = laplacian(phi, d, sD);

  const double term = p - Qxx*Qxx - Qyx*Qyx;
  const double Hxx = CC*term*Qxx + LL*del2Qxx;
  const double Hyx = CC*term*Qyx + LL*del2Qyx;

  const double mu = AA*p*(1-p)*(1-2*p) + CC*term - KK*del2p;

  // computation of sigma
  const double sigmaB = .5*AA*p*p*(1-p)*(1-p) - mu*p + (backflow_on? .5*CC*term*term : 0);
  const double sigmaF = (backflow_on? 2*xi*( (Qxx*Qxx-1.)*Hxx + Qxx*Qyx*Hyx ) : 0)
                        - zeta*Qxx*p + .5*KK*(dyPhi*dyPhi-dxPhi*dxPhi) //(preinit_flag? 0 : zeta*Qxx*p )
                        + (backflow_on? LL*(dyQxx*dyQxx+dyQyx*dyQyx-dxQxx*dxQxx-dxQyx*dxQyx) : 0);
  const double sigmaS = (backflow_on? 2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx) : 0) 
                        - zeta*Qyx*p - KK*dxPhi*dyPhi //(preinit_flag? 0 : zeta*Qyx )
                        - (backflow_on? 2*LL*(dxQxx*dyQxx+dxQyx*dyQyx) : 0);
  const double sigmaA = backflow_on? 2*(Qxx*Hyx - Qyx*Hxx) : 0;

  // transfer to arrays  
  HHxx[k]    =  Hxx;
  HHyx[k]    =  Hyx;
  dxQQxx[k]  =  dxQxx;
  dxQQyx[k]  =  dxQyx;
  dyQQxx[k]  =  dyQxx;
  dyQQyx[k]  =  dyQyx;
  sigmaXX[k] =  sigmaF + sigmaB;
  sigmaYY[k] = -sigmaF + sigmaB;
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;
  MU[k]      =  mu;

}
void FakePolarity::UpdateNematicQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNematicQuantitiesAtNode(k);
}

void FakePolarity::UpdateFluidQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);
  const double Fx = dxSxx + dySxy;
  const double Fy = dxSyx + dySyy;
  
  FFx[k] = Fx;
  FFy[k] = Fy;

  const double p   = phi[k];
  const double tau = tauNem*p + tauIso*(1-p);

  const double AA_LBF = isGuo? 0.5 : tau;
  // compute velocities
  const double nn = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  const double vx = (f[1] - f[2] + f[5] - f[6] - f[7] + f[8] + AA_LBF*FFx[k])/nn;
  const double vy = (f[3] - f[4] + f[5] - f[6] + f[7] - f[8] + AA_LBF*FFy[k])/nn;

  // transfer to arrays
  n[k]       =  nn;
  ux[k]      =  vx;
  uy[k]      =  vy;
  ux_phi[k]  =  vx*p;
  uy_phi[k]  =  vy*p;
}
void FakePolarity::UpdateFluidQuantities()
{
  double sum = 0;

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){
    sum = sum + phi[k];
    UpdateFluidQuantitiesAtNode(k);
  }
  countphi = sum;
}

void FakePolarity::UpdateNematicFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  
  const double vx = ux[k];
  const double vy = uy[k];

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  const double Hxx = HHxx[k];
  const double Hyx = HHyx[k];

  const double dxQxx = dxQQxx[k];
  const double dyQxx = dyQQxx[k];
  const double dxQyx = dxQQyx[k];
  const double dyQyx = dyQQyx[k];

  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  const double expansion = dxux + dyuy;
  const double shear     = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  //components of the Beris-Edwards equation
  const double Dxx = GammaQ*Hxx - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) +2*Qyx*shear -expansion);

  const double Dyx = GammaQ*Hyx - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*( Qyx*(expansion-traceQL) + 2*shear);
    
  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
    
  }
}
void FakePolarity::UpdateNematicFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNematicFieldsAtNode(k, first);
}

void FakePolarity::UpdateFluidFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  
  const double nn = n[k];
  const double vx = ux[k];
  const double vy = uy[k];
  const double p   = phi[k];

  // Polarisation
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double px = Px[k];
  const double py = Py[k];

  // fake polarity
  double S = sqrt(Qyx*Qyx + Qxx*Qxx);
  double nx = sqrt((1+Qxx/S)/2);
  double sgn = 1;
  if (Qyx<0){
    sgn = -1;
  }
  double ny = sgn*sqrt((1-Qxx/S)/2);

  double tp = vx*nx + vy*ny;
  double btm = sqrt(vx*vx+vy*vy)*sqrt(nx*nx+ny*ny);
  double th = acos(tp/btm);

  double dirchange = 0;
  if (th < M_PI/2 || th >= 3*M_PI/2){
    dirchange = 0; //Q = (qx,qy) is already pointing in the direction of U
  }
  else{
    dirchange = M_PI; //Q = (qx,qy) is pointing in the opposite direction of U, need to change
  }

  double Qangle = atan2(ny,nx) + dirchange;


  Px[k] = cos(Qangle);
  Py[k] = sin(Qangle);




  const double Fx = FFx[k] - friction*vx + p*alpha*cos(Qangle);
  const double Fy = FFy[k] - friction*vy + p*alpha*sin(Qangle);
  
  const double del2mu = laplacian(MU, d, sD);
  const double pFlux = derivX(ux_phi, d, sB) + derivY(uy_phi, d, sB);
  const double Dp = GammaP*del2mu - pFlux - ( conserve_phi ? (countphi-totalphi)/DomainSize : 0 );


  // LB Step for the fluid field
  const auto fe = GetEquilibriumDistribution(vx, vy, nn);
  const double tau = tauNem*p + tauIso*(1-p);
  if(first)
  {   
    for(unsigned v=0; v<lbq; ++v)
    {
      const double Si = isGuo? w[v]*(1-1./(2.*tau))*(Fx*xdir(v) + Fy*ydir(v))/cs2 : 0;

      fn[k][v] = ff[k][v] + .5*(fe[v]-ff[k][v])/tau + 0.5*Si;
      ff[k][v] = ff[k][v] +    (fe[v]-ff[k][v])/tau +     Si;
        // LB Corrected as per Kruger pp 233
    }
  }
  else
  {
    for(unsigned v=0; v<lbq; ++v){
      const double Si = isGuo? w[v]*(1-1./(2.*tau))*(Fx*xdir(v) + Fy*ydir(v))/cs2 : 0;    
      ff[k][v] = fn[k][v] + .5*(fe[v]-ff[k][v])/tau + 0.5*Si;
    }
  }

  if (first){
    phn[k]     = phi[k]  + .5*Dp;
    phi_tmp[k] = phi[k]  +    Dp;
  }
  else{
    phi_tmp[k] = phn[k]  + .5*Dp;
  }

}
void FakePolarity::UpdateFluidFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFluidFieldsAtNode(k, first);

  swap(phi.get_data(), phi_tmp.get_data());
}

void FakePolarity::Move()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<TotalSize; ++k)
  {
    for(unsigned v=0; v<lbq; ++v)
    {
      // advect particles
      ff_tmp[next(k, v)][v] = ff[k][v];
      fn_tmp[next(k, v)][v] = fn[k][v];
    }
  }
  // swap temp variables
  swap(ff.get_data(), ff_tmp.get_data());
  swap(fn.get_data(), fn_tmp.get_data());
}

void FakePolarity::BoundaryConditionsLB()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // free-slip channel
    case 1:
      ff.ApplyFreeSlipChannel();
      fn.ApplyFreeSlipChannel();
      break;
    // no-slip channel
    case 2:
      ff.ApplyNoSlipChannel();
      fn.ApplyNoSlipChannel();
      break;
    // free-slip box
    case 3:
      ff.ApplyFreeSlip();
      fn.ApplyFreeSlip();
      break;
    // no slip box
    case 4:
      ff.ApplyNoSlip();
      fn.ApplyNoSlip();
      break;
    // pbc with boundary layer
    default:
      ff.ApplyPBC();
      fn.ApplyPBC();
  }
}
void FakePolarity::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      QQxx.ApplyNeumannChannel();
      QQyx.ApplyNeumannChannel();
      phi .ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      QQxx.ApplyNeumann();
      QQyx.ApplyNeumann();
      phi .ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      QQxx.ApplyPBC();
      QQyx.ApplyPBC();
      phi .ApplyPBC();
  }
}
void FakePolarity::BoundaryConditionsFields2()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      uy     .CopyDerivativeChannel();
      ux     .CopyDerivativeChannel();
      uy_phi .ApplyDirichletChannel(0);
      ux_phi .ApplyDirichletChannel(0);

      MU     .ApplyNeumannChannel();
      sigmaXX.ApplyNeumannChannel();
      sigmaYY.ApplyNeumannChannel();
      sigmaYX.ApplyNeumannChannel();
      sigmaXY.ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      uy     .CopyDerivative();
      ux     .CopyDerivative();
      uy_phi .ApplyDirichlet(0);
      ux_phi .ApplyDirichlet(0);

      MU     .ApplyNeumann();
      sigmaXX.ApplyNeumann();
      sigmaYY.ApplyNeumann();
      sigmaYX.ApplyNeumann();
      sigmaXY.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      ux     .ApplyPBC();
      uy     .ApplyPBC();
      ux_phi .ApplyPBC();
      uy_phi .ApplyPBC();
      MU     .ApplyPBC();
      sigmaXX.ApplyPBC();
      sigmaYY.ApplyPBC();
      sigmaYX.ApplyPBC();
      sigmaXY.ApplyPBC();
  }
}

void FakePolarity::Step()
{
  // boundary conditions for primary fields
  BoundaryConditionsFields();
  // predictor step
  UpdateNematicQuantities();
  UpdateFluidQuantities();
  
  BoundaryConditionsFields2();

  // LB Step
  this->UpdateNematicFields(true); 
  this->UpdateFluidFields(true);  
  BoundaryConditionsLB();
  Move();

  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {    
    BoundaryConditionsFields();
    UpdateNematicQuantities();
    UpdateFluidQuantities();
    BoundaryConditionsFields2();
    this->UpdateNematicFields(true); 
    this->UpdateFluidFields(true);
  }
}

void FakePolarity::RuntimeChecks()
{
  // check that the sum of f is constant
  {
    double fcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);
    cout << "fcheck: " << fcheck << "/" << ftot << '\n';
    if(abs(ftot-fcheck)>1)
      throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
  }

  // check that phi is conserved
  {
    double pcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        pcheck += phi[k];
    cout << "pcheck: " << pcheck << "/" << ptot << '\n';
    if(abs(ptot-pcheck)>1)
      throw error_msg("phi is not conserved (", ptot, "/", pcheck, ")");
  }
}

option_list FakePolarity::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("GammaQ", opt::value<double>(&GammaQ),
     "Q-tensor mobility")
    ("GammaP", opt::value<double>(&GammaP),
     "binary mobility")
    ("xi", opt::value<double>(&xi),
     "tumbling/aligning parameter")
    ("tauNem", opt::value<double>(&tauNem),
     "nematic viscosity")
    ("tauIso", opt::value<double>(&tauIso),
     "isotropic viscosity")
    ("rho", opt::value<double>(&rho),
     "fluid density")
    ("friction", opt::value<double>(&friction),
     "friction from confinement")
    ("AA", opt::value<double>(&AA),
     "binary fluid bulk constant")
    ("CC", opt::value<double>(&CC),
     "coupling constant")
    ("LL", opt::value<double>(&LL),
     "elastic constant")
    ("KK", opt::value<double>(&KK),
     "binary gradient constant")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("alpha", opt::value<double>(&zeta),
     "polarity")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method")
    ("backflow_on", opt::value<bool>(&backflow_on),
     "Backflow flag")
    ("isGuo", opt::value<bool>(&isGuo),
     "LB forcing scheme")
    ("n_preinit", opt::value<int>(&n_preinit),
     "number of preinitialization steps");

  // init config options
  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("config", opt::value<string>(&init_config),
     "initial configuration")
    ("level", opt::value<double>(&level),
     "starting thickness of the nematic region")
    ("conc", opt::value<double>(&conc),
     "starting phi concentration of nematic region")
    ("radius", opt::value<double>(&radius),
     "radius of the initial circle")
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations");

  return { model_options, config_options };
}
