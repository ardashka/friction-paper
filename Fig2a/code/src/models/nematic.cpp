#include "header.hpp"
#include "models/nematic.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"


using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;
extern int t;

Nematic::Nematic(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}

void Nematic::Initialize()
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
  n.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  HHxx.SetSize(LX, LY, Type);
  HHyx.SetSize(LX, LY, Type);
  dxQQxx.SetSize(LX, LY, Type);
  dyQQxx.SetSize(LX, LY, Type);
  dxQQyx.SetSize(LX, LY, Type);
  dyQQyx.SetSize(LX, LY, Type);
  sigmaXX.SetSize(LX, LY, Type);
  sigmaYY.SetSize(LX, LY, Type);
  sigmaYX.SetSize(LX, LY, Type);
  sigmaXY.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");

  Q_fluct = (Q_kBT!=0);
  u_fluct = (u_kBT!=0);
}

void Nematic::ConfigureAtNode1(unsigned k, vector<float> Qxxinits, vector<float> Qyxinits, vector<float> Uxinits, vector<float> Uyinits)
{
  // add noise (for meta-stable configs)
  // theta is the angle of the director
  const double theta = angle + noise*M_PI*(random_real()-.5);
  //cout << "k = " << k << ", x = " << GetXPosition(k) << ", y = " << GetYPosition(k);
  QQxx[k] = Qxxinits[k];
  QQyx[k] = Qyxinits[k];
  // equilibrium dist
  ux[k] = Uxinits[k];
  uy[k] = Uyinits[k];
  //ux[k] = uy[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

void Nematic::Configure()
{
  string fileQxx = "/groups/astro/ardash/Documents/mass-initcheats/initsnematic/Qxx.txt";
  string fileQyx = "/groups/astro/ardash/Documents/mass-initcheats/initsnematic/Qxy.txt";
  string fileUx = "/groups/astro/ardash/Documents/mass-initcheats/initsnematic/ux.txt";
  string fileUy = "/groups/astro/ardash/Documents/mass-initcheats/initsnematic/uy.txt";
  string line;
  vector<float> Qxxinit = ParseInitConds(fileQxx);
  vector<float> Qyxinit = ParseInitConds(fileQyx);
  vector<float> Uxinit = ParseInitConds(fileUx);
  vector<float> Uyinit = ParseInitConds(fileUy);

  cout << "\n Done reading files\n";
  if (inits == 1){ 
    cout << "\n Using file to set initial conditions\n";
    for(unsigned k=0; k<DomainSize; ++k)
    {
      ConfigureAtNode1(k, Qxxinit, Qyxinit, Uxinit, Uyinit);
    }
  }
  else{
    cout << "\n Using default initial conditions\n";
    for(unsigned k=0; k<DomainSize; ++k){
      ConfigureAtNode(k);
    }
  }
}

void Nematic::ConfigureAtNode(unsigned k)
{
  // add noise (for meta-stable configs)
  // theta is the angle of the director
  const double theta = angle + noise*M_PI*(random_real()-.5);
  //cout << "k = " << k << ", x = " << GetXPosition(k) << ", y = " << GetYPosition(k);
  QQxx[k] = cos(2*theta);
  QQyx[k] = sin(2*theta);
  // equilibrium dist
  ux[k] = uy[k] = 0;
  n[k]  = rho;
  ff[k] = GetEquilibriumDistribution(ux[k], uy[k], n[k]);
  // compute totals for later checks
  ftot  = accumulate(begin(ff[k]), end(ff[k]), ftot);
}

vector<float> Nematic::ParseInitConds(string filename)
{
   ifstream infile;
   infile.open(filename);
   string line;
   vector<float> parsedCsv;
   while(std::getline(infile,line))
    {
      stringstream lineStream(line);
      string cell;
      size_t offset = 0;
      while(std::getline(lineStream,cell,','))
      {
        float celldbl = stof(cell,0);
        parsedCsv.push_back(celldbl);
      }
    }
   return parsedCsv; 
}


void Nematic::UpdateQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);
  const auto& f = ff[k];
  // Q-tensor
  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  // compute velocities
  const double nn = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  const double vx = (f[1] - f[2] + f[5] - f[6] - f[7] + f[8])/nn;
  const double vy = (f[3] - f[4] + f[5] - f[6] + f[7] - f[8])/nn;

  // compute derivatives etc.
  const double del2Qxx  = laplacian(QQxx,  d, sD);
  const double dxQxx    = derivX   (QQxx,  d, sB);
  const double dyQxx    = derivY   (QQxx,  d, sB);
  const double del2Qyx  = laplacian(QQyx,  d, sD);
  const double dxQyx    = derivX   (QQyx,  d, sB);
  const double dyQyx    = derivY   (QQyx,  d, sB);

  // computation of the chemical potential and molecular field...
  // ...term that couples the binary phase to the degree of nematic order
  const double term = 1. - Qxx*Qxx - Qyx*Qyx;
  // ...molecular field
  const double Hxx = CC*term*Qxx + LL*del2Qxx;
  const double Hyx = CC*term*Qyx + LL*del2Qyx;

  // computation of sigma...
  // ... on-diagonal stress components
  const double sigmaB = no_backflow*.5*CC*term*term - zetaC;
  const double sigmaF = no_backflow*2*xi*( (Qxx*Qxx-1.)*Hxx + Qxx*Qyx*Hyx )
                        - zeta*Qxx
                        + no_backflow*LL*(dyQxx*dyQxx+dyQyx*dyQyx-dxQxx*dxQxx-dxQyx*dxQyx);
  // .. off-diagonal stress components
  const double sigmaS = no_backflow*2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx) - zeta*Qyx
                        - no_backflow*2*LL*(dxQxx*dyQxx+dxQyx*dyQyx);
  const double sigmaA = no_backflow*2*(Qxx*Hyx - Qyx*Hxx);


  // transfer to arrays
  n[k]       =  nn;
  ux[k]      =  vx;
  uy[k]      =  vy;
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
}

void Nematic::UpdateFieldsAtNode(unsigned k, bool first)
{
  // pointer to neighbours
  const auto& d = get_neighbours(k);

  // store data in arrays in local varibles
  const double nn = n[k];
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

  // laplacian of the chem pot
  // derivatives of velocity
  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);
  // derivatives of stress field
  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);

  // expansion/compression term
  const double expansion = dxux + dyuy;
  // shear
  const double shear     = .5*(dxuy + dyux);
  // vorticity
  const double vorticity = .5*(dxuy - dyux);
  // trace
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  // xx component of the Beris-Edwards equation
  const double Dxx = (Gamma*Hxx - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) +2*Qyx*shear -expansion));
  // zz component of the Beris-Edwards equation
  const double Dyx = (Gamma*Hyx - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*( Qyx*(expansion-traceQL) + 2*shear));

  //anisotropic friction
  const double anifricx = vx*(friction - epsilon*Qxx) - vy*epsilon*Qyx;
  const double anifricy = -epsilon*vx*Qyx + vy*(friction + epsilon*Qxx);

  //quadrupole force
  const double qfx = -zeta2*Qxx*(dxQxx+dyQyx)-zeta2*Qyx*(dxQyx-dyQxx);
  const double qfy = -zeta2*(-Qxx)*(dxQyx-dyQxx)-zeta2*(Qyx)*(dxQxx+dyQyx);

  // forcing term
  const double Fx = dxSxx + dySxy - anifricx + quad*qfx;
  const double Fy = dxSyx + dySyy - anifricy + quad*qfy;



  // calculate the equilibrium distribution fe
  const auto fe = GetEquilibriumDistribution(vx, vy, nn);

  if(first)
  {
    double Qxx_noise = 0., Qxy_noise = 0.;
    LBNode ff_noise = {0.};

    if(Q_fluct)
    {
      static const double Q_stren = sqrt(Gamma*Q_kBT);
      Qxx_noise = Q_stren*random_real();
      Qxy_noise = Q_stren*random_real();
    }

    if(u_fluct)
    {
      const auto ff_stren = sqrt(3.*nn*u_kBT*(2.*tau-1.)/tau/tau);
      ff_noise = GenerateNoiseDistribution(ff_stren);
    }

    QNxx[k] = QQxx[k] + .5*Dxx + Qxx_noise;
    QNyx[k] = QQyx[k] + .5*Dyx + Qxy_noise;
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    for(unsigned v=0; v<lbq; ++v)
    {
      fn[k][v] = ff[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)))
          + ff_noise[v];
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
    }
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;

    for(unsigned v=0; v<lbq; ++v)
      ff[k][v] = fn[k][v] + .5*((fe[v]-ff[k][v])/tau
          + w[v]*(Fx*xdir(v) + Fy*ydir(v)));
  }
}

void Nematic::UpdateQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateQuantitiesAtNode(k);
}

void Nematic::UpdateFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFieldsAtNode(k, first);
}

void Nematic::Move()
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
  swap(ff, ff_tmp);
  swap(fn, fn_tmp);
}

void Nematic::BoundaryConditionsLB()
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

void Nematic::BoundaryConditionsFields()
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
      break;
    // box
    case 3:
    case 4:
      QQxx.ApplyNeumann();
      QQyx.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      QQxx.ApplyPBC();
      QQyx.ApplyPBC();
  }
}

void Nematic::BoundaryConditionsFields2()
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

      sigmaXX.ApplyNeumann();
      sigmaYY.ApplyNeumann();
      sigmaYX.ApplyNeumann();
      sigmaXY.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      ux     .ApplyPBC();
      uy     .ApplyPBC();
      sigmaXX.ApplyPBC();
      sigmaYY.ApplyPBC();
      sigmaYX.ApplyPBC();
      sigmaXY.ApplyPBC();
  }
}

void Nematic::Step()
{
  //cout << t;
  // boundary conditions for primary fields
  BoundaryConditionsFields();
  // predictor step
  UpdateQuantities();
  // boundary conditions for the secondary
  // fields MU, u, sigma, and phi
  BoundaryConditionsFields2();
  // update all the fields
  this->UpdateFields(true);
  // boundary conditions for
  // the flow before advection
  BoundaryConditionsLB();
  // move LB particles
  Move();
  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {
    // same thing (no advection)
    BoundaryConditionsFields();
    UpdateQuantities();
    BoundaryConditionsFields2();
    this->UpdateFields(false);
  }
}

void Nematic::RuntimeChecks()
{
  // check that the sum of f is constant
  {
    double fcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        fcheck = accumulate(begin(ff[k]), end(ff[k]), fcheck);
    if(abs(ftot-fcheck)>1)
      throw error_msg("f is not conserved (", ftot, "/", fcheck, ")");
  }
}

option_list Nematic::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("Gamma", opt::value<double>(&Gamma),
     "Mobility")
    ("xi", opt::value<double>(&xi),
     "tumbling/aligning parameter")
    ("tau", opt::value<double>(&tau),
     "viscosity")
    ("rho", opt::value<double>(&rho),
     "fluid density")
    ("Q_kBT", opt::value<double>(&Q_kBT),
     "hydrodynamic fluctations strength")
    ("u_kBT", opt::value<double>(&u_kBT),
     "nematic fluctuations strength")
    ("friction", opt::value<double>(&friction),
     "friction")
    ("CC", opt::value<double>(&CC),
     "coupling constant")
    ("LL", opt::value<double>(&LL),
     "elastic constant")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("zetaC", opt::value<double>(&zetaC),
     "isotropic activity parameter")
    ("zeta2", opt::value<double>(&zeta2),
     "activity parameter (quadrupole)")
    ("quad", opt::value<int>(&quad),
     "quadrupole force on (1) or off (0)")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method")
    ("epsilon", opt::value<double>(&epsilon),
     "anisotropic friction");

  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations")
    ("no_backflow", opt::value<int>(&no_backflow),
     "turn off backflow")
    ("inits", opt::value<int>(&inits),
     "initial conditions [inits = 0 do not read, inits = 1 read]");

  return { model_options, config_options };
}
