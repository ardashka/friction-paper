#include "header.hpp"
#include "models/incompdrypolarlyo.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"

using namespace std;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

IncompDryPolarLyo::IncompDryPolarLyo(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}
IncompDryPolarLyo::IncompDryPolarLyo(unsigned LX, unsigned LY, unsigned BC, GridType Type)
  : Model(LX, LY, BC, Type)
{}

void IncompDryPolarLyo::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  Px.SetSize(LX, LY, Type);
  Py.SetSize(LX, LY, Type);
  PNx.SetSize(LX, LY, Type);
  PNy.SetSize(LX, LY, Type);
  phi.SetSize(LX, LY, Type);
  phi_tmp.SetSize(LX, LY, Type);
  phn.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  ux_phi.SetSize(LX, LY, Type);
  uy_phi.SetSize(LX, LY, Type);
  Hx.SetSize(LX, LY, Type);
  Hy.SetSize(LX, LY, Type);
  MU.SetSize(LX, LY, Type);
  dxPx.SetSize(LX, LY, Type);
  dyPx.SetSize(LX, LY, Type);
  dxPy.SetSize(LX, LY, Type);
  dyPy.SetSize(LX, LY, Type);
  sigmaXX.SetSize(LX, LY, Type);
  sigmaYY.SetSize(LX, LY, Type);
  sigmaYX.SetSize(LX, LY, Type);
  sigmaXY.SetSize(LX, LY, Type);

  FFx.SetSize(LX, LY, Type);
  FFy.SetSize(LX, LY, Type);
  pres.SetSize(LX, LY, Type);
  pres_tmp.SetSize(LX, LY, Type);

  P_fluct = (P_kBT!=0);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");
}

void IncompDryPolarLyo::ConfigureAtNode(unsigned k)
{
  double Order = 0;
  double theta;
//  double xtemp,ytemp;
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);
//  xtemp=x;
//  ytemp=y;
  if(init_config=="circle")
  {
    if(pow(diff(x, LX/2), 2) + pow(diff(y, LY/2), 2) <= radius*radius)
      Order = init_order;
  }
  else if(init_config=="square")
  {
    if (diff(LY/2, y) < level/2 && diff(LX/2, x) < level/2)
      Order = init_order;
  }
  else if(init_config=="stripe")
  {
    if(diff(LY/2, y) < level/2) Order = init_order;
  }
  else if(init_config=="half")
  {
    if(y < level) Order = init_order;
  }
  else if(init_config=="boxatwall")
  {
    if (BC==201 or BC==501 or BC==4)
    {
      if(x < level and diff(LY/2, y) < level/2) Order = init_order;
    }
    else
    {
      if(y < level and diff(LX/2, x) < level) Order = init_order;
    }
  }
  else if(init_config=="cuttingatwall")
  {
    if (BC==201 or BC==501 or BC==4)
    {
      if(pow(x, 2) + pow(wrap(diff(y, LY/2+int(radius*1.5)), LY), 2) <= radius*radius) Order = init_order;
      if(pow(x, 2) + pow(wrap(diff(y, LY/2-int(radius*1.5)), LY), 2) <= radius*radius) Order = init_order;
    }
    else
    {
      if(y < level and diff(LX/2, x) < level) Order = init_order;
    }
  }
  else if(init_config=="circleatwall")
  {
    if (BC==201 or BC==501)
    {
      if(pow(x, 2) + pow(wrap(diff(y, LY/2), LY), 2) <= radius*radius) Order = init_order;
    }
    else
    {
      if(pow(wrap(diff(x, LX/2), LX),2) + pow(diff(y, 0*LY), 2) <= radius*radius) Order = init_order;
    }
  }
    else if(init_config=="wettedwall")
  {
    if (BC==201 or BC==501)
    {
      if(x<level) Order = init_order;
    }
    else
    {
      if(y<level) Order = init_order;
    }
  }
  else
    throw error_msg("error: initial configuration '", init_config, "' unknown.");


  theta   = angle + noise*M_PI*(random_real() - .5);
	
  Px[k] = Order*(cos(theta));
  Py[k] = Order*(sin(theta));
  phi[k]  = conc + noise*(random_real() - .5);
  totalphi += phi[k];
  pres_tmp[k]=0;
  pres[k]=0;
  ux[k] = uy[k] = ux_phi[k] = uy_phi[k] = 0;
  // compute totals for later checks
  ptot += phi[k];
}
void IncompDryPolarLyo::Configure()
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

void IncompDryPolarLyo::UpdatePolarQuantitiesAtNode(unsigned k)
{
  const auto& d = get_neighbours(k);  

  const double px = Px[k];
  const double py = Py[k]; 
  const double ps = px*px+py*py;
  const double p4 = px*px*px*px+py*py*py*py+2*px*px*py*py;

  const double dxpx = derivX(Px, d, sB);
  const double dypx = derivY(Px, d, sB);
  const double dxpy = derivX(Py, d, sB);
  const double dypy = derivY(Py, d, sB);

  const double p        = phi[k];
  const double dxPhi    = derivX   (phi, d, sB);
  const double dyPhi    = derivY   (phi, d, sB);
  const double del2p    = laplacian(phi, d, sD);

  const double term = 1 - p4;
  const double hx = CC*term*4*ps*px + Kp*laplacian(Px, d, sD) // CC*term*4*ps*px for p4
	  	  + Kn*( 1.*ps*laplacian(Px, d, sD)
         	  	+2.*px*(dxpx*dxpx+dypx*dypx-dxpy*dxpy-dypy*dypy)
         	  	+4.*py*(dxpx*dxpy+dypx*dypy))
		  - nu*dxPhi;
  const double hy = CC*term*4*ps*py + Kp*laplacian(Py, d, sD) // CC*term*4*ps*py for p4
	  	  + Kn*( 1.*ps*laplacian(Py, d, sD)
         	        +2.*py*(dxpy*dxpy+dypy*dypy-dxpx*dxpx-dypx*dypx)
         		+4.*px*(dxpx*dxpy+dypx*dypy))
		  - nu*dyPhi;

  const double mu = (surftension_on? Aphi*p*(p-1.)*(p+1.) + CC*term - Kphi*del2p : Aphi*(p-1.) + 0.0*CC*term - Kphi*del2p);// no CC term because term = 1 - ... instead of p - ...

  // computation of sigma...
  const double sigmaB = (surftension_on? .5*Aphi*p*p*(-1.+.5*p*p) - mu*p : .5*Aphi*(p-1.)*(p-1.) - mu*(p-1.) ) + (backflow_on? .5*CC*term*term : 0);
  const double sigmaF = .5*Kphi*(dyPhi*dyPhi-dxPhi*dxPhi);
  const double sigmaS = - Kphi*dxPhi*dyPhi;

  // backflow implemeted as of Giomi & Marchetti, Soft Matter (2012): https://doi.org/10.1039/C1SM06077E with \bar{lambda}=\lambda
  const double sigmaxx = -(p > 1.0 ? 1 : 0.5)*.5*zeta*(px*px-py*py) + (p > 1.0 ? 1 : 0.5)*(-zetap*(dxpx)) + ( ( backflow_on? -2.*xi*px*hx-xi*py*hy : 0) ) + sigmaB + sigmaF;
  const double sigmayy = -(p > 1.0 ? 1 : 0.5)*.5*zeta*(py*py-px*px) + (p > 1.0 ? 1 : 0.5)*(-zetap*(dypy)) + ( ( backflow_on? -2.*xi*py*hy-xi*px*hx : 0) ) + sigmaB - sigmaF;
  const double sigmaxy = -(p > 1.0 ? 1 : 0.5)*zeta*px*py + (p > 1.0 ? 1 : 0.5)*(-.5*zetap*(dxpy + dypx)) + ( ( backflow_on? .5*(1.-xi)*px*hy - .5*(1.+xi)*py*hx :0) ) + sigmaS;
  const double sigmayx = -(p > 1.0 ? 1 : 0.5)*zeta*py*px + (p > 1.0 ? 1 : 0.5)*(-.5*zetap*(dypx + dxpy)) + ( ( backflow_on? .5*(1.-xi)*py*hx - .5*(1.+xi)*px*hy :0) ) + sigmaS;

  // transfer to arrays  
  Hx[k]    =  hx;
  Hy[k]    =  hy;
  dxPx[k]  =  dxpx;
  dxPy[k]  =  dxpy;
  dyPx[k]  =  dypx;
  dyPy[k]  =  dypy;
  sigmaXX[k] =  sigmaxx;
  sigmaYY[k] =  sigmayy;
  sigmaXY[k] =  sigmaxy;
  sigmaYX[k] =  sigmayx;
  MU[k]      =  mu;

}
void IncompDryPolarLyo::UpdatePolarQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdatePolarQuantitiesAtNode(k);
}

void IncompDryPolarLyo::UpdateFluidQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);
  const double Fx = dxSxx + dySxy;
  const double Fy = dxSyx + dySyy;
  
  FFx[k] = Fx;
  FFy[k] = Fy;

  const double p   = phi[k];
  const double px = Px[k];
  const double py = Py[k];
  const double hx = Hx[k];
  const double hy = Hy[k];

  // compute velocities
  const double vx = (FFx[k] + (preinit_flag? 0 : (p > 1.0 ? 1 : 0.5)*alpha*px) - beta*hx)/friction;
  const double vy = (FFy[k] + (preinit_flag? 0 : (p > 1.0 ? 1 : 0.5)*alpha*py) - beta*hy)/friction;

  // transfer to arrays
  ux[k]      =  vx;
  uy[k]      =  vy;
  ux_phi[k]  =  1.0*vx*p + (preinit_flag? 0 : (p > 1.0 ? 1 : 0.5)*Vphi*px*p);
  uy_phi[k]  =  1.0*vy*p + (preinit_flag? 0 : (p > 1.0 ? 1 : 0.5)*Vphi*py*p);
}
void IncompDryPolarLyo::UpdateFluidQuantities()
{
  double sum = 0;

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){
    sum = sum + phi[k];
    UpdateFluidQuantitiesAtNode(k);
  }
  countphi = sum;
    Incompressibility();
}
void IncompDryPolarLyo::Incompressibility()
{
    //Make sure we always do the loop at least once
    pres_tmp_sum=pres_sum+1.0;
    int loop_pres=0;
       // cout<<"loop "<<loop_pres<<" pres_tmp_sum "<<pres_tmp_sum<<" pres_sum "<<pres_sum<<endl;

    while (abs(pres_sum-pres_tmp_sum)>accuracy_inc){
           // loop_pres=loop_pres+1;
        pres_tmp_sum=0;
        pres_sum=0;
        pres_tmp=pres;
        for(unsigned k=0; k<DomainSize; ++k){
            pres_tmp_sum=pres_tmp_sum+abs(pres_tmp[k]);
            UpdatePressureAtNode(k);
            pres_sum=pres_sum+abs(pres[k]);
        }
            //cout<<"loop "<<loop_pres<<" pres_tmp_sum "<<pres_tmp_sum<<" pres_sum "<<pres_sum<<endl;

    }

    for(unsigned k=0; k<DomainSize; ++k){
        IncompressibilityAtNode(k);
    }
}

void IncompDryPolarLyo::UpdatePressureAtNode(unsigned k)
{
        const auto& d = get_neighbours(k);

        const double pres_loc = pres[k];

        const double del2pres_tmploc  = laplacian(pres,  d, sD);
        const double dxux    = derivX   (ux,  d, sB);
        const double dyuy    = derivY   (uy,  d, sB);
        const double Source_pres = dxux + dyuy;
            //Numbers dependent on the dS.
        //const double pres_tmploc=3/10*(del2pres_tmploc+10/3*pres_loc-Source_pres);
        const double pres_tmploc=3.0 / 10.0 * (del2pres_tmploc+10.0 / 3.0 *pres_loc-Source_pres);
        pres[k]  =  pres_tmploc;
}

void IncompDryPolarLyo::IncompressibilityAtNode(unsigned k)
{
        // array placeholders for current node
        const auto& d = get_neighbours(k);
        const double presux    = derivX   (pres,  d, sB);
        const double presuy    = derivY   (pres,  d, sB);
        ux[k]=ux[k]-presux;
        uy[k]=uy[k]-presuy;
}

void IncompDryPolarLyo::UpdatePolarFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  
  const double vx = ux[k];
  const double vy = uy[k];
  const double vn = sqrt(vx*vx + vy*vy);

  const double p   = phi[k];
  const double px = Px[k];
  const double py = Py[k];
  const double hx = Hx[k];
  const double hy = Hy[k];

  const double dxpx = dxPx[k];
  const double dypx = dyPx[k];
  const double dxpy = dxPy[k];
  const double dypy = dyPy[k];

  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  // corrections to the polarisation, the beta term is alignment to velocity after Maitra et al. PRL (2020) DOI: 10.1103/PhysRevLett.124.028002
  const double Dx = hx/gamma 
	  	    - (vx + (p > 1.0 ? 1 : 0.5)*Vp*px)*dxpx - (vy + (p > 1.0 ? 1 : 0.5)*Vp*py)*dypx
                    + (flow_alignment ? (xi*(dxux*px + .5*(dyux+dxuy)*py) + .5*(dyux-dxuy)*py ) : 0 )
		    + beta*vx/vn;
  const double Dy = hy/gamma 
	  	    - (vx + (p > 1.0 ? 1 : 0.5)*Vp*px)*dxpy - (vy + (p > 1.0 ? 1 : 0.5)*Vp*py)*dypy
                    + (flow_alignment ? (xi*(dyuy*py + .5*(dxuy+dyux)*px) + .5*(dxuy-dyux)*px ) : 0 )
		    + beta*vy/vn;

  if(first)
  {
    PNx[k] = Px[k] + .5*Dx;
    PNy[k] = Py[k] + .5*Dy;

    if(P_fluct)
    {
      if(P_fluc_mode==1)
      {
        double rtheta, S, theta, thetan;
        static const double P_stren = sqrt(P_kBT/gamma);
        rtheta = P_stren*2*M_PI*random_real(-1, 1);

        S = sqrt(PNx[k]*PNx[k] + PNy[k]*PNy[k]);
        theta = atan2(PNy[k]/S, PNx[k]/S);
        thetan = theta + rtheta;
        PNx[k] = S * cos(thetan);
        PNy[k] = S * sin(thetan);

      }
      else if(P_fluc_mode==2)
      {
        static const double P_stren = sqrt(P_kBT/gamma);

        const double Px_noise = P_stren * random_real(-1, 1);
        const double Py_noise = P_stren * random_real(-1, 1);

        PNx[k] = PNx[k] + Px_noise;
        PNy[k] = PNy[k] + Py_noise;
      } 
      else
      {
        cout << "ERROR! P_fluc_mode neither 1 or 2. It is: "<< P_fluc_mode << endl;
        throw error_msg("P_fluc_mode neither 1 or 2. It is:", P_fluc_mode);
      }
    }

    Px[k] = PNx[k] + .5*Dx;
    Py[k] = PNy[k] + .5*Dy;
  }
  else
  {
    Px[k] = PNx[k] + .5*Dx;
    Py[k] = PNy[k] + .5*Dy;
    
  }
}
void IncompDryPolarLyo::UpdatePolarFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdatePolarFieldsAtNode(k, first);
}

void IncompDryPolarLyo::UpdateFluidFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  
  //const double vx = ux[k];
  //const double vy = uy[k];
  //const double px = Px[k];
  //const double py = Py[k];
  //const double hx = Hx[k];
  //const double hy = Hy[k];
  //const double p   = phi[k];

  const double del2mu = laplacian(MU, d, sD);
  const double pFlux = derivX(ux_phi, d, sB) + derivY(uy_phi, d, sB);
  const double Dp = (preinit_flag? 0 : GammaPhi*del2mu - pFlux - ( conserve_phi ? (countphi-totalphi)/DomainSize : 0 ) );


  if (first){
    phn[k]     = phi[k]  + .5*Dp;
    phi_tmp[k] = phi[k]  +    Dp;
  }
  else{
    phi_tmp[k] = phn[k]  + .5*Dp;
  }

}
void IncompDryPolarLyo::UpdateFluidFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFluidFieldsAtNode(k, first);

  swap(phi.get_data(), phi_tmp.get_data());
}

void IncompDryPolarLyo::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      Px.ApplyNeumannChannel();
      Py.ApplyNeumannChannel();
      phi .ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      Px.ApplyNeumann();
      Py.ApplyNeumann();
      phi .ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      Px.ApplyPBC();
      Py.ApplyPBC();
      phi .ApplyPBC();
  }
}
void IncompDryPolarLyo::BoundaryConditionsFields2()
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

void IncompDryPolarLyo::Step()
{
  // boundary conditions for primary fields
  BoundaryConditionsFields();
  // predictor step
  UpdatePolarQuantities();
  UpdateFluidQuantities();
  
  BoundaryConditionsFields2();

  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {    
    BoundaryConditionsFields();
    UpdatePolarQuantities();
    UpdateFluidQuantities();
    BoundaryConditionsFields2();
    this->UpdatePolarFields(true); 
    this->UpdateFluidFields(true);
  }
}

void IncompDryPolarLyo::RuntimeChecks()
{
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

option_list IncompDryPolarLyo::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("gamma", opt::value<double>(&gamma),
     "Q-tensor mobility")
    ("GammaPhi", opt::value<double>(&GammaPhi),
     "binary mobility")
    ("xi", opt::value<double>(&xi),
     "tumbling/aligning parameter")
     ("beta", opt::value<double>(&beta),
     "alignment to velocity parameter")
    ("friction", opt::value<double>(&friction),
     "friction from confinement")
    ("Aphi", opt::value<double>(&Aphi),
     "binary fluid bulk constant")
    ("CC", opt::value<double>(&CC),
     "coupling constant")
    ("Kn", opt::value<double>(&Kn),
     "nematic elastic constant")
    ("Kp", opt::value<double>(&Kp),
     "polar elastic constant")
    ("Kphi", opt::value<double>(&Kphi),
     "binary gradient constant")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("zetap", opt::value<double>(&zetap),
     "polar active stress")
    ("Vphi", opt::value<double>(&Vphi),
     "self-advection parameter in phi")
    ("Vp", opt::value<double>(&Vp),
     "self-advection parameter")
    ("nu", opt::value<double>(&nu),
     "polar alignment with density gradient")
    ("alpha", opt::value<double>(&alpha),
     "monopole activity parameter")
    ("P_kBT", opt::value<double>(&P_kBT),
     "strength of polar fluctuations")
    ("P_fluc_mode", opt::value<int>(&P_fluc_mode),
     "type of P fluctuations, 1: angle fluctuations, 2: P fluctuations")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method")
    ("backflow_on", opt::value<bool>(&backflow_on),
     "Backflow flag")
    ("surftension_on", opt::value<bool>(&surftension_on),
     "Surface tension flag")
    ("n_preinit", opt::value<int>(&n_preinit),
     "number of preinitialization steps")
    ("flow_alignment", opt::value<bool>(&flow_alignment),
     "flow_alignment flag");

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
     "size of initial variations")
    ("initial-order", opt::value<double>(&init_order),
     "initial order of the polarisation field");

  return { model_options, config_options };
}
