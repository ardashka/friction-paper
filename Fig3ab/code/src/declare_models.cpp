// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers
#include "models/minimal.hpp"
#include "models/lyotropic.hpp"
#include "models/lyoadvection.hpp"
#include "models/drylyotropic.hpp"
#include "models/nematic.hpp"
#include "models/polarlyotropic.hpp"
#include "models/drypolarlyotropic.hpp"
#include "models/fakepolarity.hpp"
#include "models/incompdrypolarlyo.hpp"
//#include "models/activemodelhplus.hpp"

void DeclareModels()
{
  declare_model<Minimal>(
     "minimal",
      "This is just an example model showing a minimal implementation. "
     "This does exactly nothing (and not particularly fast)."
    );

  declare_model<Lyotropic>(
      "lyotropic",
      "Biphasic, lyotropic, nematic model as presented as described in "
      "10.1103/PhysRevLett.113.248303. We refer the user to this reference "
      "for further information."
      );

    declare_model<Lyoadvection>(
      "lyoadvection",
      "Biphasic, lyotropic, nematic model as presented as described in "
      "10.1103/PhysRevLett.113.248303. We refer the user to this reference "
      "for further information. Modified to account for advecting density field"
      );

    declare_model<DryLyotropic>(
      "drylyotropic",
      "Biphasic, lyotropic, nematic model as presented as described in "
      "10.1103/PhysRevLett.113.248303. We refer the user to this reference "
      "for further information. Extended to dry limit and with polarity "
      "alignment to the velocity. "
      );

  declare_model<Nematic>(
      "nematic",
      "Pure nematic model with LdG free energy."
      ); 
      declare_model<FakePolarity>(
      "fakepolarity",
      "Lyotropic model with fake polarity."
      );  

  declare_model<PolarLyotropic>(
      "polarlyotropic",
      "Polar model combined with the two-phase lyotropic formulation similar to the lyotropic model"
      );  

    declare_model<DryPolarLyotropic>(
      "drypolarlyotropic",
      "Polar model combined with the two-phase lyotropic formulation similar to the lyotropic model"
      );

    declare_model<IncompDryPolarLyo>(
      "incompdrypolarlyo",
      "Polar model, incompressible. For segregation paper"
      );

   //declare_model<Polar>(
   //   "polar",
   //   "Polar model allowing for also nematic defects"
   //   );

   //declare_model<ActiveModelHPlus>(
   //   "activemodelhplus",
   //   "Active model H+,  10.1103/PhysRevLett.115.188302 with additional activity in the current density (J) "
   //   );

}
