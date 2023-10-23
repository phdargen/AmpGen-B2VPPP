#include <string>
#include <complex>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Particle.h"

using namespace AmpGen;
using namespace AmpGen::fcn; 
using namespace std::complex_literals;

DEFINE_LINESHAPE( FormFactor )
{
  auto props                  = ParticlePropertiesList::get( particleName );
  Expression radius           = Parameter( particleName + "_radius", props->radius() );
  Expression mass             = Parameter( particleName + "_mass"  , props->mass() );
  const Expression q2         = make_cse( Q2( s, s1, s2 ) );
  const Expression q20        = make_cse( Q2( mass*mass, s1, s2 ) );
  int Lp = L;
  if( lineshapeModifier == "L0" ) Lp = 0;
  if( lineshapeModifier == "L1" ) Lp = 1;
  if( lineshapeModifier == "L2" ) Lp = 2;
  if( lineshapeModifier == "L3" ) Lp = 3;

  //Expression                              FormFactor = fpow(radius,Lp)*sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, Lp ) );
  Expression                              FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, Lp ) );
  if ( lineshapeModifier == "BL" )        FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, Lp ) );
  if ( lineshapeModifier == "NFF" )       FormFactor = 1; 
  if ( lineshapeModifier == "BELLE2018" ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, Lp ) );
  
  if( L != 0 ){
    ADD_DEBUG( q2      , dbexpressions );
    ADD_DEBUG( radius  , dbexpressions );
  }
  ADD_DEBUG( FormFactor, dbexpressions );
    
  return FormFactor ; 
}

DEFINE_LINESHAPE( ExpFF )
{
  auto props                  = ParticlePropertiesList::get( particleName );
  Expression radius           = Parameter( particleName + "_radius", props->radius() );
  const Expression q2         = Q2( s, s1, s2 );
  const Expression FormFactor = Exp( -q2 * radius * radius / 2. );
  if( L != 0 ){
    ADD_DEBUG( radius, dbexpressions );
    ADD_DEBUG( s1, dbexpressions );
    ADD_DEBUG( s2, dbexpressions );
    ADD_DEBUG( s, dbexpressions );
    ADD_DEBUG( q2, dbexpressions );
  }
  ADD_DEBUG( FormFactor, dbexpressions ); 
  return FormFactor ;
}

DEFINE_LINESHAPE( None ) { return Constant( 1); }

DEFINE_LINESHAPE( BW )
{
  auto s_cse = make_cse(s);
  auto props = ParticlePropertiesList::get( particleName );
  const Expression& mass       = Parameter( particleName + "_mass", props->mass() );
  const Expression& width0     = Parameter( particleName + "_width", props->width() );
  const Expression& width0_inv = Parameter( particleName + "_width_inv", 0 );
  const Expression& radius     = Parameter( particleName + "_radius", props->radius() );
  
  Expression mass_eff = mass;
  if(lineshapeModifier.find("m0eff")!= std::string::npos){
        const Expression& mmin = Parameter( particleName + "_m0eff_min", props->mass() );
        const Expression& mmax = Parameter( particleName + "_m0eff_max", props->mass() );
        const Expression x = (mass - ((mmin+mmax)/2.))/(mmax-mmin);
        const Expression tanhterm = (exp(x)-exp(-x))/(exp(x)+exp(-x));
        mass_eff = Ternary(mass <= mmin , mmin + (mmax-mmin)*(1.+tanhterm)/2., mass);
  }
  const Expression q2          = make_cse( Abs(Q2( s_cse, s1, s2 ) ) ) ;
  const Expression q20         = make_cse( Abs(Q2( mass_eff * mass_eff, s1, s2 )) );
    
  //Expression                              FormFactor = fpow(radius,L)*sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  Expression FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if( lineshapeModifier.find("BL")!= std::string::npos )        FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  if( lineshapeModifier.find("BELLE2018")!= std::string::npos ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, L ) );
  if( lineshapeModifier.find("NFF")!= std::string::npos )       FormFactor = 1;
    
  Expression runningWidth = width( s_cse, s1, s2, mass_eff, width0, radius, L, dbexpressions ) * (mass/mass_eff) + width0_inv;
    
  const Expression BW = FormFactor / ( mass * mass - s_cse  -1i * mass * runningWidth );
  const Expression kf = kFactor( mass, width0, dbexpressions );
 
  ADD_DEBUG( s_cse, dbexpressions );
  ADD_DEBUG( s1, dbexpressions );
  ADD_DEBUG( s2, dbexpressions );
  ADD_DEBUG( FormFactor, dbexpressions );
  ADD_DEBUG( runningWidth, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  ADD_DEBUG( kf, dbexpressions );
      
  return lineshapeModifier == "BELLE2018" ? BW : kf*BW;
}

DEFINE_LINESHAPE( SBW )
{
  auto props        = ParticlePropertiesList::get( particleName );
  const Expression& mass   = Parameter( particleName + "_mass", props->mass() );
  const Expression& width0 = Parameter( particleName + "_width", props->width() );
  const Expression& radius     = Parameter( particleName + "_radius", props->radius() );

  auto s_cse = make_cse(s);
  const Expression q2          = make_cse( Abs(Q2( s_cse, s1, s2 ) ) ) ;
  const Expression q20         = make_cse( Abs(Q2( mass * mass, s1, s2 )) );
    
  Expression FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, 0, L ) );
  if( lineshapeModifier.find("BL")!= std::string::npos )        FormFactor = sqrt( BlattWeisskopf( q2 * radius * radius, L ) );
  if( lineshapeModifier.find("BELLE2018")!= std::string::npos ) FormFactor = sqrt( BlattWeisskopf_Norm( q2 * radius * radius, q20 * radius * radius, L ) );
  if( lineshapeModifier.find("NFF")!= std::string::npos )       FormFactor = 1;

  const Expression kF = kFactor( mass, width0 ) ;
  const Expression BW = ( lineshapeModifier.find("sqrtS")!= std::string::npos ) ? kF / ( mass * mass - s - 1i*mass * width0 ) : kF / ( mass * mass - s - 1i*sqrt(s) * width0 );
    
  ADD_DEBUG( kF, dbexpressions );
  ADD_DEBUG( BW, dbexpressions );
  return kF * BW;
}
