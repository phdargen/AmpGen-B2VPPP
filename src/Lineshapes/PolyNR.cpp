#include <string>

#include "AmpGen/Expression.h"
#include "AmpGen/Factory.h"
#include "AmpGen/Lineshapes.h"
#include "AmpGen/NamedParameter.h"

using namespace AmpGen;

DEFINE_GENERIC_SHAPE( PolyNR )
{
  auto p0 = *p.daughter(0);
  auto p1 = *p.daughter(1);
  auto p2 = *p.daughter(2);
  auto p01 = (p0.P() + p1.P());
  auto p02 = (p0.P() + p2.P());
  auto s01 = dot( p01, p01 );
  auto s02 = dot( p02, p02 );
  size_t degree = NamedParameter<size_t>( lineshapeModifier + "::Degree" ) +1;
  std::vector< std::vector<Parameter> > C( degree, std::vector<Parameter>(degree) );
  
  for( size_t i = 0 ; i != degree ; ++i )
  {
    for( size_t j = 0 ; j != degree ; ++j ){
      auto pname = lineshapeModifier +"_"+std::to_string(i) + "_"+std::to_string(j) ;
      C[i][j] = Parameter(pname, 0);
      if( dbexpressions != nullptr )
        dbexpressions->emplace_back( pname, C[i][j] );
    }
  }
  Expression rt; 
  for( size_t i = 0 ; i != degree; ++i )
  {
    for( size_t j = 0 ; j != degree ; ++j )
    {
      rt = rt + C[i][j] * fcn::pow(s01, i) * fcn::pow(s02, j); 
    }
  }
  return rt;
}


DEFINE_LINESHAPE( ExpNR )
{
    Expression alpha     = Parameter( particleName + "_alpha", 0 );
    Expression beta    = Parameter( particleName + "_beta", 0 );
        
    const Expression BW = Exp(-alpha * s - Constant(0,1) * beta * s);
    
    ADD_DEBUG( BW, dbexpressions );
    ADD_DEBUG( alpha, dbexpressions );
    ADD_DEBUG( beta, dbexpressions );

    return BW;
}

DEFINE_GENERIC_SHAPE( Bkg )
{
    auto p0 = *p.daughter("K+",-1);
    auto p1 = *p.daughter("pi+",-1);
    auto p2 = *p.daughter("pi-",-1);
    auto p3 = *p.daughter("mu+",-1);
    auto p4 = *p.daughter("mu-",-1);
        
    auto pR = (p3.P() + p4.P());
    auto pi = p3.P();
    auto pj = (p0.P() + p1.P() + p2.P());

    auto dot_ij = -dot(pi,pj) + dot(pi,pR) * dot(pj,pR) / dot(pR,pR);
    auto dot_i  = -dot(pi,pi) + dot(pi,pR) * dot(pi,pR) / dot(pR,pR);
    auto dot_j  = -dot(pj,pj) + dot(pj,pR) * dot(pj,pR) / dot(pR,pR);
    
    auto cosTheta = make_cse(Ternary( fcn::sqrt(dot_i * dot_j) > 1e-6, dot_ij / fcn::sqrt( dot_i * dot_j ) ,0));
    //auto cosTheta =  dot_ij / fcn::sqrt( dot_i * dot_j );
    
    Expression a     = Parameter(  "Bkg_a_cosTheta", 1 );
    Expression b     = Parameter(  "Bkg_b_cosTheta", 0.1 );
    Expression c     = Parameter(  "Bkg_c_cosTheta", -1 );
    const Expression BW = Constant(1) + a + b * cosTheta + c * fcn::pow(cosTheta,2) /(3*a+c+3);
    
    auto p012 = (p0.P() + p1.P() + p2.P());
    auto p02 = (p0.P() + p2.P());
    auto p12 = (p1.P() + p2.P());
    auto p3412 = (p3.P() + p4.P() + p1.P() + p2.P());
    auto p341 = (p3.P() + p4.P() + p1.P());
    auto p342 = (p3.P() + p4.P() + p2.P());
    auto p340 = (p3.P() + p4.P() + p0.P());
    auto p3402 = (p3.P() + p4.P() + p0.P() + p2.P());
    auto p3401 = (p3.P() + p4.P() + p0.P() + p1.P());
    auto p01 = (p0.P() + p1.P());

    auto m012 = fcn::sqrt( dot( p012, p012 ) );
    auto m02 = fcn::sqrt( dot( p02, p02 ) );
    auto m12 = fcn::sqrt( dot( p12, p12 ) );
    auto m3412 = fcn::sqrt( dot( p3412, p3412 ) );
    auto m341 = fcn::sqrt( dot( p341, p341 ) );
    auto m342 = fcn::sqrt( dot( p342, p342 ) );
    auto m340 = fcn::sqrt( dot( p340, p340 ) );
    auto m3402 = fcn::sqrt( dot( p3402, p3402 ) );
    auto m3401 = fcn::sqrt( dot( p3401, p3401 ) );
    auto m01 = fcn::sqrt( dot( p01, p01 ) );
    
    Expression a_m012     = Parameter(  "Bkg_a_m012" , 0 );
    Expression a_m02      = Parameter(  "Bkg_a_m02"  , 0 );
    Expression a_m12      = Parameter(  "Bkg_a_m12"  , 0 );
    Expression a_m3412    = Parameter(  "Bkg_a_m3412", 0 );
    Expression a_m341     = Parameter(  "Bkg_a_m341" , 0 );
    Expression a_m342     = Parameter(  "Bkg_a_m342" , 0 );
    Expression a_m340     = Parameter(  "Bkg_a_m340" , 0 );
    Expression a_m3402    = Parameter(  "Bkg_a_m3402", 0 );
    Expression a_m3401    = Parameter(  "Bkg_a_m3401", 0 );
    Expression a_m01      = Parameter(  "Bkg_a_m01"  , 0 );
    
    Expression b_m012     = Parameter(  "Bkg_b_m012" , 0 );
    Expression b_m02      = Parameter(  "Bkg_b_m02"  , 0 );
    Expression b_m12      = Parameter(  "Bkg_b_m12"  , 0 );
    Expression b_m3412    = Parameter(  "Bkg_b_m3412", 0 );
    Expression b_m341     = Parameter(  "Bkg_b_m341" , 0 );
    Expression b_m342     = Parameter(  "Bkg_b_m342" , 0 );
    Expression b_m340     = Parameter(  "Bkg_b_m340" , 0 );
    Expression b_m3402    = Parameter(  "Bkg_b_m3402", 0 );
    Expression b_m3401    = Parameter(  "Bkg_b_m3401", 0 );
    Expression b_m01      = Parameter(  "Bkg_b_m01"  , 0 );
    
    Expression c_m012     = Parameter(  "Bkg_c_m012" , 0 );
    Expression c_m02      = Parameter(  "Bkg_c_m02"  , 0 );
    Expression c_m12      = Parameter(  "Bkg_c_m12"  , 0 );
    Expression c_m3412    = Parameter(  "Bkg_c_m3412", 0 );
    Expression c_m341     = Parameter(  "Bkg_c_m341" , 0 );
    Expression c_m342     = Parameter(  "Bkg_c_m342" , 0 );
    Expression c_m340     = Parameter(  "Bkg_c_m340" , 0 );
    Expression c_m3402    = Parameter(  "Bkg_c_m3402", 0 );
    Expression c_m3401    = Parameter(  "Bkg_c_m3401", 0 );
    Expression c_m01      = Parameter(  "Bkg_c_m01"  , 0 );
    
    Expression BW2 = BW * Exp(-a_m012 * m012 - a_m02 * m02 - a_m12 * m12 - a_m3412 * m3412 - a_m341 * m341 - a_m342 * m342 - a_m340 * m340 - a_m3402 * m3402 - a_m3401 * m3401 - a_m01 * m01);
    
    Expression BW3 = BW2 * ( Constant(1) + b_m340  * m340  + c_m340  * fcn::pow(m340 ,2) )
                         * ( Constant(1) + b_m3412 * m3412 + c_m3412 * fcn::pow(m3412,2) );

    ADD_DEBUG( cosTheta, dbexpressions );
    ADD_DEBUG( BW, dbexpressions );
    ADD_DEBUG( BW2, dbexpressions );

    return BW3  ;
}
