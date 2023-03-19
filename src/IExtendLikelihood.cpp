#include "AmpGen/IExtendLikelihood.h"

#include <ostream>

#include "AmpGen/Factory.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/PolarisedSum.h"
using namespace AmpGen;

double GaussianConstraint::getVal() 
{
  return -( m_param->mean() - m_mean ) * ( m_param->mean() - m_mean ) / ( 2 * m_sigma * m_sigma );
}

void GaussianConstraint::configure( const std::string& configString, 
                                    PolarisedSum& pdf,
                                    const MinuitParameterSet& mps )
{
  auto tokens            = split( configString, ' ' );
  const std::string name = tokens[1];
  if ( tokens.size() == 4 ) {
    m_mean  = stod( tokens[2] );
    m_sigma = stod( tokens[3] );
  } else if ( tokens.size() == 3 ) {
    if ( tokens[2] == "PDG" ) {
      auto split_tokens        = split( name, '_' );
      auto particle_properties = ParticlePropertiesList::get( split_tokens[0] );
      if ( split_tokens[1] == "mass" ) {
        m_mean  = particle_properties->mass();
        m_sigma = particle_properties->mErrPlus();
      } else if ( split_tokens[1] == "width" ) {
        m_mean  = particle_properties->width();
        m_sigma = particle_properties->wErrPlus();
      } else {
        ERROR( "PDG property " << split_tokens[1] << " not recognised" );
      }
    }
  }
  m_param = mps[name];
  if ( m_param == nullptr ) {
    ERROR( "Parameter - " << name << " not found in MPS" );
  }
}

void GaussianConstraint::configure( const std::string& configString, 
                                   const MinuitParameterSet& mps )
{
    m_param = mps[configString];
    if ( m_param == nullptr ) {
        ERROR( "Parameter - " << configString << " not found in MPS" );
    }
    else {
        INFO("Adding gaussian constraint for parameter " << configString << " with mean=" <<  m_param->mean() << " and sigma=" << m_param->stepInit());
    }

    m_mean  = m_param->mean();
    m_sigma = m_param->stepInit();
}


double ExtendedNormalisation::getVal()
{
    return m_pdf->norm();
}

void ExtendedNormalisation::configure( const std::string& configString, 
                       PolarisedSum& pdf,
                      const MinuitParameterSet& mps )
{
    m_pdf       = &pdf;
}

REGISTER( IExtendLikelihood, ExtendedNormalisation );
REGISTER( IExtendLikelihood, GaussianConstraint );
