#include <complex>
#include <string>
#include <vector>

#include "AmpGen/Factory.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Types.h"

using namespace AmpGen;

double LASSO::getVal() 
{
  double sum( 0 );
  for ( unsigned int i = 0; i < m_pdf->numAmps(); ++i ) {
    std::complex<double> c_i = (*m_pdf)[i].coefficient;
      sum += sqrt( std::norm(c_i) * m_pdf->norm(i, i).real() );
  }
  return m_lambda * sum;
}

void LASSO::configure( const std::string& configString,
                       PolarisedSum& pdf,
                       const MinuitParameterSet& mps )
{
  m_pdf       = &pdf;
  auto tokens = split( configString, ' ' );
  m_lambda    = stod( tokens[1] );
  INFO("Configured LASSO with lambda = " << m_lambda << " m_pdf->size() = " << m_pdf->size() << " norm = " << m_pdf->norm());
  INFO("getVal = " << getVal());

    for ( unsigned int i = 0; i < m_pdf->numAmps(); ++i ){
        INFO((*m_pdf)[i].name() << " c = " << (*m_pdf)[i].coefficient << " norm = " << m_pdf->norm(i, i).real());
    }

}

REGISTER( IExtendLikelihood, LASSO );
