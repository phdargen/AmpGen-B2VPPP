#include "AmpGen/BkgPDF.h"

#include <memory.h>
#include <iomanip>
#include <memory>
#include <ostream>

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/EventList.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AmpGen;

BkgPDF::BkgPDF(const EventType& type)
  : m_eventType (type)
{
}

void BkgPDF::reset(bool resetEvents){
  //m_prepareCalls=0;
  //m_lastPrint=0;
  if( resetEvents ){
    m_events = 0;
    m_integrator = 0 ;
  };
}

void BkgPDF::setEvents( EventList_type& events )
{
  reset();
  if( m_events != nullptr ) delete m_events;
  m_events = &events;
}

void BkgPDF::setMC(EventList_type& events ){
  if( m_integrator == &events ) return ;
  reset();
  m_integrator = &events;
  m_norm = calculate_norm();
}

void BkgPDF::setWeight( MinuitProxy param ){ m_weight = param; }
double BkgPDF::getWeight() const { return m_weight ; }

double BkgPDF::calculate_norm(){
  double norm = 0;
  double sample_norm = 0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction (+: norm, sample_norm)
  #endif
  for ( unsigned int i = 0; i < m_integrator->size() ; ++i ) {
    const auto& event = m_integrator->at(i);
    double v = prob_unnormalisedNoCache(event);
    double w = event.weight() / event.genPdf();
    norm        += w * v;
    sample_norm += w;
  }
  return norm / sample_norm;
}

double BkgPDF::norm() const
{
  return m_norm;
}

std::function<real_t(const Event&)> BkgPDF::evaluator(const EventList_type* ievents) const
{
    auto events = ievents == nullptr ? m_integrator : ievents;
    std::vector<double> values ( events->size() );
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (unsigned int i = 0 ; i < events->size(); ++i ){ values[i] = operator()(events->at(i)); }
    
    return [values](const Event& event) -> real_t {return values[event.index()];};
}
