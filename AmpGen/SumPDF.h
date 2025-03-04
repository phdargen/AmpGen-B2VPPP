#ifndef AMPGEN_SUMPDF_H
#define AMPGEN_SUMPDF_H

#include "AmpGen/IExtendLikelihood.h"
#include "AmpGen/MetaUtils.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/KeyedFunctors.h"
#include <tuple>

#if ENABLE_AVX
  #include "AmpGen/simd/utils.h"
#endif

namespace AmpGen
{
  class EventList; 
  class EventListSIMD;
  /** @class SumPDF
      @brief A pdf that contains one or more terms.

      A pdf with a probability of the form 
      @f[
        P(\psi) = \sum_{j} \mathcal{P}_j (\psi),
      @f] 
      where @f$ \mathcal{P}_j(\psi) @f$ are some normalised probability density functions 
      as a function of position in the phase space @f$ \psi @f$ 
      , and the sum is over the different terms, typically a signal term and then a number of background terms. 
      The pdf is also equipped with a log-likelihood of the form:
      @f[ 
        -2 \mathcal{L} = - 2 \sum_{i} \log \left( \sum_{j} \mathcal{P}_j \left(\psi_i\right) \right)
      @f]
      and the sum over @f$ i @f$ is over some dataset. 
      This combined functionality is largely historical and the two roles should be separated at some point in the future. 
      The sum is variadically unrolled at compile time, i.e. the wrapper
      is the same for 1..N pdfs. The unrolling should be properly inlined,
      hence N can be reasonably large with out afflicting either
      compile time or binary size. It isn't primarily used as PDF, as its primary function is 
      as a likelihood via function getVal(). 
      Typically constructed using either the make_pdf helper function or make_likelihood helper function.  */
  template <class eventListType, class... pdfTypes>
  class SumPDF
  {
  private:
    typedef typename eventListType::value_type eventValueType; ///< The value type stored in the eventListType
    std::tuple<pdfTypes...> m_pdfs;                            ///< The tuple of probability density functions
    const eventListType*    m_events = {nullptr};              ///< The event list to evaluate likelihoods on

  public:
    /// Default Constructor
    SumPDF() = default; 

    /// Constructor from a set of PDF functions
    SumPDF( const pdfTypes&... pdfs ) : m_pdfs( std::tuple<pdfTypes...>( pdfs... ) ) {}

    /// Returns negative twice the log-likelihood for this PDF and the given dataset.     

    double getVal()
    {
      if constexpr( std::is_same<eventListType,EventList>::value )
      {
        double LL = 0;
        double sumW = 0;
        double sumW2 = 0;
        for_each( m_pdfs, []( auto& f ) { f.prepare(); } );
        #pragma omp parallel for reduction( +: LL,sumW,sumW2 )
        for ( unsigned int i = 0; i < m_events->size(); ++i ) {
          auto prob = ((*this))(( *m_events)[i] );
          auto w = (*m_events)[i].weight();
          //w = (abs(w)<3) ? w : 0;
          prob = (w < 0 && prob < 0.001) ? 0.001 : prob;
          //if(w < 0 prob<0.0001)INFO(prob);
          auto val = w*log(prob);
          //INFO(val);
          //LL += std::isnan(val) ? 0. : val;
          LL += val;
          sumW += w;
          sumW2 += w*w;
        }
        //double N = (double)m_events->size();
        //for_each( m_pdfs, [&LL,&N]( auto& f ) { LL += N * log(f.norm()) - f.norm() ; } );
        //INFO(sumW/sumW2);
        return -2 * sumW/sumW2 * LL;
//        return -2 * sumW/sumW2 * LL;
      }
      #if ENABLE_AVX 
      if constexpr( std::is_same<eventListType, EventListSIMD>::value )
      {
        float_v LL = 0.f;
        for_each( m_pdfs, []( auto& f ) { f.prepare(); } );
        #pragma omp parallel for reduction( +: LL )
        for ( size_t block = 0; block < m_events->nBlocks(); ++block ) 
        {
          LL += m_events->weight(block) * AVX::log(this->operator()(m_events->block(block), block));
        }
        return -2 * utils::sum_elements(LL);
      }
      #endif
    } 
    /// Returns the probability for the given event. 
    float_v operator()( const float_v* evt , const unsigned block)
    {
      float_v prob = 0.f;
      for_each( this->m_pdfs, [&prob, &evt,block]( const auto& f ) { prob += f(evt, block); } );
      return prob;
    }
    /// Returns the probability for the given event. 
    double operator()( const eventValueType& evt )
    {
      double prob = 0;
      for_each( this->m_pdfs, [&prob, &evt]( const auto& f ) { prob += f(evt); } );
      return prob;
    }

    /// Sets the events to be summed over in the likelihood
    void setEvents( eventListType& events )
    {
      m_events = &events;
      for_each( m_pdfs, [&events]( auto& f ) { f.setEvents( events ); } );
    }
    
    /// Returns the number of PDFs contained by this function 
    std::size_t nPDFs() const { return sizeof...(pdfTypes); }

    /// Returns the tuple of PDFs used by this function
    std::tuple<pdfTypes...> pdfs() const { return m_pdfs; }

    std::function<double(const eventValueType&)> evaluator(const eventListType* events) const
    {     
      std::vector<double> values( events->size() );
      for_each( this->m_pdfs, [events, &values](const auto& pdf ) mutable { 
        auto eval = pdf.evaluator(events);
        for( unsigned i = 0; i != events->size(); ++i ) values[i] += eval( events->at(i) ); 
      } );
      return arrayToFunctor<double, typename eventListType::value_type>(values);
    }
    KeyedFunctors<double(eventValueType)> componentEvaluator(const eventListType* events) const
    { 
      KeyedFunctors<double(eventValueType)> view;
      for_each( this->m_pdfs, [&view, &events]( const auto& pdf) mutable { 
        auto eval = pdf.evaluator(events);
        view.add([eval](const auto& event){ return eval(event) ; } , type_string(pdf), "" ); 
      } );
      return view; 
    }
  };

  /** @function make_pdf 

    Usage is 
    \code{cpp}
      auto pdf = make_pdf( signal, bkg1, ... );
    \endcode
    which returns a PDF that is the sum of the signal and bkg1 etc. The sum is also equipped with a likelihood, which can be used by setting 
    the data for the PDF.
    Therefore, named SumPDF, it is useful to use this function to get a likelihood for a PDF containing a single term (i.e. signal or background only).  
    */
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_pdf( pdfTypes&&... pdfs )
  {
    //return SumPDF<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    return SumPDF<eventListType, pdfTypes...>( pdfs... );
  }
  
  template <class eventListType = EventList, class... pdfTypes> 
  auto make_likelihood( eventListType& events, pdfTypes&&... pdfs )
  {
    auto rt = SumPDF<eventListType, pdfTypes...>( std::forward<pdfTypes>( pdfs )... );
    rt.setEvents(events);
    return rt; 
  }
} // namespace AmpGen

#endif
