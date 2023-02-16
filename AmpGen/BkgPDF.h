#ifndef AMPGEN_BKGPDF_H
#define AMPGEN_BKGPDF_H

#include <complex>
#include <iomanip>
#include <string>
#include <vector>

#include "AmpGen/EventList.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Integrator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/Event.h"
#include "AmpGen/Types.h"

namespace AmpGen
{
  class MinuitParameterSet;
  class Event;
  class EventList;
  class MinuitProxy;

  class BkgPDF
  {
    
    public:
      #if ENABLE_AVX
        using EventList_type  = EventListSIMD;
      #else
        using EventList_type  = EventList;
      #endif

      BkgPDF() = default;
      BkgPDF(const EventType&);
      void prepare(){
          m_weight.update();
          m_nCalls++;
      };
      void setEvents(EventList_type&);
      void setMC(EventList_type&);
      #if ENABLE_AVX
        void setEvents(EventList& evts){ m_ownEvents = true; setEvents( *new EventList_type(evts)) ; };
        void setMC(EventList& evts){ setMC( *new EventList_type(evts)) ; };
        double operator()(const double*, const unsigned) const;
      #endif
      //float_v operator()(const float_v*, const unsigned) const;
      real_t  operator()(const Event& evt) const {return (m_weight/m_norm) * evt[evt.size()-1] ; };
      void reset(bool resetEvents = false);
      void debug(const Event&){};
      void setWeight(MinuitProxy);
      double getWeight() const;

      real_t norm() const;
      double calculate_norm();
      real_t prob_unnormalisedNoCache(const Event& evt) const {return evt[evt.size()-1]; };

      void transferParameters(){};

      std::function<real_t(const Event&)> evaluator(const EventList_type* = nullptr) const;
      //KeyedFunctors<double(Event)> componentEvaluator(const EventList_type* = nullptr) const;
      EventType eventType() const{ return m_eventType; }
        
  private:
      size_t                        m_nCalls      = {0};
      real_t                        m_norm        = {1};
      EventList_type*               m_events      = {nullptr};
      EventList_type*               m_integrator  = {nullptr};
      MinuitParameterSet*           m_mps         = {nullptr};
      MinuitProxy                   m_weight      = {nullptr,1};
      bool                          m_verbosity   = {0};
      bool                          m_debug       = {0};
      EventType                     m_eventType;
  };
} // namespace AmpGen

#endif
