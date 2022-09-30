#ifndef AMPGEN_STORE_H
#define AMPGEN_STORE_H

#include "AmpGen/simd/utils.h"
#include "AmpGen/EventList.h"
#ifdef _OPENMP 
#include <omp.h>
#endif

namespace AmpGen {

  enum Alignment {
    SoA, AoS
  };

  template <typename stored_type, Alignment align = SoA> class Store 
  {
    public:
      Store( const size_t& nEntries=0, const size_t& nFields=0) : 
        m_nEntries(nEntries), 
        m_nBlocks(utils::aligned_size<stored_type>( nEntries ) / utils::size<stored_type>::value ),
        m_nFields(nFields),
        m_store(m_nBlocks * m_nFields) {}
      
      template <typename functor_type> 
      void addFunctor( const functor_type& functor )
      {
        if( m_index.count(functor.name()) != 0 ) return; 
        unsigned vsize = functor.returnTypeSize() / sizeof(stored_type);
        DEBUG("Registering: " << functor.name() << " field: " << m_nFields << " " << functor.returnTypeSize() << " / " << sizeof(stored_type)  ); 
        std::vector<unsigned> offsets ( vsize ) ;
        for ( unsigned i = 0 ; i != vsize; ++i ) offsets[i] = m_nFields + i; 
        m_index[ functor.name() ] = offsets; 
        m_nFields += vsize;
      }
      template <typename functor_type> void allocate( const size_t& nEntries, const std::vector<functor_type>& functors)
      {
        for(const auto& functor : functors) addFunctor( functor );
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_store.resize(m_nBlocks * m_nFields); 
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type > 
        void allocate( const size_t& nEntries, const functor_type& functor)
      {
        addFunctor(functor);
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_store.resize(m_nBlocks * m_nFields); 
      }
       
      template <typename functor_type> Store( const size_t& nEntries, const std::vector<functor_type>& functors)
      {
        allocate(nEntries, functors);
      }
      template <typename functor_type, typename = typename std::enable_if<!std::is_integral<functor_type>::value>::type>
      Store( const size_t& nEntries, const functor_type& functor)
      {
        allocate( nEntries, {functor});
      }

      inline stored_type operator[]( const size_t& index ) const { return m_store[index]; }
      inline stored_type& operator[]( const size_t& index ) { return m_store[index]; }
      template <typename T> unsigned find( const T& t ) const { return m_index.find( t.name() )->second[0]; }

      inline size_t size()         const { return m_nEntries; }
      inline size_t size_raw()     const { return m_store.size(); }
      inline size_t nBlocks()      const { return m_nBlocks; }
      inline size_t nFields()      const { return m_nFields; }
      inline size_t aligned_size() const { return m_nBlocks * utils::size<stored_type>::value ; }
      inline const stored_type& operator()(const size_t& index, const size_t& field) const
      {
        if constexpr( align == Alignment::SoA ) return m_store[ field * m_nBlocks + index] ; 
        else return m_store[index*m_nFields+field]; 
      }
      template <typename return_type> 
      inline const return_type get(const size_t& index, const size_t& field ) const 
      {
        return utils::at( operator()( index / utils::size<stored_type>::value, field ), index % utils::size<stored_type>::value );
      }
      inline const stored_type* data() const { return m_store.data(); }
      inline stored_type* data() { return m_store.data() ;}
      inline stored_type& operator()(const size_t& index, const size_t& field)
      {
        if constexpr( align == Alignment::SoA ) return m_store[ field * m_nBlocks + index] ; 
        else return m_store[index*m_nFields+field]; 
      }
      
      void resize(const size_t& nEntries, const size_t& nFields )
      {
        m_nEntries = nEntries;
        m_nBlocks  = utils::aligned_size<stored_type>(nEntries)/utils::size<stored_type>::value;
        m_nFields  = nFields;
        m_store.resize(m_nBlocks * m_nFields); 
        m_index.clear();
      }
      void clear() { m_store.clear(); m_index.clear() ; }
      void store( const size_t& event0, const unsigned* index, const stored_type* item, const unsigned N = 1 )
      {
        for( unsigned i = 0 ; i != N; ++i )
          (*this)(event0, index[i] ) = item[i]; 
      }

      template <typename functor_type, typename input_type> void update(const Store<input_type, Alignment::AoS>& is, const functor_type& fcn)
      {
        auto f = m_index.find( fcn.name() ); 
        if( f == m_index.end() ) FATAL("Expression: " << fcn.name() << " is not registed");
        DEBUG("Updating: "  << fcn.name() 
           << " index: {"    << vectorToString( f->second, ",")  <<"}"
           << " on store: " << is.size() 
           << " blocks: "   << is.nBlocks() 
           << " fields: "   << is.nFields () ); 
        auto stagger      = align == Alignment::AoS ? 1 : m_nBlocks;
        auto fieldStagger = align == Alignment::AoS ? m_nFields : 1;
        std::vector<size_t> offsets( f->second.size() );
        for( unsigned int i = 0 ; i != offsets.size(); ++i ) offsets[i] = f->second[i] * stagger; 
        
        if constexpr( std::is_same< typename functor_type::return_type, void >::value )
        {
          fcn.batch(aligned_size(), 
                    is.nFields(), 
                    fieldStagger,
                    nullptr, 
                    m_store.data(), 
                    offsets.data(),
                    fcn.externBuffer().data(), 
                    is.data());  
        }
        else 
        {
          fcn.batch(aligned_size(), is.nFields(), fieldStagger, m_store.data() + f->second[0] *stagger, fcn.externBuffer().data(), is.data() ); 
        }
      }
      template <typename functor_type> void update( const EventList& events, const functor_type& fcn )
      {  
        auto f = m_index.find( fcn.name() ); 
        if( f == m_index.end() ) FATAL("Expression: " << fcn.name() << " is not registed");
        auto p0 = f->second[0];
        auto s  = f->second.size(); 
        auto stagger      = align == Alignment::AoS ? 1 : m_nBlocks;
        std::vector<size_t> offsets( s );
        std::iota( offsets.begin(), offsets.end(), 0 );
        DEBUG("Updating: " << fcn.name() << " index = " << p0 << " size_of = " << s << " on store: " << size() << " blocks = " << nBlocks() << " fields = " << nFields () ); 
        if constexpr( std::is_same< typename functor_type::return_type, void >::value ) 
        {          
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t evt = 0; evt < events.size(); ++evt )
          { 
            std::vector<stored_type> buffer(s);
            fcn(buffer.data(), offsets.data(), fcn.externBuffer().data(), events[evt].address() ); 
            store(evt, f->second.data(), buffer.data(), s ); 
          }
        }
        else {
          #ifdef _OPENMP
          #pragma omp parallel for
          #endif
          for ( size_t evt = 0; evt < events.size(); ++evt )
          {
            auto tmp = fcn( events[evt].address() );
            store( evt, f->second.data(), &tmp, s);
          }       
        }
      }

    private:
      size_t m_nEntries{0};                                 ///< Number of entries, i.e. number of events 
      size_t m_nBlocks {0};                                 ///< Number of blocks, i.e. number of entries aligned to the size, divided by block size. 
      size_t m_nFields {0};                                 ///< Number of fields per entry 
      std::vector<stored_type> m_store;                     ///< Actual store of values  
      std::map<std::string, std::vector<unsigned>> m_index; ///< Index between name of functors and the stored values 
  };
}
#if DEBUG_LEVEL ==1 
using aos_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::AoS>;
using soa_store = AmpGen::Store<AmpGen::complex_v, AmpGen::Alignment::SoA>;

ENABLE_DEBUG(aos_store)
ENABLE_DEBUG(soa_store)
#endif

#endif
