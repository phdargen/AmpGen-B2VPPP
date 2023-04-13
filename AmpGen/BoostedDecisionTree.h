#include "AmpGen/MsgService.h"
#include "AmpGen/BinDT.h"
#include <algorithm>
#include "AmpGen/ProfileClock.h"
#include <random>

struct Store {
  Store( const size_t& dim ) : m_dim(dim) {} ;
  template < class CONTAINER, class FUNCTOR > Store( const CONTAINER& container , const FUNCTOR& functor, const unsigned int& dim ) : m_dim(dim) {
    for( auto& item : container ) push( functor( item ) );
  }
  void push( const std::vector<double>& f ){
    store.insert( store.end(), f.begin(), f.end() ); 
  }
  size_t m_dim;
  std::vector<double> store;
  std::vector<double*> make_pointers() {
    std::vector<double*> rt( store.size() / ( m_dim + 1) );
    for( unsigned int i = 0 ; i < store.size()/(m_dim+1) ; ++i ) rt[i] = &( store[i*(m_dim+1)] );
    return rt;
  }
  const double* evt ( const unsigned int& index ) const { return & ( store[(m_dim+1)*index]  ); }
  double weight( const unsigned int& index ) const { return store[( m_dim +1 ) * index  + m_dim ] ; }
  void setWeight( const unsigned int& index , const double& weight ){
    store[( m_dim +1 ) * index  + m_dim ] = weight ; 
  } 
  size_t n() const { return store.size()/(m_dim+1) ; }
  void print( const unsigned int& index){
    auto E = evt(index);
    for( unsigned int i = 0 ; i < (m_dim+1) ; ++i ){
      INFO( "store["<<i<<"] = " << E[i] );
    }
  }
  size_t dim() const { return m_dim ; } 
};

struct WeightedDecisionTree {

  AmpGen::BinDT       m_binning; 
  std::vector<double> m_leafWeights;
  WeightedDecisionTree() {};
  WeightedDecisionTree( AmpGen::BinDT binning ) : m_binning(binning) {} 
  void generateLeafWeights( const int& dim, 
      const std::vector<double*>& sourcePointers, 
      const std::vector<double*>& targetPointers ){
    m_binning.makeNodes( sourcePointers, targetPointers );
    std::vector<double> weight_source(m_binning.const_nodes().size(),0);
    std::vector<double> weight_target(m_binning.const_nodes().size(),0);
    for( auto& evt : targetPointers ) weight_target[m_binning.getBinNumber(evt)] += evt[dim];
    for( auto& evt : sourcePointers ) weight_source[m_binning.getBinNumber(evt)] += evt[dim];
    m_leafWeights.resize( m_binning.const_nodes().size() );
    double weights = 0 ; 
    for( unsigned int i = 0 ; i < m_leafWeights.size(); ++i ){
      double leafWeight = weight_target[i] / weight_source[i];
      m_leafWeights[i] = leafWeight;
      weights         += leafWeight;
      DEBUG("Leaf["<<i<<"] weight = " << m_leafWeights[i] << " " << weight_target[i] << " / " << weight_source[i] ); 
    }
  }
  double getEventWeight( const double* event ) const {
    return m_leafWeights[ m_binning.getBinNumber( event ) ] ; 
  }
  void write( std::ofstream& output ){
    INFO("Writing: " << m_leafWeights.size() << " nodes");
    for( auto& w : m_leafWeights ) output << w << " ";
    output << "\n";
    m_binning.serialize(output);
  } 
  WeightedDecisionTree( std::ifstream& input ){
    std::string weights; 
    getline( input, weights );
    auto tokens = AmpGen::split( weights, ' ');
    for( auto& weight : tokens ) m_leafWeights.push_back( stod( weight ) );
    m_binning.readFromStream(input);
  }
};


class BoostedDecisionTree {
  private: 
    std::vector< WeightedDecisionTree > m_trees;
    unsigned int                        m_nTrees; 
    double                              m_learningRate; 
    double                              m_globalWeight; 
  public:
    BoostedDecisionTree( 
        Store source, 
        Store target, 
        const unsigned int& nTrees=50, 
        const unsigned int& maxDepth=10, 
        const double& learningRate=0.1,
        const size_t& minEvents=100) :
      m_nTrees(nTrees),
      m_learningRate( learningRate )
  {
    double w_total_target = 0 ; 
    for( unsigned int i = 0; i < target.n(); ++i) w_total_target += target.weight(i);
    double w_total_source = 0;

    m_globalWeight = double( target.n() ) / double( source.n() );
    INFO("Global weight = " << m_globalWeight );
    for( unsigned int i = 0; i < source.n(); ++i){
      source.setWeight(i, m_globalWeight);
      w_total_source += source.weight(i);
    }
    std::vector< unsigned > things( source.dim() );
    std::iota( things.begin(), things.end(), 0 );
    std::random_device rd;
    std::mt19937 g(rd());

    for( unsigned int tree = 0 ; tree < nTrees; ++tree ){ 
      AmpGen::ProfileClock pc; 
      std::shuffle( things.begin(), things.end(), g );
      m_trees.emplace_back( AmpGen::BinDT( AmpGen::MaxDepth(maxDepth),
            AmpGen::Dim(source.dim()),
            AmpGen::MinEvents(minEvents) ) );
      m_trees[tree].m_binning.setQueueOrdering( things );
      m_trees[tree].generateLeafWeights(source.dim(), source.make_pointers(), target.make_pointers() );
      double w_total = 0 ;
      for( unsigned int i = 0 ; i < source.n(); ++i ){
        double current_weight = source.weight(i);
        double new_weight = m_trees[tree].getEventWeight(source.evt(i));
        double updated_weight = current_weight * ( 1 + m_learningRate * ( new_weight - 1 ) );
        source.setWeight(i,updated_weight );
        w_total += updated_weight ; 
      }
      pc.stop();
      INFO("Generating DT  : " << tree << " / " << nTrees << " source weight = " << w_total << " weight[0] = " << source.weight(0) << " t = " << pc );
    }
  }
    double weight( const double* event ) const {
      double w = m_globalWeight;
      for( auto& t : m_trees ) w *= ( 1 + m_learningRate * ( t.getEventWeight( event ) - 1 ) ); 
      return w;
    }
    void writeToFile( const std::string& fname ){
      std::ofstream output(fname);
      output << "BDT " << m_nTrees << " " << m_globalWeight << " " << m_learningRate << "\n";
      for( auto& tree : m_trees ){ tree.write( output ); output << "\n"; }
      output.close();
    }

    BoostedDecisionTree( const std::string& filename ){
      std::ifstream input(filename);
      /// get header //
      std::string header;
      getline( input, header);
      auto tokens = AmpGen::split( header, ' ' );
      m_nTrees       = stod( tokens[1]  );
      m_globalWeight = stod( tokens[2]  );
      m_learningRate = stod( tokens[3]  );
      for( unsigned int i = 0 ; i < m_nTrees ; ++i ) m_trees.emplace_back( input );
      input.close();
    }
};

