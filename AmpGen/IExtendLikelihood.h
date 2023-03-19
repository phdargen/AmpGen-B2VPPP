#ifndef AMPGEN_IEXTENDLIKELIHOOD_H
#define AMPGEN_IEXTENDLIKELIHOOD_H

#include <string>
#include <vector>

namespace AmpGen
{
  class MinuitParameterSet;
  class MinuitParameter;
  //class CoherentSum;
  class PolarisedSum;

  class IExtendLikelihood
  {
  public:
    virtual ~IExtendLikelihood() = default;
    virtual double getVal()  = 0;
    virtual void configure( const std::string& configString, AmpGen::PolarisedSum& pdf,
                            const MinuitParameterSet& mps ) = 0;
    virtual IExtendLikelihood* create()                     = 0;
  };

  class GaussianConstraint : public IExtendLikelihood
  {
  public:
    double getVal()  override;
    void configure( const std::string& configString, AmpGen::PolarisedSum& pdf,
                    const MinuitParameterSet& mps ) override;
    void configure( const std::string& configString, const MinuitParameterSet& mps );

    IExtendLikelihood* create() override { return new GaussianConstraint(); }
    static std::string _id;

  private:
    MinuitParameter* m_param;
    double m_mean;
    double m_sigma;
  };

  class PartialWidthConstraint : public IExtendLikelihood
  {
  public:
    double getVal()  override;
    void configure( const std::string& configString, AmpGen::PolarisedSum& pdf,
                    const AmpGen::MinuitParameterSet& mps ) override;
    IExtendLikelihood* create() override { return new PartialWidthConstraint(); }
    static std::string _id;

  private:
    PolarisedSum* m_pdf;
    double m_ratio;
    double m_weight;
    std::vector<unsigned int> m_denComponents;
    std::vector<unsigned int> m_numComponents;
  };

  class LASSO : public IExtendLikelihood
  {
  public:
    double getVal()  override;
    void configure( const std::string& configString, AmpGen::PolarisedSum& pdf,
                    const MinuitParameterSet& mps ) override;
    IExtendLikelihood* create() override { return new LASSO(); }
    static std::string _id;

  private:
    double m_lambda;
    PolarisedSum* m_pdf;
  };
    
  class ExtendedNormalisation : public IExtendLikelihood
    {
    public:
        double getVal()  override;
        void configure( const std::string& configString, AmpGen::PolarisedSum& pdf,
                       const MinuitParameterSet& mps ) override;
        IExtendLikelihood* create() override { return new ExtendedNormalisation(); }
        static std::string _id;
        
    private:
      PolarisedSum* m_pdf;
    };

    
} // namespace AmpGen

#endif
