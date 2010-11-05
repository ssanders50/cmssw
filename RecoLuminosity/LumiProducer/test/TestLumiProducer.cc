
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <cmath>
#include <iostream>


namespace edm {
  class EventSetup;
}

using namespace std;
using namespace edm;

namespace edmtest
{

  class TestLumiProducer : public edm::EDAnalyzer
  {
  public:

    explicit TestLumiProducer(edm::ParameterSet const&);
    virtual ~TestLumiProducer();

    virtual void analyze(edm::Event const& e, edm::EventSetup const& c);
    virtual void endLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& c);
  };

  // -----------------------------------------------------------------

  TestLumiProducer::TestLumiProducer(edm::ParameterSet const& ps)
  {
  }

  // -----------------------------------------------------------------

  TestLumiProducer::~TestLumiProducer()
  {
  }

  // -----------------------------------------------------------------

  void TestLumiProducer::analyze(edm::Event const& e,edm::EventSetup const&)
  {
    LuminosityBlock const& lumiBlock = e.getLuminosityBlock();

    Handle<LumiDetails> lumiDetails;
    lumiBlock.getByLabel("lumiProducer", lumiDetails);

    Handle<LumiSummary> lumiSummary;
    lumiBlock.getByLabel("lumiProducer", lumiSummary);
  }

  // -----------------------------------------------------------------

  void TestLumiProducer::endLuminosityBlock(LuminosityBlock const& lumiBlock, EventSetup const& c) {

    Handle<LumiSummary> lumiSummary;
    lumiBlock.getByLabel("lumiProducer", lumiSummary);
    if(lumiSummary->isValid()){
      std::cout << *lumiSummary << "\n";
    }else{
      std::cout << "no valid lumi summary data" <<std::endl;
    }
    Handle<LumiDetails> lumiDetails;
    lumiBlock.getByLabel("lumiProducer", lumiDetails);
    if(lumiDetails->isValid()){
      //std::cout << *lumiDetails << "\n";
      std::cout<<"lumivalue 5 "<< lumiDetails->lumiValue("OCC1",5)<<std::endl;
      std::cout<<"lumivalue 200 "<< lumiDetails->lumiValue("OCC1",200)*6.37<<std::endl;
      std::cout<<"lumivalue 880 "<< lumiDetails->lumiValue("OCC1",880)*6.37<<std::endl;
      std::cout<<"lumivalue 890 "<< lumiDetails->lumiValue("OCC1",890)*6.37<<std::endl;
      std::cout<<"lumivalue 900 "<< lumiDetails->lumiValue("OCC1",900)*6.37<<std::endl;
      std::cout<<"lumivalue 950 "<< lumiDetails->lumiValue("OCC1",950)*6.37<<std::endl;
      std::cout<<"lumivalue 1500 "<< lumiDetails->lumiValue("OCC1",1500)*6.37<<std::endl;
      std::cout<<"lumivalue 2000 "<< lumiDetails->lumiValue("OCC1",2000)*6.37<<std::endl;
      std::cout<<"lumivalue 2500 "<< lumiDetails->lumiValue("OCC1",2500)*6.37<<std::endl;
      std::cout<<"lumivalue 3500 "<< lumiDetails->lumiValue("OCC1",3500)*6.37<<std::endl;
    }else{
      std::cout << "no valid lumi detail data" <<std::endl;
    }
    // We know the content we put into the objects in the
    // configuration, manually check to see that we can
    // retrieve the same values.

    // A small value to allow for machine precision variations when
    // comparing small numbers.
    /*
    double epsilon = 0.001;

    for (int i = 0; i < 5; ++i) {

      if ( (fabs(lumiDetails->lumiEtSum(i)     - (100.0 + i)) > epsilon) ||
           (fabs(lumiDetails->lumiEtSumErr(i)  - (200.0 + i)) > epsilon) ||
           (lumiDetails->lumiEtSumQual(i) != (300 + i)) ||
           (fabs(lumiDetails->lumiOcc(i)       - (400.0 + i)) > epsilon) ||
           (fabs(lumiDetails->lumiOccErr(i)    - (500.0 + i)) > epsilon) ||
           (lumiDetails->lumiOccQual(i) != (600 + i)) ) {
	std::cerr << "TestLumiProducer: Values read from LumiDetails object do not match input values (1)\n";
        abort();
      }    
      if ( (fabs(lumiDetails->lumiEtSum()[i]     - (100.0 + i)) > epsilon) ||
           (fabs(lumiDetails->lumiEtSumErr()[i]  - (200.0 + i)) > epsilon) ||
           (lumiDetails->lumiEtSumQual()[i] != (300 + i)) ||
           (fabs(lumiDetails->lumiOcc()[i]       - (400.0 + i)) > epsilon) ||
           (fabs(lumiDetails->lumiOccErr()[i]    - (500.0 + i)) > epsilon) ||
           (lumiDetails->lumiOccQual()[i] != (600 + i)) ) {
	std::cerr << "TestLumiProducer: Values read from LumiDetails object do not match input values (2)\n";
        abort();
      }    

      if ( (lumiSummary->l1RateCounter(i)  != (10 + i)) ||
           (lumiSummary->l1Scaler(i)       != (20 + i)) ||
           (lumiSummary->hltRateCounter(i) != (30 + i)) ||
           (lumiSummary->hltScaler(i)      != (40 + i)) ||
           (lumiSummary->hltInput(i)       != (50 + i)) ) {
        std::cerr << "TestLumiProducer: Values read from LumiSummary object do not match input values (1)\n";
        abort();
      }

      if ( (lumiSummary->l1RateCounter()[i]  != (10 + i)) ||
           (lumiSummary->l1Scaler()[i]       != (20 + i)) ||
           (lumiSummary->hltRateCounter()[i] != (30 + i)) ||
           (lumiSummary->hltScaler()[i]      != (40 + i)) ||
           (lumiSummary->hltInput()[i]       != (50 + i)) ) {
        std::cerr << "TestLumiProducer: Values read from LumiSummary object do not match input values (2)\n";
        abort();
      }
    }

    if ( (fabs(lumiSummary->avgInsDelLumi()    - 1.0)  > epsilon) ||
         (fabs(lumiSummary->avgInsDelLumiErr() - 2.0)  > epsilon) ||
         (fabs(lumiSummary->avgInsRecLumi()    - 1.0 * 0.95)  > epsilon) ||
         (fabs(lumiSummary->avgInsRecLumiErr() - 2.0 * 0.95)  > epsilon) ||
         (fabs(lumiSummary->deadFrac()      - 0.05) > epsilon) ||
         (fabs(lumiSummary->liveFrac()      - 0.95) > epsilon) ||
         (lumiSummary->lumiSecQual() != 3) ||
         (lumiSummary->lsNumber() != 5) ) {
      std::cerr << "TestLumiProducer: Values read from LumiSummary object do not match input values (3)\n";
      abort();
    }
    */
  }
}

using edmtest::TestLumiProducer;

DEFINE_FWK_MODULE(TestLumiProducer);
