#include <Rcpp.h>
using namespace Rcpp;

#include <string>
#include <iostream>
using namespace std;

//' @export PCRSim
class PCRSim {
private:
  // Number of loci. The loci of alleles being used in the system
  int m_nNumLoci;
  // Number of cells to be analysed in this case
  int m_nNumCells;
  // Number of PCR cycles. By default we use 28 cycles.
  int m_nNumPCRCycles;
   // Efficiency of the extraction process. Constrained to be between zero and one.
  double m_dExtractionEfficiency;
  // Efficiency of the PCR process. Constrained to be between zero and one.
  double m_dPCREfficiency;
  // Proportion of an aliquot that contains m_nNumCells that is taken forward to PCR.
  // Constrained to be between zero and one.
  double m_dAliquot;
  // Is this a haploid case? false by default.
  bool m_bHaploid;
  // Are we employing a stutter model? false by default.
  bool m_bStutter;
  double m_dAlpha;
  double m_dBeta;
public:
  // Number of simulations to run
  int m_nNumSims;

private:
  NumericVector rbinom(const int n, const double size, const double p){
    NumericVector r;

    for(int i = 0; i < n; i++){
      r.push_back(::Rf_rbinom(size, p));
    }

    return r;
  }

  NumericVector rbinom(const NumericVector& size, const double p){
    NumericVector r;

    for(NumericVector::const_iterator i = size.begin(); i != size.end(); i++){
      r.push_back(::Rf_rbinom(*i, p));
    }

    return r;
  }

  NumericVector rbinom(const NumericVector& size, const NumericVector& p){
    NumericVector r;
    NumericVector::const_iterator n = size.begin();
    NumericVector::const_iterator pn = p.begin();

    while(n != size.end() && pn != p.end()){
      r.push_back(::Rf_rbinom(*n, *pn));
      n++;
      pn++;
    }

    return r;
  }

  // // Calculates the heterozygous balance and the proportion of dropout events
  // double CalculateHb(NumericVector& alleleCount, NumericVector& Hb){
  //   int nLoc;
  //   int nA, nB;
  //   double dA,dB;
  //   int nDropout = 0;
  //
  //   for(nLoc=0;nLoc<m_nNumLoci;nLoc++)
  //   {
  //     nA = alleleCount[2 * nLoc];
  //     nB = alleleCount[2 * nLoc + 1];
  //
  //     if(nA==0||nB==0) // dropout
  //     {
  //       Hb[nLoc]=0;
  //       nDropout++;
  //     }
  //     else
  //     {
  //       dA = (double)nA;
  //       dB = (double)nB;
  //
  //       if(dA>dB)
  //         Hb[nLoc] = dB/dA;
  //       else
  //         Hb[nLoc] = dA/dB;
  //     }
  //   }
  //
  //   return (double)nDropout/(double)(2*m_nNumLoci);
  // }
  // Simulates the extraction processes.
  NumericVector Extract(void){
    NumericVector survived = rbinom(2 * m_nNumLoci, m_nNumCells, m_dExtractionEfficiency);
    return rbinom(survived, m_dAliquot);
  }

  void PCR(NumericVector& alleleCount, NumericVector& stutterAlleleCount){
    stutterAlleleCount = NumericVector::create(2 * m_nNumLoci, 0);

    for(int nCycle = 0; nCycle < m_nNumPCRCycles; nCycle++){
      alleleCount += rbinom(alleleCount, m_dPCREfficiency);

      if(m_bStutter){ // only do this if we're using a stutter model
        /* by PCRing stutter product first ensures that PCR only happens on product
        previous PCR cycle */
        stutterAlleleCount += rbinom(stutterAlleleCount, m_dPCREfficiency);

        NumericVector probStutter = rbeta(2 * m_nNumLoci, m_dAlpha, m_dBeta);
        stutterAlleleCount += rbinom(alleleCount, probStutter);
      }
    }
  }

  void PCR(NumericVector& alleleCount){
    for(int nCycle = 0; nCycle < m_nNumPCRCycles; nCycle++){
      alleleCount += rbinom(alleleCount, m_dPCREfficiency);
    }
  }

  void printProf(const NumericVector& alleleCount){
    ostringstream oss;

    for(int n = 0; n < m_nNumLoci; n++){
      oss << n + 1 << ":" << alleleCount[2 * n] << " " << alleleCount[2 * n + 1] << endl;
    }

    Rprintf("%s\n", oss.str().c_str());
  }


public:
  // Default constructor
  /* Default constructor.*/
  PCRSim(void){
    m_nNumCells = 0;
    m_nNumLoci = 0;
    m_nNumPCRCycles = 28;
    m_nNumSims = 0;
    m_dExtractionEfficiency = 0;
    m_dPCREfficiency = 0;
    m_dAliquot = 0;
    m_bHaploid = false;
    m_bStutter = false;
  }

  //  Properties
  double getAliquotPpn(void){
    return m_dAliquot;
  }

  void setAliquotPpn(double dAliquot){
    m_dAliquot = dAliquot;
  }

  NumericVector getAlpha(void){
    return NumericVector::create(m_dAlpha, m_dBeta);
  }

  void setAlpha(const NumericVector& alpha){
    m_dAlpha = alpha[0];
    m_dBeta = alpha[1];
  }

  double getExtractionEfficiency(void){
    return m_dExtractionEfficiency;
  }

  void setExtractionEfficiency(double dExtractionEfficiency){
    m_dExtractionEfficiency = dExtractionEfficiency;
  }

  bool getHaploid(void){
    return m_bHaploid;
  }

  void setHaploid(bool bHaploid){
    m_bHaploid = bHaploid;
  }

  int getNumCells(void){
    return m_nNumCells;
  }

  void setNumCells(int numCells){
    m_nNumCells = numCells;
  }

  int getNumLoci(void){
    return m_nNumLoci;
  }

  void setNumLoci(int numLoci){
    m_nNumLoci = numLoci;
  }

  int getNumPCRCycles(void){
    return m_nNumPCRCycles;
  }

  void setNumPCRCycles(int nNumPCRCycles){
    m_nNumPCRCycles = nNumPCRCycles;
  }

  double getPCREfficiency(void){
    return m_dPCREfficiency;
  }

  void setPCREfficiency(double dPCREfficiency){
    m_dPCREfficiency = dPCREfficiency;
  }

  bool getStutter(void){
    return m_bStutter;
  }

  void setStutter(bool bStutter){
    m_bStutter = bStutter;
  }

  void print(void){
    ostringstream Rout;
    Rout << "Simulation parameters" << endl;
    Rout << "----------------------------" << endl;
    Rout << "Number of cells: " << m_nNumCells << endl;
    Rout << "Extraction efficiency: " << m_dExtractionEfficiency << endl;
    Rout << "Aliquot: " << m_dAliquot << endl;
    Rout << "Number of loci: " << m_nNumLoci << endl;
    Rout << "Number of PCR cycles: " << m_nNumPCRCycles << endl;
    Rout << "PCR Efficiency: " << m_dPCREfficiency << endl;
    Rout << "Stutter: " << ( m_bStutter ? "true" : "false") << endl;

    if(m_bStutter){
      double dMean = m_dAlpha/(m_dAlpha+m_dBeta);
      double dVar = m_dAlpha*m_dBeta/((m_dAlpha+m_dBeta)*(m_dAlpha+m_dBeta)*(m_dAlpha+m_dBeta+1));
      Rout << "Stutter mean: " << dMean << endl;
      Rout << "Stutter variance: " << dVar << endl;
    }
    Rout << "Haploid: " << ( m_bHaploid ? "true" : "false") << endl;
    Rout << "-----" << endl;

    Rprintf("%s", Rout.str().c_str());
  }

  // This is where the simulation happens
  List simulate(void){
    NumericVector alleleCount(2 * m_nNumLoci, 0.0);
    NumericVector stutterAlleleCount(2 * m_nNumLoci, 0.0); //currently assuming back stutter only one repeat
    NumericVector Hb(m_nNumLoci, 0.0);

    List results;
    NumericMatrix Alleles(m_nNumSims, 2 * m_nNumLoci);
    NumericMatrix Stutters(m_nNumSims, 2 * m_nNumLoci);

    for(int nSim = 0; nSim < m_nNumSims; nSim++){
      alleleCount = Extract();

      if(m_bStutter){
        PCR(alleleCount, stutterAlleleCount);
        Alleles(nSim,_) = alleleCount;
        Stutters(nSim,_) = stutterAlleleCount;
      }else{
        PCR(alleleCount);
        Alleles(nSim,_) = alleleCount;
      }
      if(nSim % 100 == 0)
        cout << nSim << endl;
    }

    results["alleles"] = Alleles;

    if(m_bStutter){
      results["stutters"] = Stutters;
    }


    return results;
  }

  // List SimulatePCROnly(int n0, const string& strResultsFileName){
  //   NumericVector alleleCount(2 * m_nNumLoci, 0);
  //   List results;
  //
  //   for(int nSim = 0; nSim < m_nNumSims; nSim++){
  //     PCR(alleleCount);
  //    }
  //
  //   return results;
  // }
};


// Expose the classes
RCPP_MODULE(PCRSim) {
  using namespace Rcpp;

  class_<PCRSim>( "PCRSim")
    .default_constructor("Default constructor") // This exposes the default constructor

    .field( "numSims", &PCRSim::m_nNumSims) // number of simulations

    .method("getAlpha", &PCRSim::getAlpha)
    .method("setAlpha", &PCRSim::setAlpha)
    .method("print", &PCRSim::print) // this exposes the print method
    .method("simulate", &PCRSim::simulate) // This exposes the simulate method

    .property("aliqot", &PCRSim::getAliquotPpn, &PCRSim::setAliquotPpn)
    .property("extract", &PCRSim::getExtractionEfficiency, &PCRSim::setExtractionEfficiency)
    .property("haploid", &PCRSim::getHaploid, &PCRSim::setHaploid)
    .property("numCells", &PCRSim::getNumCells, &PCRSim::setNumCells)
    .property("numLoci", &PCRSim::getNumLoci, &PCRSim::setNumLoci)
    .property("numCycles", &PCRSim::getNumPCRCycles, &PCRSim::setNumPCRCycles)
    .property("PCReff", &PCRSim::getPCREfficiency, &PCRSim::setPCREfficiency)
    .property("stutter", &PCRSim::getStutter, &PCRSim::setStutter)
  ;
}


