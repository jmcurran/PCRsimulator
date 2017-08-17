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
  // Number of simulations to run
  int m_nNumSims;
  // Efficiency of the extraction process. Constrained to be between zero and one.
  double m_dExtractionEfficiency;
  // Efficiency of the PCR process. Constrained to be between zero and one.
  double m_dPCREfficiency;
  // Proportion of an aliquot that contains m_nNumCells that is taken forward to PCR.
  // Constrained to be between zero and one.
  double m_dAliquot;
  // Probability of an allele stuttering on a given PCR cycle.
  double m_dStutter;
  // Is this a haploid case? false by default.
  bool m_bHaploid;
  // Are we employing a stutter model? false by default.
  bool m_bStutter;
  double m_dAlpha;
  double m_dBeta;

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
  // Standard Constructor
  /* Standard Constructor
  \param nNumCells - the number of cells to be analysed in this case
  \param nNumAlleles - the number of alleles in the system
  \param dPCREfficiency - the efficiency ofthe PCR process
  \param dAliquot - the proportion of an aliquot that contains m_nNumCells that is taken forward to PCR
  \param nNumPCRCycles - the number of PCRCycles
  \param bStutter - are we using a stutter model.
  \param bHaploid - is this a haploid case?.*/
  PCRSim(int nNumCells, int nNumLoci, double dExtractionEfficiency, double dPCREfficiency, double dAliquot,
         double dStutter, int nNumPCRCycles = 28, bool bStutter = false, bool bHaploid = false){
    m_nNumCells = nNumCells;
    m_nNumLoci = nNumLoci;
    m_dExtractionEfficiency = dExtractionEfficiency;
    m_dPCREfficiency = dPCREfficiency;
    m_dAliquot = dAliquot;
    m_dStutter = dStutter;
    m_nNumPCRCycles = nNumPCRCycles;
    m_bHaploid = bHaploid;
    m_bStutter = bStutter;
  }
  // This is where the simulation happens
  List simulate(){
    ostringstream Rout;
    Rout << "Simulation parameters" << endl;
    Rout << "Number of cells: " << m_nNumCells << endl;
    Rout << "Extraction efficiency: " << m_dExtractionEfficiency << endl;
    Rout << "Aliquot: " << m_dAliquot << endl;
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

    NumericVector alleleCount(2 * m_nNumLoci, 0.0);
    NumericVector stutterAlleleCount(2 * m_nNumLoci, 0.0); //currently assuming back stutter only one repeat
    NumericVector Hb(m_nNumLoci, 0.0);

    List results;
    results["alleles"] = NumericMatrix::create(m_nNumSims, 2 * m_nNumLoci);

    if(m_bStutter){
      results["stutters"] = NumericMatrix(m_nNumSims, 2 * m_nNumLoci);
    }

    for(int nSim = 0; nSim < m_nNumSims; nSim++){
      alleleCount = Extract();

      if(m_bStutter){
        PCR(alleleCount, stutterAlleleCount);
       as<NumericMatrix>(results["alleles"])(nSim,_) = alleleCount;
       as<NumericMatrix>(results["stutters"])(nSim,_) = stutterAlleleCount;
      }else{
        PCR(alleleCount);
      }
      if(nSim % 100 == 0)
        cout << nSim << endl;
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
    .constructor<PCRSim>("Constructor with an argument") // This exposes the other constructor
    .method("Simulate", &PCRSim::simulate) // This exposes the simulate method
  ;
}


