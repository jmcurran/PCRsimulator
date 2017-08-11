#include <Rcpp.h>
using namespace Rcpp;

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class PCRSim {
private:
  // Number of loci. The loci of alleles being used in the system
  int m_nNumLoci;
  // Number of cells to be analysed in this case
  int m_nNumCells;
  // Number of PCR cycles
  /* Number of PCR cycles. By default we use 28 cycles.*/
  int m_nNumPCRCycles;
  // Number of simulations to run
  /* Number of simulations to run*/
  int m_nNumSims;
  // Efficiency of the extraction process
  /* Efficiency of the extraction process. Constrained to be between zero and one */
  double m_dExtractionEfficiency;
  // Efficiency of the PCR process
  /* Efficiency of the PCR process. Constrained to be between zero and one.*/
  double m_dPCREfficiency;
  // Proportion of an aliquot that contains m_nNumCells that is taken forward to PCR
  /* Proportion of an aliquot that contains m_nNumCells that is taken forward to PCR.
  Constrained to be between zero and one.*/
  double m_dAliquot;
  // Probability of an allele stuttering on a given PCR cycle
  /* Probability of an allele stuttering on a given PCR cycle.*/
  double m_dStutter;
  // Is this a haploid case?
  // Is this a haploid case? \a false by default
  bool m_bHaploid;
  // Are we employing a stutter model?
  /* Are we employing a stutter model? \a false by default */
  bool m_bStutter;
  string m_strResFileName;
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

  // Calculates the heterozygous balance and the proportion of dropout events
  double CalculateHb(NumericVector& vAlleleCount, double *pHb){
    int nLoc;
    int nA, nB;
    double dA,dB;
    int nDropout = 0;

    for(nLoc=0;nLoc<m_nNumLoci;nLoc++)
    {
      nA = pAlleleCount[nLoc][0];
      nB = pAlleleCount[nLoc][1];

      if(nA==0||nB==0) // dropout
      {
        pHb[nLoc]=0;
        nDropout++;
      }
      else
      {
        dA = (double)nA;
        dB = (double)nB;

        if(dA>dB)
          pHb[nLoc] = dB/dA;
        else
          pHb[nLoc] = dA/dB;
      }
    }

    return (double)nDropout/(double)(2*m_nNumLoci);
  }
  // Simulates the extraction processes
  /* Simulations the extraction process.
  \param vAlleleCount */
  NumericVector Extract(void){
    NumericVector survived = rbinom(2 * m_nNumLoci, m_nNumCells, m_dExtractionEfficiency);
    return rbinom(survived, m_dAliquot);
  }

  void PCR(NumericVector& alleleCount, NumericVector& stutterAlleleCount){

    if(m_dPCREfficiency < 1){
      int nLoc;
      int nCycle;

      if(m_bStutter){
        stutterAlleleCount = NumericVector(2 * m_nNumLoci, 0);
      }

      for(nCycle = 0; nCycle < m_nNumPCRCycles; nCycle++)
        alleleCount += rbinom(alleleCount, m_dPCREfficiency);

        if(m_bStutter){ // only do this if we're using a stutter model
            /* by PCRing stutter product first ensures that PCR only happens on product
            previous PCR cycle */
            stutterAlleleCount += rbinom(stutterAlleleCount, m_dPCREfficiency);
            double dStutter = rbeta(m_dAlpha, m_dBeta);
            double dStutter = rbeta(m_dAlpha, m_dBeta);
            if(pAlleleCount[nLoc][1]>33554432L) // 2^25 somewhat arbitrary use a normal approx.
              pStutterAlleleCount[nLoc][1]+= rbinomapprox(pAlleleCount[nLoc][1], dStutter);
            else
              pStutterAlleleCount[nLoc][1]+= rbinom(pAlleleCount[nLoc][1], dStutter);
            if(f1!=NULL)
              *f1 << string64(pStutterAlleleCount[nLoc][1]) << ',' << dStutter << endl;
          }
        }
    }
  }

  void Write(ofstream &f1, NumericVector& vAlleleCount, double *pHb, double **pStutterAlleleCount = NULL);
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
  /* This is where the simulation happens.
  \param nSimulations - number of simulations to carry out*/
  bool Simulate(){
    ofstream f1(m_strResFileName.c_str());
    ofstream f2("pcrresults.csv");

    if(f1.is_open())
    {

      time_t s = time(NULL);

      srand(s);

      init_generator(rand(),rand());
      initnorm(rand());

      Rout << "Simulation parameters" << endl;
      Rout << "Number of cells: " << m_nNumCells << endl;
      Rout << "Extraction efficiency: " << m_dExtractionEfficiency << endl;
      Rout << "Aliquot: " << m_dAliquot << endl;
      Rout << "PCR Efficiency: " << m_dPCREfficiency << endl;
      Rout << "Stutter: " << ( m_bStutter ? "true" : "false") << endl;

      if(m_bStutter){
        double dMean = m_dAlpha/(m_dAlpha+m_dBeta);
        double dVar = m_dAlpha*m_dBeta/((m_dAlpha+m_dBeta)*(m_dAlpha+m_dBeta)*(m_dAlpha+m_dBeta+1));
        cout << "Stutter mean: " << dMean << endl;
        cout << "Stutter variance: " << dVar << endl;
      }
      Rout << "Haploid: " << ( m_bHaploid ? "true" : "false") << endl;
      Rout << "-----" << endl;

      int nLoc, nSim;
      NumericVector alleleCount(2 * m_nNumLoci);
      NumericVector stutterAlleleCount(2 * m_nNumLoci); // currently assuming back stutter only one repeat
      NumericVector Hb(m_nNumLoci);
      double dProbDropOut;


      for(nSim = 1; nSim <= m_nNumSims; nSim++){
        alleleCount = Extract();

        if(m_bStutter){
          PCR(alleleCount, stutterAlleleCount);
        }else{
          PCR(alleleCount);
        }

        //dProbDropOut = CalculateHb(pAlleleCount ,pHb);

        /* if(m_bStutter)
          Write(f1, pAlleleCount, pHb, pStutterAlleleCount);
        else
          Write(f1, pAlleleCount, pHb);

        if(nSim%100==0)
          cout << nSim << endl;
        */
      }
  }

  double SimulatePCROnly(int n0, const string& strResultsFileName){
    ofstream f1(strResultsFileName.c_str());

    if(f1.is_open()){
      int nLoc, nSim;
      NumericVector vAlleleCount(m_nNumLoci);
      double *pHb = new double[m_nNumLoci];
      double dProbDropOut;

      for(nLoc=0;nLoc<m_nNumLoci;nLoc++)
      {
        pAlleleCount[nLoc] = new double[2];
        pAlleleCount[nLoc][0] = n0;
        pAlleleCount[nLoc][1] = n0;
      }

      for(nSim=0;nSim<m_nNumSims;nSim++)
      {
        PCR(pAlleleCount);
        CalculateHb(pAlleleCount, pHb);
        Write(f1, pAlleleCount, pHb);
      }

      f1.close();

      delete [] pHb;

      for(nLoc=0;nLoc<m_nNumLoci;nLoc++)
        delete [] pAlleleCount[nLoc];
      delete [] pAlleleCount;

      return 0;
    }
    return 1;
  }
};


