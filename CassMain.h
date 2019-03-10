#if !defined(__CassMain_h__)
#define __CassMain_h__

# include "Cassandra.h"

namespace Cassandra
{
  // This is an entire set of data - doesn't need to understand details

  class DataManager
  {
  public:
    // Unknown must be last
    enum NodeType : int {
      BCDMS, SLAC, JLAB, BoNuS, Unknown };

  protected:
    DataSet     * m_DataSets[Unknown];
    ModelParams & m_Params;

  public:
    DataManager(ModelParams & mp, int iNumReplicas);
    ~DataManager();
    DataSet & operator[]( NodeType nt ) { return * m_DataSets[nt]; }

  public:
    size_t size(void) const noexcept;
    bool MakeOnePrediction(DataSet & ds,
			const char * pPDFSet, const char * pMassScheme,
			const char * pszOutFilePrefix,
			int iSeq, double W2Min, double Q2Min, bool bStrictCut,
			size_t SampleSize, bool bShowCut,
			double xBinSize, double Q2BinSize);
    bool LoadData( const char * fileName );
    void ApplyCutoff( double W2, double Q2, bool bStrictCut,
		      size_t SampleSize, bool bShowCut );
    void ComputeOverallStats( size_t SampleSize );
    bool InitPrediction(const DataSet & ds, const char * pszOutFilePrefix);
    bool MakePrediction(const char * pPDFSet, const char * pMassScheme,
			const char * pszOutFilePrefix,
			bool bStrictFirstCut,
			int SampleSize, bool bDumpOnly, bool bShowCut,
			double W2Min, double W2Max,
			double Q2Min, double Q2Max, int NumSteps,
			double xBinSize, double Q2BinSize);
  };
}
#endif
