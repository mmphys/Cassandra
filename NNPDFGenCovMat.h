/**
 * Generate covariance matrix and inverse
 */
void Experiment::GenCovMat()
{
  fCovMat.clear();
  fSqrtCov.clear();
  fCovMat.resize(fNData, fNData, 0);
  fSqrtCov.resize(fNData, fNData, 0);

  for (int i = 0; i < fNData; i++)
    // Diagonal case
    fCovMat(i, i) = fStat[i] * fStat[i]; // stat error
  for (int l = 0; l < fNSys; l++)
    {
      auto &compsys = fSys[0][l];
      if (compsys.name == "SKIP") {continue;}
      bool iscorrelated =
        (compsys.name != "UNCORR" && compsys.name != "THEORYUNCORR");
      for (int i = 0; i < fNData; i++)
	{
	  double diagsig = 0.0;
	  double diagsignor = 0.0;
	  auto &sys = fSys[i][l];
	  switch (compsys.type)
	    {
	    case ADD:
	      diagsig += sys.add * sys.add;
	      break; // additive systematics
	    case MULT:
	      diagsignor += sys.mult * sys.mult;
	      break; // multiplicative systematics
	    }
	  fCovMat(i,i) += diagsig + diagsignor * fT0Pred[i] * fT0Pred[i] * 1e-4;

	  // No need to loop over the nondiagonal parts
	  if (!iscorrelated) {continue;}
	  for (int j = 0; j < i; j++)
	    {
	      auto &othersys = fSys[j][l];
	      // Hopefully easy enough for the compiler to fuse this up
	      decltype(sys.add) res;
	      switch (compsys.type)
		{
		case ADD:
		  res = sys.add * othersys.add;
		  break; // additive systematics
		case MULT:
		  res = sys.mult * othersys.mult * fT0Pred[i] * fT0Pred[j] * 1e-4;
		  break; // multiplicative systematics
		default:
		  throw NNPDF::RuntimeException("Experiment::GenCovMat", "sys type not recognized");
		  break;
		}
	      fCovMat(i, j) += res;
	      fCovMat(j, i) += res;
	    }
	}
    }
  CholeskyDecomposition(fCovMat, fSqrtCov);
}
