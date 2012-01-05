#include "MiniTree/Selection/interface/SampleInfo.h"

ClassImp(SampleInfo)

SampleInfo::SampleInfo()
{
  sampleName="";
  mcEvtType = 0;
  pileup.clear();
  puWeights.clear();
}

SampleInfo::~SampleInfo()
{
}
