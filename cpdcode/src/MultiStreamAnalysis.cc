#include "MultiStreamAnalysis.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperPtr.h"
#include "PodioDataSvc/DataHandle.hh"
#include "EventDataModel/CaloSimCellCollection.h"

using namespace edm;

DECLARE_ALGORITHM(MultiStreamAnalysis);

MultiStreamAnalysis::MultiStreamAnalysis(const std::string& name)
: AlgBase(name)
{
    m_iEvt = 0;
}

MultiStreamAnalysis::~MultiStreamAnalysis()
{
}

bool MultiStreamAnalysis::initialize()
{
    LogInfo << " initialized successfully" << std::endl;
    return true;
}

bool MultiStreamAnalysis::execute()
{
    LogDebug << "Processing event " << m_iEvt << std::endl;
    ++m_iEvt;

    // Read in calohits in two streams (defined in two sub-tasks)
    auto CaloSimHits_1 = getROColl(CaloSimCellCollection, "File1:CaloHitCol");
    auto CaloSimHits_2 = getROColl(CaloSimCellCollection, "File2:CaloHitCol");

    // Do some analysis here
    if (CaloSimHits_1 and CaloSimHits_2) {
        LogDebug << "size of File1:CaloHitCol " << CaloSimHits_1->size() << std::endl;
        LogDebug << "size of File2:CaloHitCol " << CaloSimHits_2->size() << std::endl;
    } 

    return true;
}

bool MultiStreamAnalysis::finalize()
{
   LogInfo << " finalized successfully" << std::endl;
   return true;
}
