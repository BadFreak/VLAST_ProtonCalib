#ifndef MultiStreamAnalysis_H_
#define MultiStreamAnalysis_H_

#include "SniperKernel/AlgBase.h"

class MultiStreamAnalysis : public AlgBase
{
   public :

	MultiStreamAnalysis(const std::string& name);
	virtual ~MultiStreamAnalysis();

	virtual bool initialize();
	virtual bool execute();
	virtual bool finalize();

    private :
	int   m_iEvt;
};

#endif
