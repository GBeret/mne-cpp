#ifndef FWD_H
#define FWD_H

#include "fwdrt_global.h"
#include "../../MNE/mne/mne.h"


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE FWDRTLIB
//=============================================================================================================

namespace FWDRTLIB
{

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace MNELIB;


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================


class FWDRTSHARED_EXPORT Fwd
{
public:
    Fwd();

    static bool clusterRois(const MNEForwardSolution* p_fwdIn, MNEForwardSolution* p_fwdOut, qint32 p_iClusterSize);

};

} // NAMESPACE

#endif // FWD_H
