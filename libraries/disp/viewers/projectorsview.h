//=============================================================================================================
/**
* @file     projectorsview.h
* @author   Lorenz Esch <lesch@mgh.harvard.edu>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     July, 2018
*
* @section  LICENSE
*
* Copyright (C) 2018, Lorenz Esch and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    Declaration of the ProjectorsView Class.
*
*/

#ifndef PROJECTORSVIEW_H
#define PROJECTORSVIEW_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>
#include <QMap>


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

class QCheckBox;

namespace FIFFLIB {
    class FiffInfo;
}


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE DISPLIB
//=============================================================================================================

namespace DISPLIB
{


//*************************************************************************************************************
//=============================================================================================================
// DISPLIB FORWARD DECLARATIONS
//=============================================================================================================


//=============================================================================================================
/**
* DECLARE CLASS ProjectorsView
*
* @brief The ProjectorsView class provides a view to select between different modalities
*/
class DISPSHARED_EXPORT ProjectorsView : public QWidget
{
    Q_OBJECT

public:    
    typedef QSharedPointer<ProjectorsView> SPtr;              /**< Shared pointer type for ProjectorsView. */
    typedef QSharedPointer<const ProjectorsView> ConstSPtr;   /**< Const shared pointer type for ProjectorsView. */

    //=========================================================================================================
    /**
    * Constructs a ProjectorsView which is a child of parent.
    *
    * @param [in] parent        parent of widget
    */
    ProjectorsView(QWidget *parent = 0,
                Qt::WindowFlags f = Qt::Widget);

    //=========================================================================================================
    /**
    * Update the selection.
    *
    * @param [in] state    The state (unused).
    */
    void init(const QSharedPointer<FIFFLIB::FiffInfo> pFiffInfo);

protected:
    //=========================================================================================================
    /**
    * Slot called when user enables/disables all projectors
    */
    void onEnableDisableAllProj(bool status);

    //=========================================================================================================
    /**
    * Slot called when the projector check state changes
    */
    void onCheckProjStatusChanged(bool state);

    QList<QCheckBox*>                                   m_qListProjCheckBox;            /**< List of projection CheckBox. */
    QCheckBox*                                          m_pEnableDisableProjectors;     /**< Holds the enable disable all check box. */

    QSharedPointer<FIFFLIB::FiffInfo>                   m_pFiffInfo;                    /**< Connected fiff info. */

signals:
    //=========================================================================================================
    /**
    * Emit this signal whenever the user changes the projections.
    */
    void projSelectionChanged();

};

} // NAMESPACE

#endif // PROJECTORSVIEW_H
