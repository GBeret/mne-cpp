//=============================================================================================================
/**
* @file     tfplotscene.h
* @author   Daniel Knobl <daniel.knobl@tu-ilmenau.de>;
*           Martin Henfling <martin.henfling@tu-ilmenau.de>;
* @version  1.0
* @date     april, 2016
*
* @section  LICENSE
*
* Copyright (C) 2014, Daniel Knobl and Martin Henfling. All rights reserved.
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
* @brief    Contains the declaration of the TFPlotScene class.
*
*/
#ifndef TFPLOTSCENE_H
#define TFPLOTSCENE_H

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../../disp_global.h"
#include "layoutscene.h"
#include "tfplotsceneitem.h"
#include <fiff/fiff.h>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QGraphicsScene>
#include <QWidget>
#include <QMutableListIterator>


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE DISPLIB
//=============================================================================================================

namespace DISPLIB
{

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================


//=============================================================================================================
/**
* TFPlotScene...
*
* @brief The TFPlotScene class provides a reimplemented QGraphicsScene for 2D layout plotting.
*/
class DISPSHARED_EXPORT TFPlotScene : public LayoutScene
{
    Q_OBJECT

public:

    typedef TFPlotSceneItem TFSceneItem;

    //=========================================================================================================
    /**
    * Constructs a TFPlotScene.
    */
    explicit TFPlotScene(QGraphicsView* view, QObject *parent = 0);

    //=========================================================================================================
    /**
    * Updates layout data.
    *
    * @param [in] layoutMap layout data map.
    * @param [in] bad channel list.
    */
    void repaintItems(QList<TFPlotItemStruct> tfPlotStructList, QStringList badChannels);

    void fitInView();



signals:

    void current_item_dbclicked(TFSceneItem * item);

protected:
     QPointF                         m_mousePressPosition;           /**< The current mouse press location.*/

    //=========================================================================================================
    /**
    * Reimplemented double mouse press event.
    */
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent);

    int         m_iChannelTypeMode;
};

} // NAMESPACE DISPLIB

#endif // TFPLOTSCENE_H
