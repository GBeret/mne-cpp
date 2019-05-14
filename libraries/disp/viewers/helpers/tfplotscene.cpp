//=============================================================================================================
/**
* @file     tfplotscene.cpp
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
* @brief    Contains the implementation of the TFPlotScene class.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "tfplotscene.h"
#include <QMouseEvent>

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISPLIB;
using namespace std;

//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

TFPlotScene::TFPlotScene(QGraphicsView* view, QObject* parent)
: LayoutScene(view, parent)
, m_iChannelTypeMode(FIFFV_MEG_CH)
{
}

//*************************************************************************************************************

void TFPlotScene::repaintItems(QList<TFPlotItemStruct> tfPlotStructList, QStringList badChannels)
{
    this->clear();

    for(qint32 i = 0; i < tfPlotStructList.length(); i++)
    {
        TFPlotSceneItem* TFPlotSceneItemTemp;
        if(tfPlotStructList.at(i).channelName.contains("EEG"))
            TFPlotSceneItemTemp = new TFPlotSceneItem(tfPlotStructList.at(i).channelName,
                                                      0,
                                                      tfPlotStructList.at(i).coordinates,
                                                      tfPlotStructList.at(i).tfPlotImage,
                                                      FIFFV_EEG_CH,
                                                      FIFF_UNIT_T_M,
                                                      Qt::blue,
                                                      badChannels.contains(tfPlotStructList.at(i).channelName));
        else
            TFPlotSceneItemTemp = new TFPlotSceneItem(tfPlotStructList.at(i).channelName,
                                                      0,
                                                      tfPlotStructList.at(i).coordinates,
                                                      tfPlotStructList.at(i).tfPlotImage,
                                                      FIFFV_MEG_CH,
                                                      FIFF_UNIT_T_M,
                                                      Qt::blue,
                                                      badChannels.contains(tfPlotStructList.at(i).channelName));
        this->addItem(TFPlotSceneItemTemp);

    }
}

//*************************************************************************************************************

void TFPlotScene::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
    QGraphicsItem *item = this->itemAt(mouseEvent->buttonDownScenePos(Qt::LeftButton), QTransform());
    if(item)
    {
        TFPlotSceneItem *tfItem = dynamic_cast<TFPlotSceneItem *>(item);
        if(tfItem)
           emit current_item_dbclicked(tfItem);
    }
    else if(mouseEvent->button() == Qt::LeftButton)
        m_qvView->fitInView(this->itemsBoundingRect(), Qt::KeepAspectRatio);

    QGraphicsScene::mouseDoubleClickEvent(mouseEvent);
}

//*************************************************************************************************************

void TFPlotScene::fitInView()
{
     m_qvView->fitInView(this->itemsBoundingRect(), Qt::KeepAspectRatio);
}




























