//=============================================================================================================
/**
* @file     tfplotsceneitem.cpp
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
* @brief    Contains the implementation of the TFPlotSceneItem class.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "TFPlotSceneItem.h"

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

TFPlotSceneItem::TFPlotSceneItem(QString channelName, int channelNumber, QPointF channelPosition, QImage *tfImage, int channelKind, int channelUnit, QColor channelColor, bool bIsBadChannel)
: m_sChannelName(channelName)
, m_iChannelNumber(channelNumber)
, m_qpChannelPosition(channelPosition)
, m_tfImage(tfImage)
, m_cChannelColor(channelColor)
, m_bHighlightItem(false)
, m_iChannelKind(channelKind)
, m_iChannelUnit(channelUnit)
, m_bIsBadChannel(bIsBadChannel)
{
    this->setAcceptHoverEvents(true);
    this->setFlag(QGraphicsItem::ItemIsSelectable, true);
}

//*************************************************************************************************************

QRectF TFPlotSceneItem::boundingRect() const
{
    return QRectF(-25, -30, 50, 50);
}


//*************************************************************************************************************

void TFPlotSceneItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    Q_UNUSED(option);
    Q_UNUSED(widget);

    this->setPos(10*m_qpChannelPosition.x(), -10*m_qpChannelPosition.y());

    if(m_tfImage->isNull())
    {
        // Plot shadow
        painter->setPen(Qt::NoPen);
        painter->setBrush(Qt::darkGray);
        painter->drawEllipse(-12, -12, 30, 30);
        //Plot red if bad
        if(m_bIsBadChannel) {
            painter->setBrush(Qt::red);
            painter->drawEllipse(-15, -15, 30, 30);
        }
        else {
            painter->setBrush(m_cChannelColor);
            painter->drawEllipse(-15, -15, 30, 30);
        }

        //Plot selected item
        if(this->isSelected()){
            //painter->setPen(QPen(QColor(255,84,22), 5));
            painter->setPen(QPen(Qt::red, 5));
            painter->drawEllipse(-15, -15, 30, 30);
        }

        // Plot electrode name
        painter->setPen(QPen(Qt::black, 1));
        QStaticText staticElectrodeName = QStaticText(m_sChannelName);
        QSizeF sizeText = staticElectrodeName.size();
        painter->drawStaticText(-15+((30-sizeText.width())/2), -32, staticElectrodeName);
    }
    else
    {
        painter->drawImage(QRectF(-30, -15, 60, 30), *m_tfImage);
        // Plot electrode name
        painter->setPen(QPen(Qt::black, 1));
        QStaticText staticElectrodeName = QStaticText(m_sChannelName);
        QSizeF sizeText = staticElectrodeName.size();
        painter->drawStaticText(-30+((60-sizeText.width())/2), -32, staticElectrodeName);
    }
}
