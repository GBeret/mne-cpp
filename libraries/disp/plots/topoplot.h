//=============================================================================================================
/**
* @file     topoplot.cpp
* @author   Martin Henfling <martin.henfling@tu-ilmenau.de>;
*           Daniel Knobl <daniel.knobl@tu-ilmenau.de>;
* @version  1.0
* @date     September, 2015
*
* @section  LICENSE
*
* Copyright (C) 2014, Martin Henfling and Daniel Knobl. All rights reserved.
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
* @brief    Declaration of topo plot class.
*/

#ifndef TOPOPLOT_H
#define TOPOPLOT_H

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"
#include <disp/plots/helpers/colormap.h>

//*************************************************************************************************************
//=============================================================================================================
// Qt INCLUDES
//=============================================================================================================

#include <QImage>
#include <QGridLayout>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include <QGraphicsSceneResizeEvent>
#include <QtConcurrent>
//#include <QList>
//#include <QSize>

//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Splines>

namespace DISPLIB
{

struct TopoPlotInputData
{
    Eigen::MatrixXd signalMatrix;
    quint32 iRangeLow;
    quint32 iRangeHigh;
    qint32 dampingFactor;
    ColorMaps colorMap;
    QMap<QString, QPoint> topoMap;
    QSize topoMatrixSize;
    QSize imageSize;
    QList<Eigen::MatrixXd> topoMatrixList;
};

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

class DISPSHARED_EXPORT TopoPlot : public QThread
{

    Q_OBJECT

    typedef QMap<QString, QPointF> channelMap;
    typedef Eigen::MatrixXd MatrixXd;
    typedef QList<QImage> imageList;
    typedef ColorMaps colorMaps;

public:
    TopoPlot();
    //~TopoPlot();
    static QList<Eigen::MatrixXd> createTopoPlotMatrix(const TopoPlotInputData& inputData);
    static QList<QImage> createTopoPlotImages(const TopoPlotInputData& inputData);
    static void reduceMatrix(QList<Eigen::MatrixXd> &resultData, const QList<Eigen::MatrixXd> &data);
    static void reduceImages(QList<QImage> &resultData, const QList<QImage> &data);
    QList<QImage> createTopoPlotImageList(const MatrixXd signalMatrix, const channelMap layoutMap, const QSize topoMatrixSize, const QSize imageSize, const colorMaps cmap, const qint32 dampingFactor);

private:

    static QMap<QString,QPoint> createMapGrid(const QMap<QString,QPointF> layoutMap, const QSize topo_matrix_size);
    static Eigen::MatrixXd normSignal(Eigen::MatrixXd signalMatrix);
    static Eigen::MatrixXd createGridPointMatrix(const Eigen::MatrixXd normSignal, const QMap<QString, QPoint> mapGrid, const QSize gridPointMatrixSize, const qint32 timeSample);
    static QImage * creatPlotImage(const MatrixXd topoMatrix, const QSize imageSize, const ColorMaps cmap);
    static Eigen::MatrixXd calcNearestNeighboursInterpolation(Eigen::MatrixXd topoMatrix, const QMap<QString,QPoint> mapGrid);
    static Eigen::MatrixXd calcBilinearInterpolation(const MatrixXd topoMatrix, const QMap<QString, QPoint> mapGrid, const qint32 dampingFactor);

public slots:
    void recieveInputStartCalculation(const MatrixXd signalMatrix, const channelMap layoutMap, const QSize topoMatrixSize, const QSize imageSize, const colorMaps cmap, const qint32 dampingFactor);

signals:
    void sendResult(imageList topoPlotImageList, bool finished);
};

}

#endif // TOPOPLOT_H
