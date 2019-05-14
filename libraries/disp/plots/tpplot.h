//=============================================================================================================
/**
* @file     tpplot.cpp
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

#ifndef TPPLOT_H
#define TPPLOT_H

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

struct InterpolationInputData {
    Eigen::MatrixXd maInputData;
    quint32 iRangeLow;
    quint32 iRangeHigh;
    QList<QPoint> coors;
    qint32 x;
    qint32 y;

    //qint32 window_size;
};

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace Eigen;

class DISPSHARED_EXPORT Tpplot
{
public:
    Tpplot();

    QMap<QString,QPoint> createMapGrid(QMap<QString,QPointF> layoutMap, QSize topo_matrix_size);
    MatrixXd normSignal(MatrixXd signalMatrix);
    MatrixXd createGridPointMatrix(MatrixXd normSignal, QMap<QString, QPoint> mapGrid, QSize gridPointMatrixSize, qint32 timeSample);
    QImage * creatPlotImage(MatrixXd topoMatrix, QSize imageSize, ColorMaps cmap, bool nomalization);
    MatrixXd calcNearestNeighboursInterpolation(MatrixXd topoMatrix, QMap<QString,QPoint> mapGrid);
    MatrixXd calcBilinearInterpolation(MatrixXd topoMatrix, QMap<QString, QPoint> mapGrid);

    //=========================================================================================================
    /**
    * Calculates the spectogram matrix for a given input data matrix.
    *
    * @param[in] data       The input data.
    *
    * @return               The spectogram matrix.
    */
    static Eigen::MatrixXd computeyxf(const InterpolationInputData& data);

    static Eigen::MatrixXd computexyf(const InterpolationInputData& data);

    //=========================================================================================================
    /**
    * Sums up (reduces) the in parallel processed spectogram matrix.
    *
    * @param[out] resultData    The result data.
    * @param[in]  data          The incoming, temporary result data.
    */
    static void reduce(Eigen::MatrixXd &resultData,
                       const Eigen::MatrixXd &data);
};

}

#endif // TPPLOT_H
