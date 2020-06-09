//=============================================================================================================
/**
* @file     topoplot.cpp
* @author   Martin Henfling <martin.henfling@tu-ilmenau.de>;
*           Daniel Knobl <daniel.knobl@tu-ilmenau.de>;
* @version  1.0
* @date     March, 2019
*
* @section  LICENSE
*
* Copyright (C) 2019, Martin Henfling and Daniel Knobl. All rights reserved.
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

//#include <QtConcurrent>
#include <QtConcurrent/QtConcurrent>

//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================

#include <QImage>
#include <Eigen/Core>

namespace DISPLIB
{

struct TopoPlotInputData
{
    typedef ColorMaps colorMaps;

    Eigen::MatrixXd signalMatrix;
    quint32 iRangeLow;
    quint32 iRangeHigh;
    qint32 dampingFactor;
    qint32 colorMap;
    QMap<QString, QPoint> topoMap;
    QSize topoMatrixSize;
    QSize imageSize;
    QList<Eigen::MatrixXd> topoMatrixList;
};

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace Eigen;

class DISPSHARED_EXPORT TopoPlot : public QThread
{
    Q_OBJECT

    typedef QMap<QString, QPointF> channelMap;
    typedef Eigen::MatrixXd MatrixXd;
    typedef QList<QImage> imageList;
    typedef ColorMaps colorMaps;

public:
    TopoPlot();
    ~TopoPlot();

    //==========================================================================================================
    /**
    * TopoPlot_createTopoPlotImageList
    *
    * ### Topoplot public methode ###
    *
    * calculated foreach timesample and selected channel the topoplot images
    *
    * @param[in]    signalMatrix        signal with all channels
    * @param[in]    layoutMap           layout of the electrodes
    * @param[in]    topoMatrixSize      size of the new topomatrix
    * @param[in]    imageSize           size of the topoplot image
    * @param[in]    cmap                Colormap of the topoplot image
    * @param[in]    dampingFactor       damping factor of the interpolation
    *
    * @return List of topoplot images
    */
    QList<QImage> createTopoPlotImageList(const MatrixXd signalMatrix, const channelMap layoutMap, const QSize topoMatrixSize, const QSize imageSize, const colorMaps cmap, const qint32 dampingFactor);

protected:

    //==========================================================================================================
    /**
    * TopoPlot_createTopoPlotMatrix
    *
    * ### Topoplot protected methode ###
    *
    * calculated the interpolation for topoplot in matrix
    *
    * @param[in]    inputData        all datas fpr interpolation
    *
    *
    * @return List of topoplot matrix
    */
    static QList<Eigen::MatrixXd> createTopoPlotMatrix(const TopoPlotInputData& inputData);

    //==========================================================================================================
    /**
    * TopoPlot_createTopoPlotImages
    *
    * ### Topoplot protected methode ###
    *
    * tranfer topoplot matrix to image with colormap
    *
    * @param[in]    inputData        all datas topoplot interpolation
    *
    * @return List of topoplot matrix
    */
    static QList<QImage> createTopoPlotImages(const TopoPlotInputData& inputData);

    //=========================================================================================================
    /**
    * Sums up (reduces) the in parallel processed topoplot matrix.
    *
    * ### Topoplot protected methode ###
    *
    * @param[out] resultData    The result data.
    * @param[in]  data          The incoming, temporary result data.
    *
    * @return void
    */
    static void reduceMatrix(QList<Eigen::MatrixXd> &resultData, const QList<Eigen::MatrixXd> &data);

    //=========================================================================================================
    /**
    * Sums up (reduces) the in parallel processed topoplot images.
    *
    * ### Topoplot protected methode ###
    *
    * @param[out] resultData    The result data.
    * @param[in]  data          The incoming, temporary result data.
    *
    * @return void
    */
    static void reduceImages(QList<QImage> &resultData, const QList<QImage> &data);

    //==========================================================================================================
    /**
    * TopoPlot_createMapGrid
    *
    * ### Topoplot protected methode ###
    *
    * tranfer the layout in a spezific layoutsize
    *
    * @param[in]    layoutMap            layout of the electrodes
    * @param[in]    topo_matrix_size     size of the new matrix
    *
    * @return map of gridpoints
    */
    static QMap<QString,QPoint> createMapGrid(const QMap<QString,QPointF> layoutMap, const QSize topo_matrix_size);

    //==========================================================================================================
    /**
    * TopoPlot_createMapGrid
    *
    * ### Topoplot protected methode ###
    *
    * calculated the the value of each channel to the gridpointmap in a matrix
    *
    * @param[in]    signal                  signal with all channels
    * @param[in]    mapGrid                 map of gridpoints
    * @param[in]    gridPointMatrixSize     size gridpoint matrix
    * @param[in]    timeSample              the timesample of signal wich should calculat
    *
    * @return grid point matrix
    */
    static Eigen::MatrixXd createGridPointMatrix(const Eigen::MatrixXd signal, const QMap<QString, QPoint> mapGrid, const QSize gridPointMatrixSize, const qint32 timeSample);

    //==========================================================================================================
    /**
    * TopoPlot_creatPlotImage
    *
    * ### Topoplot protected methode ###
    *
    * tranfer topoplot matrix to image with colormap
    *
    * @param[in]    topoMatrix      the interpolated matrix
    * @param[in]    imageSize       size of the new image
    * @param[in]    cmap            colormap of the new image
    *
    * @return topoplot image
    */
    static QImage * creatPlotImage(const MatrixXd topoMatrix, const QSize imageSize, const qint32 cmap);

    //==========================================================================================================
    /**
    * TopoPlot_calcNearestNeighboursInterpolation
    *
    * ### Topoplot protected methode ###
    *
    * interpolat by  NearestNeighbours interpolation the grip point matrix to topoplot matrix
    *
    * @param[in]    topoMatrix      the interpolated matrix
    * @param[in]    mapGrid         map of gridpoints
    *
    * @return topoplot matrix
    */
    static Eigen::MatrixXd calcNearestNeighboursInterpolation(Eigen::MatrixXd topoMatrix, const QMap<QString,QPoint> mapGrid);

    //==========================================================================================================
    /**
    * TopoPlot_calcBilinearInterpolation
    *
    * ### Topoplot protected methode ###
    *
    * interpolate bilinear interpolation the grip point matrix to topoplot matrix
    *
    * @param[in]    topoMatrix      the interpolated matrix
    * @param[in]    mapGrid         map of gridpoints
    * @param[in]    dampingFactor   damping factor of the interpolation
    *
    * @return topoplot matrix
    */
    static Eigen::MatrixXd calcBilinearInterpolation(const MatrixXd topoMatrix, const QMap<QString, QPoint> mapGrid, const qint32 dampingFactor);

public slots:
    //==========================================================================================================
    /**
    * TopoPlot_recieveInputStartCalculation
    *
    * ### Topoplot public slot ###
    *
    * recieve the input and start calculation for topoplot
    *
    * @param[in]    signalMatrix        signal with all channels
    * @param[in]    layoutMap           layout of the electrodes
    * @param[in]    topoMatrixSize      size of the new topomatrix
    * @param[in]    imageSize           size of the topoplot image
    * @param[in]    cmap                Colormap of the topoplot image
    * @param[in]    dampingFactor       damping factor of the interpolation
    *
    * @return void
    */
    void recieveInputStartCalculation(const MatrixXd signalMatrix, const channelMap layoutMap, const QSize topoMatrixSize, const QSize imageSize, const ColorMaps cmap, const qint32 dampingFactor);

signals:

    //==========================================================================================================
    /**
    * TopoPlot_sendResult
    *
    * ### Topoplot signal ###
    *
    * send the list of topoplot images to the call function
    *
    * @param[in]    topoPlotImageList  result of the calculation (list of topoplot images)
    * @param[in]    finished           true, if calculation ready
    *
    * @return void
    */
    void sendResult(imageList topoPlotImageList, bool finished);
};

}

#endif // TOPOPLOT_H
