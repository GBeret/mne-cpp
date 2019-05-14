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
* @brief    Implementation of topo plot class.
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "tpplot.h"
#include "math.h"
#include <limits>
#include <iostream>
#include <QDebug>
#include <QElapsedTimer>
#include <QThread>
#include <QtConcurrent>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISPLIB;

Tpplot::Tpplot() {}

//*************************************************************************************************************

QMap<QString,QPoint> Tpplot::createMapGrid(QMap<QString,QPointF> layoutMap, QSize topo_matrix_size)
{    
    QMap<QString, QPoint> layoutMapGrid;
    qreal minXCoor = std::numeric_limits<int>::max();
    qreal maxXCoor = std::numeric_limits<int>::min();
    qreal minYCoor = std::numeric_limits<int>::max();
    qreal maxYCoor = std::numeric_limits<int>::min();
    qreal factorXCoor = 0;
    qreal factorYCoor = 0;
    QMapIterator<QString, QPointF> coor(layoutMap);

    // find min and max values
    while (coor.hasNext())
    {
        coor.next();
        if(coor.value().x() < minXCoor)
            minXCoor = coor.value().x();
        if(coor.value().x() > maxXCoor)
            maxXCoor = coor.value().x();
        if(coor.value().y() < minYCoor)
            minYCoor = coor.value().y();
        if(coor.value().y() > maxYCoor)
            maxYCoor = coor.value().y();
    }

    // 2px at the left edge
    minXCoor -= 2;
    minYCoor -= 2;
                                                                // 2px at the right edge
    factorXCoor = (maxXCoor - minXCoor) / (qreal)(topo_matrix_size.width() - 2);
    factorYCoor = (maxYCoor - minYCoor) / (qreal)(topo_matrix_size.height() - 2);

    // reset iterator
    coor.toFront();
    while (coor.hasNext())
    {
        coor.next();
        layoutMapGrid[coor.key()] = QPoint(qRound((coor.value().x() - minXCoor) / factorXCoor), qRound((coor.value().y() - minYCoor) / factorYCoor));

        // just safty :)
        if(layoutMapGrid[coor.key()].x() >= topo_matrix_size.width())
            layoutMapGrid[coor.key()].setX(topo_matrix_size.width() - 1);
        if(layoutMapGrid[coor.key()].y() >= topo_matrix_size.height())
            layoutMapGrid[coor.key()].setY(topo_matrix_size.height() - 1);
        if(layoutMapGrid[coor.key()].x() < 0)
            layoutMapGrid[coor.key()].setX(0);
        if(layoutMapGrid[coor.key()].y() < 0)
            layoutMapGrid[coor.key()].setY(0);
    }
    return  layoutMapGrid;
}

//*************************************************************************************************************

MatrixXd Tpplot::normSignal(MatrixXd signalMatrix)
{
    //normalisation foreach channel
    for(qint32 chnI = 0; chnI < signalMatrix.cols(); chnI++)   // over all Channels
    {
        VectorXd channel = signalMatrix.col(chnI);

        qreal max = channel.maxCoeff();
        qreal min = channel.minCoeff();
        for(qint32 i = 0; i < channel.rows(); i++)
            channel(i) -= min;
        channel /= ((max - min) / 1.0);
        signalMatrix.col(chnI) = channel;
    }
    return signalMatrix;
}

//*************************************************************************************************************

MatrixXd Tpplot::createGridPointMatrix(MatrixXd signal, QMap<QString,QPoint> mapGrid, QSize gridPointMatrixSize, qint32 timeSample)
{
    qint32 channel = 0;   
    MatrixXd tp_map =  MatrixXd::Zero(gridPointMatrixSize.height(), gridPointMatrixSize.width());

    // count of elektrods != count auf channels
    if(mapGrid.count() != signal.cols())
        return tp_map;

    QMapIterator<QString, QPoint> coor(mapGrid);
    while (coor.hasNext())
    {
        coor.next();
        qreal value = signal(timeSample, channel);
        tp_map(mapGrid[coor.key()].y(), mapGrid[coor.key()].x()) = value;
        channel++;
    }
    return tp_map;
}

//*************************************************************************************************************

QImage * Tpplot::creatPlotImage(MatrixXd topoMatrix, QSize imageSize, ColorMaps cmap, bool nomalization)
{
    //normalisation of the tp-matrix
    if(nomalization)
    {
        topoMatrix.array() += abs(topoMatrix.minCoeff());
        topoMatrix /= topoMatrix.maxCoeff();
    }

    qint32 y_factor =  topoMatrix.rows() / imageSize.height();
    qint32 x_factor =  topoMatrix.cols() / imageSize.width();
    if(y_factor == 0) y_factor = 1;
    if(x_factor == 0) x_factor = 1;

    //setup image
    QImage * topoImage = new QImage(imageSize, QImage::Format_RGB32);

    //setup pixelcolors in image
    QColor color;
    qint32 ximage = 0;
    qint32 yimage = 0;

    for ( qint32 y = 0; y < topoMatrix.rows(); y = y + y_factor)
    {
        ximage = 0;
        for ( qint32 x = 0; x < topoMatrix.cols(); x = x + x_factor )
        {
            switch  (cmap)
            {
                case Jet:
                    color.setRgb(ColorMap::valueToJet((topoMatrix(y, x))));
                    break;
                case Hot:
                    color.setRgb(ColorMap::valueToHot(topoMatrix(y, x)));
                    break;
                case HotNeg1:
                    color.setRgb(ColorMap::valueToHotNegative1(topoMatrix(y, x)));
                    break;
                case HotNeg2:
                    color.setRgb(ColorMap::valueToHotNegative2(topoMatrix(y, x)));
                    break;
                case Bone:
                    color.setRgb(ColorMap::valueToBone(topoMatrix(y, x)));
                    break;
                case RedBlue:
                    color.setRgb(ColorMap::valueToRedBlue(topoMatrix(y, x)));
                    break;
            }
            if(ximage < topoImage->width() && yimage < topoImage->height())
                topoImage->setPixel(ximage, topoImage->height() - 1 - yimage,  color.rgb());
            ximage++;
        }
        yimage++;
    }
    return topoImage;
}

//*************************************************************************************************************

MatrixXd Tpplot::calcNearestNeighboursInterpolation(MatrixXd topoMatrix, QMap<QString, QPoint> mapGrid)
{

    QList<QPoint> coors = mapGrid.values();

    //nearest
    for(qint32 y_axis = 0; y_axis < topoMatrix.rows(); y_axis++) //y_axis to interpolate among this axis
    {
        for(qint32 x_axis = 0; x_axis < topoMatrix.cols(); x_axis++) //x_axis to interpolate among this axis
        {

            //qint32 pointcounter = 0; //number of known points used to calc interpolation
            qint32 x = 0;
            qint32 y = 0;
            qint32 x_max = 0;
            qint32 y_max = 0;
            qreal temp_scalar = 0;

            for(qint32 i = 0; i < coors.length(); i++)
            {
                x = coors[i].x(); //known Coordinate x
                y = coors[i].y(); // known coordinate y

                if(!(x==x_axis && y==y_axis))
                {
                    qreal scalar = 1 - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) / sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) ;

                    if (scalar >= temp_scalar)
                    {
                        temp_scalar = scalar;
                        x_max = x;
                        y_max = y;
                        //result_yxf(y_axis, x_axis) += scalar*topoMatrix(y,x);
                        //qreal test = topoMatrix(y,x);
                        //pointcounter++;
                    }
                }
            }
            topoMatrix(y_axis, x_axis) = topoMatrix(y_max, x_max);
            //if(pointcounter != 0)
            //   topoMatrix(y_axis, x_axis) /= pointcounter; //normalisation to number of known points
        }
    }


//-------------------------------------------------
    /*//y-x backward
    for(qint32 y_axis = topoMatrix.rows()-1; y_axis > 0; y_axis--) //y_axis to interpolate among this axis
    {
        for(qint32 x_axis = topoMatrix.cols()-1; x_axis > 0; x_axis--) //x_axis to interpolate among this axis
        {
            qint32 rows = topoMatrix.rows();
            qint32 cols = topoMatrix.cols();

            qint32 pointcounter = 0; //number of known points used to calc interpolation
            qint32 y= 0;
            qint32 x = 0;

            for(qint32 i = 0; i < coors.length(); i++)
            {
                x = coors[i].x(); //known Coordinate x
                y = coors[i].y(); // known coordinate y
                if(!(x==x_axis && y==y_axis))
                {
                    qreal scalar = 1 - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) / sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) ;

                    if (abs(scalar) > 0.85)
                    {
                        result_yxb(y_axis, x_axis) += scalar*topoMatrix(y,x);
                        qreal test = topoMatrix(y,x);
                        pointcounter++;
                    }
                }
            }
            if(pointcounter != 0)
               topoMatrix(y_axis, x_axis) /= pointcounter; //normalisation to number of known points
        }
    }

    //x-y backward
    for(qint32 x_axis = topoMatrix.cols()-1; x_axis > 0; x_axis--) //x_axis to interpolate among this axis
    {
        for(qint32 y_axis = topoMatrix.rows()-1; y_axis > 0 ; y_axis--) //y_axis to interpolate among this axis
        {
            qint32 rows = topoMatrix.rows();
            qint32 cols = topoMatrix.cols();

            qint32 pointcounter = 0; //number of known points used to calc interpolation
            qint32 y = 0;
            qint32 x = 0;

            for(qint32 i = 0; i < coors.length(); i++)
            {
                x = coors[i].x(); //known Coordinate x
                y = coors[i].y(); // known coordinate y
                if(!(x==x_axis && y==y_axis))
                {
                    qreal scalar = 1 - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) / sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) ;

                    if (abs(scalar) > 0.85)
                    {
                        result_xyb(y_axis, x_axis) += scalar*topoMatrix(y,x);
                        qreal test = topoMatrix(y,x);
                        pointcounter++;
                    }
                }
            }
            if(pointcounter != 0)
               topoMatrix(y_axis, x_axis) /= pointcounter; //normalisation to number of known points
        }
    }
*/
    return topoMatrix;
}

//*************************************************************************************************************

MatrixXd Tpplot::calcBilinearInterpolation(MatrixXd gridPointMatrix, QMap<QString, QPoint> mapGrid)
{
    QList<QPoint> coors = mapGrid.values();
    MatrixXd topoMatrix = MatrixXd::Zero(gridPointMatrix.rows(), gridPointMatrix.cols());

    /*
    QList<InterpolationInputData> lData;
    int iThreadSize = QThread::idealThreadCount()*2;
    int iStepsSize = topoMatrix.rows()/iThreadSize;
    int iResidual = topoMatrix.rows()%iThreadSize;

    InterpolationInputData dataTemp;
    dataTemp.maInputData = topoMatrix;
    dataTemp.coors = coors;

    for (int i = 0; i < iThreadSize; ++i)
    {
        dataTemp.iRangeLow = i*iStepsSize;
        dataTemp.iRangeHigh = i*iStepsSize+iStepsSize;
        lData.append(dataTemp);
    }

    dataTemp.iRangeLow = iThreadSize*iStepsSize;
    dataTemp.iRangeHigh = iThreadSize*iStepsSize+iResidual;
    lData.append(dataTemp);

    QFuture<MatrixXd> resultMat = QtConcurrent::mappedReduced(lData, computeyxf, reduce);
    resultMat.waitForFinished();
    topoMatrix = resultMat.result();
    */

    //qint32 size = gridPointMatrix.rows() * gridPointMatrix.cols();
    //        qint32 y= 0;
    //        qint32 x = 0;
    //somehow bilininear
    /*
    for(qint32 i = 0; i < size; i++) //y_axis to interpolate among this axis
    {
        qint32 x_axis = i % gridPointMatrix.cols();
        qint32 y_axis = i /gridPointMatrix.cols();
        if(x_axis == 0)
        {
           y= 0;
           x = 0;
        }

        //don´t interpolate already given points
        //if(coors.indexOf(QPoint(x_axis, y_axis)) > 0)
        //    continue;

        for(qint32 i = 0; i < coors.length(); i++)
        {
            x = coors[i].x(); //known Coordinate x
            y = coors[i].y(); // known coordinate y

            qreal scalar = 1 - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) / sqrt( pow(gridPointMatrix.cols(), 2) + pow(gridPointMatrix.rows(), 2)) ;
            //qreal scalar = sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) ;

            topoMatrix(y_axis, x_axis) += pow(scalar, 20) * gridPointMatrix(y,x);
            //yxf(y_axis, x_axis) += pow(1-scalar,-2)*topoMatrix(y,x); decreasing is too fast
            //yxf(y_axis, x_axis) += 1/pow(scalar,2)  *topoMatrix(y,x);
            // * exp(-10*pow((scalar-1),2))
            //qreal test = topoMatrix(y,x);
        }
        //if(coors.length() != 0)
        //   topoMatrix(y_axis, x_axis) /= coors.length(); //normalisation to number of known points
    }
    */

    for(qint32 y_axis = 0; y_axis < gridPointMatrix.rows(); y_axis++) //y_axis to interpolate among this axis
    {
        for(qint32 x_axis = 0; x_axis < gridPointMatrix.cols(); x_axis++) //x_axis to interpolate among this axis
        {           
            qint32 y= 0;
            qint32 x = 0;

            //don´t interpolate already given points
            //if(coors.indexOf(QPoint(x_axis, y_axis)) > 0)
            //    continue;

            for(qint32 i = 0; i < coors.length(); i++)
            {
                x = coors[i].x(); //known Coordinate x
                y = coors[i].y(); // known coordinate y

                qreal scalar = 1 - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) / sqrt( pow(gridPointMatrix.cols(), 2) + pow(gridPointMatrix.rows(), 2)) ;
                //qreal scalar = sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) ;

                topoMatrix(y_axis, x_axis) += pow(scalar, 20) * gridPointMatrix(y,x);
                //yxf(y_axis, x_axis) += pow(1-scalar,-2)*topoMatrix(y,x); decreasing is too fast
                //yxf(y_axis, x_axis) += 1/pow(scalar,2)  *topoMatrix(y,x);
                // * exp(-10*pow((scalar-1),2))
                //qreal test = topoMatrix(y,x);
            }
            //if(coors.length() != 0)
            //   topoMatrix(y_axis, x_axis) /= coors.length(); //normalisation to number of known points
        }
    }
    //qreal minCoeff = topoMatrix.minCoeff();
    //topoMatrix.array() += abs(topoMatrix.minCoeff());// + topoMatrix;

    return topoMatrix;
}

//*************************************************************************************************************

MatrixXd Tpplot::computeyxf(const InterpolationInputData& inputData)
{
    MatrixXd topoMatrix = MatrixXd::Zero(inputData.maInputData.rows(), inputData.maInputData.cols()); //inputData.maInputData;
    //somehow bilininear
    for(qint32 y_axis = inputData.iRangeLow; y_axis < inputData.iRangeHigh; y_axis++) //y_axis to interpolate among this axis
    {
        for(qint32 x_axis = 0; x_axis < inputData.maInputData.cols(); x_axis++) //x_axis to interpolate among this axis
        {
            qint32 y= 0;
            qint32 x = 0;
            //don´t interpolate already given points
            //if(coors.indexOf(QPoint(x_axis, y_axis)) > 0)
            //    continue;

            for(qint32 i = inputData.iRangeLow; i < inputData.iRangeHigh; i++)
            {
                x = inputData.coors[i].x(); //known Coordinate x
                y = inputData.coors[i].y(); // known coordinate y

                qreal scalar = 1 - sqrt(pow((x - inputData.x), 2) + pow((y - inputData.y), 2)) / sqrt( pow(inputData.maInputData.cols(), 2) + pow(inputData.maInputData.rows(), 2)) ;
                //qreal scalar = sqrt( pow(topoMatrix.cols(), 2) + pow(topoMatrix.rows(), 2)) - sqrt(pow((x - x_axis), 2) + pow((y - y_axis), 2)) ;
                topoMatrix(inputData.y, inputData.x) += pow(scalar, 20) * inputData.maInputData(y,x);
                //yxf(y_axis, x_axis) += pow(1-scalar,-2)*topoMatrix(y,x); decreasing is too fast
                //yxf(y_axis, x_axis) += 1/pow(scalar,2)  *topoMatrix(y,x);
                // * exp(-10*pow((scalar-1),2))
                //qreal test = topoMatrix(y,x);
            }
            //if(coors.length() != 0)
            //   topoMatrix(y_axis, x_axis) /= coors.length(); //normalisation to number of known points
        }
    }
    //qreal minCoeff = topoMatrix.minCoeff();
    //topoMatrix.array() += abs(topoMatrix.minCoeff());// + topoMatrix;

    return topoMatrix;
}

//*************************************************************************************************************

void Tpplot::reduce(MatrixXd &resultData, const MatrixXd &data)
{
    if(resultData.size() == 0)
        resultData = data;
    else
        resultData += data;
}

















