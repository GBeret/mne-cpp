
//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtGui>
#include <QApplication>
#include <QCommandLineParser>
#include <QWidget>
#include <QLabel>
#include <QVBoxLayout>

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <fiff/fiff.h>
#include <utils/layoutloader.h>
#include <disp/plots/topoplot.h>

//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace FIFFLIB;
using namespace DISPLIB;
using namespace UTILSLIB;

//*************************************************************************************************************
//=============================================================================================================
// MAIN
//=============================================================================================================



//=============================================================================================================
/**
* The function main marks the entry point of the program.
* By default, main has the storage class extern.
*
* @param [in] argc (argument count) is an integer that indicates how many arguments were entered on the command line when the program was started.
* @param [in] argv (argument vector) is an array of pointers to arrays of character objects. The array objects are null-terminated strings, representing the arguments that were entered on the command line when the program was started.
* @return the value that was set to exit() (which is 0 if exit() is called via quit()).
*/
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QWidget w;
    QLabel * topoLabel = new QLabel();
    w.resize(800, 800);
    QVBoxLayout * layout = new QVBoxLayout();
    layout->addWidget(topoLabel);
    w.setLayout(layout);
    w.show();

    // Command Line Parser
    QCommandLineParser parser;
    parser.setApplicationDescription("Topoplot Example");
    parser.addHelpOption();

    QCommandLineOption evokedFileOption("ave", "Path to the evoked/average <file>.", "file", QCoreApplication::applicationDirPath() + "/MNE-sample-data/sample_audvis-ave.fif");
    QCommandLineOption evokedIdxOption("aveIdx", "The average <index> to choose from the average file.", "index", "1");

    parser.addOption(evokedFileOption);
    parser.addOption(evokedIdxOption);

    parser.process(a);

    //generate FiffEvoked object
    QFile t_sampleFile(parser.value(evokedFileOption));
    FiffEvoked p_FiffEvoked(t_sampleFile,QVariant(parser.value(evokedIdxOption)));

    MatrixXd signaltemp  = MatrixXd::Zero(p_FiffEvoked.data.cols(), p_FiffEvoked.data.rows());
    for(qint32 channels = 0; channels <  p_FiffEvoked.data.rows(); channels++)
        signaltemp.col(channels) = p_FiffEvoked.data.row(channels);

    MatrixXd signal = signaltemp.block(0, 0, 1, 306);
    QSize topoSize(64, 64); //(256, 256);//(64, 64);

    QString selectionName("Vectorview-all.lout");
    QString path = QCoreApplication::applicationDirPath() + selectionName.prepend("/Resources/2DLayouts/");
    QMap<QString,QPointF> layoutMap;
    LayoutLoader::readMNELoutFile(path, layoutMap);

    //topoplot
    TopoPlot *create_plot = new TopoPlot();
    QList<QImage> imageList = create_plot->createTopoPlotImageList(signal, layoutMap, topoSize, QSize(800, 800), Jet, 20);
    topoLabel->setPixmap(QPixmap::fromImage(imageList.at(0)));
    w.repaint();

    return a.exec();
}
