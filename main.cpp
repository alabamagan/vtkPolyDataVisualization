/*
 * vtk poly data visualization module
 */
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include "vtkCamera.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkPolyDataMapper.h"
#include "vtkCellData.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "vtkActor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkWindowToImageFilter.h"
using namespace std;

static void PrintUsage(char* argv[]) {
    std::cerr << "Usage: " << argv[0] << " -i <vtkfile> [-c <camera>|-a <azimuth>|-e <elevation>|-s 400x400|-t <x,y>|-o <./tmp.jpg>]\n"
    << "Options:\n"
    << "\t-h,--help\t\tShow this help message\n"
    << "\t-i,--input\t\tDirectory to input vtk data file.\n"
    << "\t-a,--azimuth\t\tAzimuthal rotation applied to the model along the axis that run in the direction of view up angle through the focal point\n"
    << "\t-e,--elevation\t\tElevation rotation applied to the model.\n"
    << "\t-c,--camera\t\tPoints to the file that stores previous camera information\n"
    << "\t-t,--translation\tTranslation factor, inform of 'X,Y', defualt to 0,0, which corresponds to translation along x-axis and y-axis. Rightwards and Upwards is positive. Note that translation is done after rotation.\n"
    << "\t-o,--output\t\tOutput image directory, default to [./tmp.jpg]\n"
    << std::endl;
};

static vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

//double* crossProduct(double* x, double* y) {
//    double a = x[1] * y[2] - x[2] * y[1];
//    double b = x[2] * y[0] - x[0] * y[2];
//    double c = x[0] * y[1] - x[1] * y[0];
//    cout <<"CrossProduct: " << a <<","<< b <<"," <<c << endl;
//    double* out = (double*)malloc(3*sizeof(double));
//    out[0] = a;
//    out[1] = b;
//    out[2] = c;
//    return out;
//}
//
//double* dotProduct(double* x, double* y) {
//    double a = x[0]*y[0];
//    double b = x[1]*y[1];
//    double c = x[2]*y[2];
//    cout <<"dotProduct: " << a <<","<< b <<"," <<c << endl;
//    double* out = (double*)malloc(3*sizeof(double));
//    out[0] = a;
//    out[1] = b;
//    out[2] = c;
//    return out;
//}
//
//double* dotScalarToNormal(double* x, double C) {
//    double a = x[0]*C;
//    double b = x[1]*C;
//    double c = x[2]*C;
//    cout <<"dotScalarToNormal: " << a <<","<< b <<"," <<c << endl;
//    double* out = (double*)malloc(3*sizeof(double));
//    out[0] = a;
//    out[1] = b;
//    out[2] = c;
//    return out;
//}
//
//double* addVectors(double* x, double* y) {
//    double a = x[0] + y[0];
//    double b = x[1] + y[1];
//    double c = x[2] + y[2];
//    cout <<"addVectors: " << a <<","<< b <<"," <<c << endl;
//    double* out = (double*)malloc(3*sizeof(double));
//    out[0] = a;
//    out[1] = b;
//    out[2] = c;
//    return out;
//}
//
//static void CameraTrans(vtkCamera* camera, char* translationString) {
//    double* curpos = camera->GetPosition();
//    double* curViewPlaneNormal = camera->GetViewPlaneNormal();
//    double* curViewUp = camera->GetViewUp();
//    double *normX, *normY, *transXV, *transYV, *outpos;
//    double transX, transY;
//
//    cout << curpos[0] <<","<< curpos[1] <<"," <<curpos[2] << endl;
//    cout << curViewPlaneNormal[0] <<","<< curViewPlaneNormal[1] <<"," <<curViewPlaneNormal[2] << endl;
//    cout << curViewUp[0] <<","<< curViewUp[1] <<"," <<curViewUp[2] << endl;
//
//    vector<string> translationVect;
//
//    normX = crossProduct(curViewUp, curViewPlaneNormal);
//    cout << "normX: " << normX[0] <<","<< normX[1] <<"," <<normX[2] << endl;
//    normY = crossProduct(normX, curViewPlaneNormal);
//    cout << "normY: " << normY[0] <<","<< normY[1] <<"," <<normY[2] << endl;
//    split(translationString, ',', translationVect);
//    transX = atof(translationVect[0].c_str());
//    transY = atof(translationVect[1].c_str());
//    transXV = dotScalarToNormal(normX, transX);
//    transYV = dotScalarToNormal(normY, transY);
//
//    outpos = addVectors(transXV, curpos);
//    outpos = addVectors(transYV, outpos);
//
//    cout << "transX,Y: " << transX <<","<<transY << endl;
//    cout << "transXV: " << transXV[0] <<","<< transXV[1] <<"," <<transXV[2] << endl;
//    free(normX);
//    free(normY);
//    free(transXV);
//    free(transYV);
//    camera->SetPosition(outpos);
//}

static void CameraTrans(vtkCamera* camera, double transX, double transY) {
//    double transX, transY;
//    vector<string> translationVect;
//    split(translationString, ',', translationVect);
//    transX = atof(translationVect[0].c_str());
//    transY = atof(translationVect[1].c_str());

    vtkMatrix4x4* transMatrix = vtkMatrix4x4::New();
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                transMatrix->SetElement(i,j,1.);
            } else {
                transMatrix->SetElement(i,j,0.);
            }

        }
    }

    transMatrix->SetElement(0,3, -transX);
    transMatrix->SetElement(1,3, -transY);

    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix(transMatrix);

    camera->ApplyTransform(transform);
}

static void paramFileParser(char* filename, double* output) {
    string line;
    vector<string> l1, l2;
    ifstream paramFile;
    paramFile.open(filename);
    if (paramFile.is_open()) {
        while (getline(paramFile, line)) {
            l1.push_back(line);
        }
    } else {
        cerr << "Unable to open file" << endl;
    }
    split(l1[0], ',',l2); // position
    split(l1[1], ',',l2); // focal point
    split(l1[2], ',',l2); // view up
    split(l1[3], ',',l2); // translation

    for(int i = 0; i < 11;i++) {
        output[i] = atof(l2[i].c_str());
    }
    output[11] = atof(l1[4].c_str());
}

int main(int argc, char* argv[]) {
    char* inputFile = (char*) "/home/lwong/Documents/netbeanProjects/vtkVisualization/TestData/tract5000.vtk";
    char* cameraParametersFile = NULL;
    char* cameraParametersOutDir = (char*) "./camera.txt";
    char* sizeString = (char*) "400x400";
    char* translationString = (char*) "0,0";
    char* outputName = (char*) "tmp.jpg";
    int sizeX, sizeY;
    double zoomFactor = 1.0;
    double azimuth = 0;
    double elevation = 0;
    double output[12];
    double transX, transY;

    if (argc < 1) {
        PrintUsage(argv);
        return 1;
    }
    for (int i = 0; i < argc; i++) {
        char* arg = argv[i];
        if (strcmp(arg, "-h") == 0|| strcmp(arg,"--help")==0) {
            PrintUsage(argv);
            return 1;
        } else if (strcmp(arg, "-i")==0 || strcmp(arg, "--input") == 0) {
            if (i+1 < argc) {
                inputFile = argv[++i];
            } else {
                cerr << "No arguments to option -i" << endl;
            }
        } else if (strcmp(arg, "-c")==0 || strcmp(arg, "--camera")==0) {
            if (i+1 < argc) {
                cameraParametersFile = argv[++i];
            } else {
                cerr << "No arguments to option -c/--camera" << endl;
            }
        } else if (strcmp(arg, "-s")==0 || strcmp(arg, "--size")==0) {
            if (i+1 < argc) {
                sizeString = argv[++i];
            } else {
                cerr << "No arguments to option -s/--size" << endl;
            }
        } else if (strcmp(arg, "-z")==0 || strcmp(arg, "--zoom")==0) {
            if (i+1 < argc) {
                zoomFactor = atof(argv[++i]);
            } else {
                cerr << "No arguments to option -z/--zoom" << endl;
            }
        } else if (strcmp(arg, "-a")==0 || strcmp(arg, "--azimuthRotation")==0) {
            if (i+1 < argc) {
                azimuth = atof(argv[++i]);
            } else {
                cerr << "No arguments to option -a/--AzimuthRotation" << endl;
            }
        } else if (strcmp(arg, "-e")==0 || strcmp(arg, "--elevationRotation")==0) {
            if (i+1 < argc) {
                elevation = atof(argv[++i]);
            } else {
                cerr << "No arguments to option -z/--zoom" << endl;
            }
        } else if (strcmp(arg, "-o")==0 || strcmp(arg, "--output")==0) {
            if (i+1 < argc) {
                outputName = argv[++i];
            } else {
                cerr << "No arguments to option -o/--output" << endl;
            }
        } else if (strcmp(arg, "-t")==0 || strcmp(arg, "--translation")==0) {
            if (i + 1 < argc) {
                translationString = argv[++i];
            } else {
                cerr << "No arguments to option -t/--translation" << endl;
            }

        }
    };

    // Initialize parameter from user input
    vector<string> elems;
    split(sizeString, 'x', elems);
    split(translationString, ',', elems);
    sizeX = (int) atof(elems[0].c_str());
    sizeY = (int) atof(elems[1].c_str());
    transX = atof(elems[2].c_str());
    transY = atof(elems[3].c_str());


    if (cameraParametersFile != NULL) {
        paramFileParser(cameraParametersFile, output);
    }

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(inputFile);
    reader->GetOutput()->GetCellData()->SetActiveScalars("Average Direction");
    
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection(reader->GetOutputPort());
    
    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);
    
    vtkRenderer* renderer = vtkRenderer::New();
    renderer->AddActor(actor);

    vtkRenderWindow* renWin = vtkRenderWindow::New();
    renWin->AddRenderer(renderer);
    renWin->SetOffScreenRendering(1);
    renWin->SetSize(sizeX, sizeY);
    renWin->Render();
    renWin->SetAAFrames(10);

    vtkCamera* camera =renderer->GetActiveCamera();
    if (cameraParametersFile != NULL) {
        // Translation handling

        camera->SetPosition(output[0], output[1], output[2]);
        camera->SetFocalPoint(output[3], output[4], output[5]);
        camera->SetViewUp(output[6],output[7], output[8]);
        transX += output[9];
        transY += output[10];
        zoomFactor = zoomFactor*output[11];

    }
    camera->Zoom(zoomFactor);
    camera->Azimuth(azimuth);
    camera->Elevation(elevation);
    camera->OrthogonalizeViewUp();
    renWin->Render();

    if (translationString != "0,0") {
        CameraTrans(camera, transX, transY);
    }

    vtkWindowToImageFilter* wintoim = vtkWindowToImageFilter::New();
    wintoim->SetInput(renWin);
    wintoim->Update();

    vtkJPEGWriter* writer = vtkJPEGWriter::New();
    writer->SetFileName(outputName);
    writer->SetInputConnection(wintoim->GetOutputPort());
    writer->Write();

    // Write CameraFile
    double* t1 = camera->GetPosition();
    double* t2 = camera->GetFocalPoint();
    double* t3 = camera->GetViewUp();
    char s1[100], s2[100], s3[100], s4[100], s5[100];
    sprintf(s1, "%f,%f,%f", t1[0], t1[1], t1[2]);
    sprintf(s2, "%f,%f,%f", t2[0], t2[1], t2[2]);
    sprintf(s3, "%f,%f,%f", t3[0], t3[1], t3[2]);
    sprintf(s4, "%f,%f", transX, transY);
    sprintf(s5, "%f", zoomFactor);
    ofstream paramOut;
    paramOut.open(cameraParametersOutDir);
    if (paramOut.is_open()) {
        paramOut << s1 <<endl;
        paramOut << s2 <<endl;
        paramOut << s3 <<endl;
        paramOut << s4 << endl;
        paramOut << s5 << endl;
    }
    paramOut.close();

    return 0;
    }
