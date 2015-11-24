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
#include "vtkActor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkWindowToImageFilter.h"
using namespace std;

static void PrintUsage(char* argv[]) {
    std::cerr << "Usage: " << argv[0] << " -i <vtkfile> [-c <camera>|-a <azimuth>|-e <elevation>|-s 400x400|-o <./tmp.jpg>]\n"
    << "Options:\n"
    << "\t-h,--help\t\tShow this help message\n"
    << "\t-i,--input\t\tDirectory to input vtk data file.\n"
    << "\t-a,--azimuth\t\tAzimuthal rotation applied to the model along the axis that run in the direction of view up angle through the focal point\n"
    << "\t-e,--elevation\t\tElevation rotation applied to the model.\n"
    << "\t-c,--camera\t\tPoints to the file that stores previous camera information\n"
    << "\t-o,--output\t\tOutput image directory, default to [./tmp.jpg]"
    << std::endl;
};

static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
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

    for(int i = 0; i < 9;i++) {
        output[i] = atof(l2[i].c_str());
    }
    output[9] = atof(l1[3].c_str());
}

int main(int argc, char* argv[]) {
    char* inputFile = "/home/lwong/Documents/netbeanProjects/vtkVisualization/TestData/tract5000.vtk";
    char* cameraParametersFile = NULL;
    char* cameraParametersOutDir = "./camera.txt";
    char* sizeString = "400x400";
    char* outputName = "tmp.jpg";
    int sizeX, sizeY;
    double zoomFactor = 1.0;
    double azimuth = 0;
    double elevation = 0;
    double output[10];

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
        }


    };

    // Initialize parameter
    vector<string> elems;
    split(sizeString, 'x', elems);
    sizeX = (int) atof(elems[0].c_str());
    sizeY = (int) atof(elems[1].c_str());

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
        camera->SetPosition(output[0], output[1], output[2]);
        camera->SetFocalPoint(output[3], output[4], output[5]);
        camera->SetViewUp(output[6],output[7], output[8]);
        zoomFactor = zoomFactor*output[9];
    }
    camera->Zoom(zoomFactor);
    camera->Azimuth(azimuth);
    camera->Elevation(elevation);
    camera->OrthogonalizeViewUp();
//    renWin->Render();

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
    char s1[100];
    char s2[100];
    char s3[100];
    char s4[100];
    sprintf(s1, "%f,%f,%f", t1[0], t1[1], t1[2]);
    sprintf(s2, "%f,%f,%f", t2[0], t2[1], t2[2]);
    sprintf(s3, "%f,%f,%f", t3[0], t3[1], t3[2]);
    sprintf(s4, "%f", zoomFactor);
    ofstream paramOut;
    paramOut.open(cameraParametersOutDir);
    if (paramOut.is_open()) {
        paramOut << s1 <<endl;
        paramOut << s2 <<endl;
        paramOut << s3 <<endl;
        paramOut << s4 << endl;
    }
    paramOut.close();

    return 0;
    }
