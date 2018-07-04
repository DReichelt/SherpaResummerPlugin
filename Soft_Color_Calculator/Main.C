#include <termios.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include <time.h>
#include "Color.H"
#include "Tools.H"
#include "CBasis.H"
#include "CMatrix.C"
#include "Auto_Basis.C"
#include "Custom_Basis.C"
#ifdef MALLOC_TRACE
#include <mcheck.h>
#endif

using namespace ATOOLS;

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int main(int argc, char **)
{
    argc = argc + 0;
    clock_t t1,t2;
    t1=clock();

    //--------------------------------------------
    //basis input parameters
    double Nc = 3.;
    int num_quark_pairs = 0;
    int num_gluons = 4;
    //--------------------------------------------

    //--------------------------------------------
    //output info
    int anti_quarks_plus_quarks = 2*num_quark_pairs;
    int numlegs = num_gluons + anti_quarks_plus_quarks;
    std::string filename = "";
    for(unsigned i = 0; i < double(num_quark_pairs); i++) filename = filename + "qqb";
    for(unsigned i = 0; i < num_gluons; i++) filename = filename + "g";
    //--------------------------------------------

    //--------------------------------------------
    print_cl_info(filename);
    //--------------------------------------------

    //automatic basis contructor, 'true' gives adjoint (hep-ph/9910563) basis
    //true gives trace basis
    //--------------------------------------------
    abasis NBasis;
    Custom_Basis cb;
    //--------------------------------------------
    //Enter custom basis information here
    //--------------------------------------------
    //NBasis = cb.manual_Cbasis();
    //--------------------------------------------
    //--------------------------------------------
    NBasis = auto_Cbasis(num_gluons,num_quark_pairs,false);
    int bdim = NBasis.size();
    normalize_basis(NBasis, Nc);
    printBasis(NBasis);
    //--------------------------------------------


    //compute color matrices
    //--------------------------------------------
    //--------------------------------------------
    CMatrix Cmetric;
    std::vector< CMatrix > Tprods;
    int numprods = numlegs*(numlegs-1)/2;
    compute_metric(Cmetric,NBasis);
    compute_color_product(NBasis,Tprods,num_quark_pairs,num_gluons,numprods,Nc);
    //--------------------------------------------
    //--------------------------------------------


    //print to command line
    //--------------------------------------------
    //--------------------------------------------
    std::cout << "Metric:" << endl;
    Cmetric.printMat(bdim);
    std::cout << "Tprods:" << endl;
    int lin_iterator = 0;
    for(unsigned i = 0; i<numlegs; i++)
    {
        for(unsigned j = i+1; j<numlegs; j++)
        {
            cout << "T" << i+1 << ".T" << j+1 << std::endl;
            Tprods[lin_iterator].printMat(bdim);
            lin_iterator++;
        }
    }
    //--------------------------------------------
    //--------------------------------------------


    //--------------------------------------------
    //--------------------------------------------
    check_color_conservation(Tprods,Cmetric,num_gluons,num_quark_pairs,bdim);
    std::cout << "Dimension of basis: " << bdim << std::endl;
    //--------------------------------------------
    //--------------------------------------------


    //--------------------------------------------
    //--------------------------------------------
    std::ofstream outfileMetric;
    outfileMetric.open(filename+"_met.dat", std::ios_base::out);
    Cmetric.printFile(outfileMetric,bdim);
    outfileMetric << bdim << std::endl;
    outfileMetric.close();
    //--------------------------------------------
    //--------------------------------------------


    //--------------------------------------------
    //--------------------------------------------
    std::ofstream outfile;
    outfile.open(filename+".dat", std::ios_base::out);
    for(int i(0); i<numprods ; i++) Tprods[i].printFile(outfile,bdim);
    outfile << bdim << std::endl;
    outfile.close();
    //--------------------------------------------
    //--------------------------------------------


    //--------------------------------------------
    Cmetric.cleanup();
    for(int i(0); i < numprods; i++) Tprods[i].cleanup();
    //--------------------------------------------


    t2=clock();
    float diff ((float)t2-(float)t1);
    float seconds = diff / CLOCKS_PER_SEC;
    std::cout << "time: " << seconds << std::endl;

    return 0;
}
//--------------------------------------------------------------------
//--------------------------------------------------------------------




