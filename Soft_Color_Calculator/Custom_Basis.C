#ifndef __Custom_Basis_C
#define __Custom_Basis_C
#include "Custom_Basis.H"



    Custom_Basis::Custom_Basis(){}
    Custom_Basis::~Custom_Basis(){}


abasis Custom_Basis::manual_Cbasis()
{

    //enter custom basis information here
    //and comment out get_Cbasis line previous
    //do not forget to set the flavors accordingly
    //so that the color product calculation knows 
   	//which exchange operators to assign
    
    abasis results;
    CBasis<int> qel(1,0);
    
    qel.Tadd(1,51,1);
    qel.Tadd(2,1,52);
    add_element(qel,results);
    
    qel.Tadd(2,51,1);
    qel.Tadd(1,1,52);
    add_element(qel,results);
    
    qel.Dadd(51,52);
    qel.Tadd(1,1,2);
    qel.Tadd(2,2,1);
    add_element(qel,results);
    
    return results;
}

void add_element(CBasis<int>& toadd, abasis& res){

    std::vector< CBasis<int> > tbas;
    tbas.push_back(toadd);
    res.push_back(tbas);
    tbas.clear(); 
    toadd.clear();    

}

#endif
