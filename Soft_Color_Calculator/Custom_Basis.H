#ifndef __Custom_Basis_H
#define __Custom_Basis_H

#include <vector>
#include <iostream>
#include "CBasis.H"
#include "CBasis.C"
#include "CMatrix.H"

class Custom_Basis
{
private:

public:
	Custom_Basis();
	~Custom_Basis();
	abasis manual_Cbasis();
};

void add_element(CBasis<int>& toadd, abasis& res);



#endif
