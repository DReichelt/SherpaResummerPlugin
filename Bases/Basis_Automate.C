#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include "Tools/CBasis.H"
#include "Bases/Perm_Helpers.C"


using namespace RESUM;

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// start of basis functions
// input - Number of gluons, quark, a-quarks
// output - an ABasis representing the normal
// ordered basis for this process...ie with all 
// quarks to left of aquarks left of gluons
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//Inlines
//------------------------------------------------------------------------

inline std::vector< CBasis <int> > combineA(std::vector< CBasis <int> > rhs1 , std::vector< CBasis <int> > rhs2 ){
  std::vector< CBasis <int> > results;
  CBasis <int> plan(0,0);
  
  for(unsigned i = 0; i < rhs1.size(); i++){
    for(unsigned j = 0; j < rhs2.size(); j++){
      plan = rhs1[i].cadd(rhs2[j]);
      results.push_back(plan);  
    }
  }
  return results;
}

//Fold 2 ABasis elements into 1
inline ABasis fold(ABasis rhs1, ABasis rhs2){
 ABasis results;
 if(rhs1.size() == 0) return rhs2;
 if(rhs2.size() == 0) return rhs1;
 
 int ind_hold = 0;
 results.resize(rhs1.size()*rhs2.size());
 
 for(unsigned i1 = 0; i1 < rhs1.size(); i1++){
    for(unsigned i2 = 0; i2 < rhs2.size(); i2++){
          results[ind_hold]=combineA(rhs1[i1],rhs2[i2]);
	  ind_hold++;
       }
    }
 return results;
}

//Returns an adjoint delta function
inline CBasis<int>   AdjointDelta(int a, int b){
    CBasis<int> results(1,0);
    results.Tadd(a,a,b);
    results.Tadd(b,b,a);
    return results;
}


//Get symmetry factor of a partition
inline void  get_partition_symmetry_factor(std::vector< int > partition, size_t &sym_factor, size_t &sym_s){
  int last = 0;
  for(unsigned it = 0; it < partition.size(); ++it ){
        if(last == partition[it]){sym_factor++; sym_s = partition[it];}
        last = partition[it];
        }
}


//get groupings
inline std::vector< std::string > arrange_groups( std::vector<int> parts, std::string labels, int is_quarks){
  std::vector< std::string > results;
  size_t sum_of_c = 0, last , sym_s = 0; 
      size_t sym_fact=1;
      std::vector<int> cumsum_of_c;

  //get all index groupings
  if(is_quarks == 0) get_partition_symmetry_factor( parts , sym_fact, sym_s );

  //determine the symmetry factor for the partition
  for(unsigned it = 0; it < parts.size(); ++it ){
  sum_of_c += parts[it];
  cumsum_of_c.push_back(sum_of_c);
}

  //Get all index groupings for given partitions
  std::vector<int> order(sum_of_c); 
  std::vector<char> usage(sum_of_c);
  allo_discon_groups_impl(usage, order, cumsum_of_c, 0, 0, 0, sym_fact, sym_s, results, labels);

  return results;
}


//------------------------------------------------------------------------
//End of Inlines
//------------------------------------------------------------------------

CBasis<int> ToAdjoint(std::string rhs){
  CBasis <int>  results(1,0);
  int int_start=501;
  int nsize = rhs.size();
  results.Fadd(rhs[0]-'0',rhs[1]-'0',int_start);
  results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());

        //loop over internal f's (n > 4)
	for(unsigned m = 0; m < nsize-4; m++){
	results.Fadd(int_start,rhs[m+2]-'0',int_start+1);
	results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());
	int_start++;
        }

  results.Fadd(int_start,rhs[nsize-2]-'0',rhs[nsize-1]-'0');
  results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());

  return results;
}



std::vector < CBasis<int> > ToFundamental(char n, std::string rhs, int nquarks, std::string qs,  std::string aqs){
  std::vector< CBasis <int> > results;
  int ssize = rhs.size();
  int index;
  int sign = pow(-1,ssize); //Relative sign in Reflection term
  CBasis<int> flow(2*sqrt(2),0),flowR(2*sign*sqrt(2),0);
  flow.Tin.resize(ssize);
  flowR.Tin.resize(ssize);
  
  //Gluon indices
  for(unsigned i = 0; i<ssize; i++){
    flow.Tin[i].resize(3); 
    flow.Tin[i][1] = rhs[i]-'0'; flow.Tin[i][0] = rhs[i]-'0';
    flow.Tin[(ssize-1+i)%ssize].resize(3); 
    flow.Tin[(ssize-1+i)%ssize][2] = (rhs[i]-'0');
    index = (i == 0) ? 0 : ssize-i;
    flowR.Tin[index].resize(3); 
    flowR.Tin[index][1] = rhs[i]-'0'; flowR.Tin[index][0] = rhs[i]-'0';
    flowR.Tin[(ssize-1+index)%ssize].resize(3); 
    flowR.Tin[(ssize-1+index)%ssize][2] = (rhs[i]-'0');    
    }
   
  //Add quark indices
  if(nquarks == 2)
    {
    int qstart = 51;
    flow.Tin[0][1] = qstart;
    flow.Tin[ssize-1][2] = -qstart-1;
    qstart++;
    }

  results.push_back(flow);
  if(nquarks == 0) results.push_back(flowR);
  
  return results;
}



//Returns f-basis elements
ABasis AdjointConGluons(std::string rhs)
{
   ABasis results;
   std::vector< CBasis<int> > presults;
   std::vector< std::string > Perms;
   std::string base, baseT;  
   int nsize = rhs.size();
   baseT = rhs;
           	
	   //f-basis element (Ng >= 4)
	   //permute (n-2) string
           base = baseT.substr(1,nsize-2);

	   //get (n-1)! permutations in lexographic order
	   Perms = perms(base);
	     
	   //sum over permutations   
	   int perm_size = factorial(nsize-2);
	   for(unsigned j = 0; j < perm_size; j++){
	     base = Perms[j]; 
	     baseT = baseT[0]+base+baseT[nsize-1];

	     //Fbasis elements with internal indices open
	     presults.push_back(ToAdjoint(baseT));
	     results.push_back(presults);

	     presults.clear();
	     }
  
    return results;
}


//Returns trace-basis elements
ABasis TraceConGluons( std::string rhs, int num_quarks, std::string ql, std::string aql)
{
   ABasis results;
   std::vector< CBasis<int> > presults;
   std::vector< std::string > Perms;
   std::string base, baseT;  
   int nsize = rhs.size();
   baseT = rhs;
           
   if(rhs.size() == 0 && num_quarks>0){ 
         CBasis<int> qdis(1,0);
         int qstart = 50;
    	 presults.push_back(qdis.Dadd(qstart+(ql[0]-'0'),qstart+(aql[0]-'0') ) );
	 qstart++;
         results.push_back(presults);
	 return results;
         }

	   //trace-basis element (Ng >= 4)
	   //permute (n-1) string
           if(num_quarks == 0) base = baseT.substr(1,nsize-1);
	      else{base = baseT;}
	   //get (n-1)! (or n!) permutations in lexographic order
	   Perms = perms(base);
	     
	   //sum over permutations   
	   int perm_size = factorial(base.size());
	   for(unsigned j = 0; j < perm_size; j++){
	   base = Perms[j]; 
	   if(num_quarks == 0) baseT = baseT[0]+base;
	      else{baseT = base;}

	   //Cbasis elements with indices contracted
	   //conditional assures that only terms which 
	   //are not reflection related are included 
	   //in the sum.
	   
	   if(baseT[1] < baseT[baseT.size()-1] && num_quarks == 0){
	     presults = ToFundamental('n',baseT,num_quarks,"","");
	   results.push_back(presults);}
	   else if(num_quarks > 0){
	     presults = ToFundamental('n',baseT,num_quarks,ql,aql);
	     results.push_back(presults);}

	   presults.clear();
	   }
    
    return results;
}


ABasis gluon_groups( int ngluon, std::string gluon_labels )
{
    //Pure gluons component of basis depending on
    //returns set of basis elements permuted and 
    //distributed over integer partitions
  
    ABasis results, iresults, iresults2, iresultsT;
    std::vector< CBasis<int> > presults;
    size_t place;

    //get partition of gluons integers
    std::vector< std::vector<int> > numglu = gluon_partitions(ngluon);
    
    //for each partition get elements
    for(unsigned g = 0; g < numglu.size(); g++){
      int Ngroups = numglu[g].size();
      std::vector< std::string > strt;
      std::string base, baseT;
 
      //get all index groupings
      strt = arrange_groups(numglu[g] , gluon_labels , 0);  

      //For each grouping get the corresponding elements
      for(unsigned i = 0; i < strt.size(); i++) { 
        place = 0;
       
      //sum over groups in partition
      for(unsigned k = 0; k < Ngroups; k++){
	int nsize = numglu[g][k];
	baseT = strt[i].substr(place,nsize);
	
	//adjoint delta function for ng = 2 partition
        if(nsize == 2){
	  presults.push_back(AdjointDelta(baseT[0]-'0',baseT[1]-'0'));
	  iresults.push_back(presults);
        }
	else if(nsize == 3){
	  iresults =  TraceConGluons(baseT,0,"","");
            }
        else{
	   //Trace basis connected gluons
	   //if(nsize < 6) 
	   results = TraceConGluons(baseT,0,"","");

	   //adjoint basis gluons (for pure gluons only)
	   //if(nsize == 6) 
	   //iresults = AdjointConGluons(baseT);
	   }
      
        //Fold together the partial elements
        if(k >= 1){
	  iresultsT = fold(iresults,iresults2);
	  iresults2 = iresultsT;
        }
        //the first or only element
        else{iresults2 = iresults;}
        
	//book-keeping
	place += numglu[g][k];
        presults.clear(); iresults.clear();  
      }
    
      //insert into the basis
      results.insert(results.end(), iresults2.begin(), iresults2.end());
      }
    }
    
  return results;
}




ABasis qqng_groups(int ngluon, int nquark, int naquark, std::string gluon_labels, bool only_connected){
    
    ABasis results,iresults,iresultsK,iresultsC,iresultsD,iresultsg;
    std::string quark_labels = "12345678", aquark_labels = "2345", aquark_copy; 
    std::vector< CBasis<int> > presults;
    std::vector< std::string > PermsQ,PermsAQ,PermsG;
    int con_to_q;
    int ncon_to_q;
    
    if(nquark+naquark > 2){
    aquark_labels =  quark_labels.substr((nquark+naquark)/2., (nquark+naquark)/2);
    quark_labels =  quark_labels.substr(0,(nquark+naquark)/2.);
    PermsQ = perms(quark_labels);
    PermsAQ = perms(aquark_labels);
    int qstart = 50;
    gluon_labels =  gluon_labels.substr(0,ngluon);
    if(ngluon > 0) gluon_labels.append("0");
    PermsG = perms(gluon_labels);
    
    int keep = 0;
    int ng,nq;
    int f1,af1,adj;
    for(unsigned i = 0; i < PermsG.size(); i++){
      for(unsigned k = 0; k < PermsQ.size(); k++){
         CBasis<int> quark_elms(1,0);
         keep = 0;
	 ng = ngluon;
	 nq = (nquark+naquark)/2;
	 
	 for(unsigned m = 0; m < PermsG[i].size()+1; m++){
	     aquark_copy = aquark_labels;
	     adj =  (PermsG[i][m]-'0');
	     for(unsigned j = m; j < PermsQ[k].size(); j++){
	       f1  =  PermsQ[k][j]-'0';
	       af1 =  aquark_labels[j] -'0';
	       if(ngluon > 0 && adj == 0 && nq > 0  ){ quark_elms.Dadd(qstart+PermsQ[k][j]-'0', qstart+aquark_labels[j] -'0');
		 nq = nq - 2;
	       }
	       if(ngluon > 0 && adj != 0 && adj != 9){
		 quark_elms.Tadd(adj,qstart+PermsQ[k][j]-'0',qstart+aquark_labels[j] -'0');
		 adj = 9;
	       }
	       if(ngluon == 0) quark_elms.Dadd(qstart+PermsQ[k][j]-'0', qstart+aquark_labels[j] -'0');
	       }
	 }
	 if(quark_elms.get_Tdim() < ngluon){
	   //Add Gluon
	   quark_elms.Tadd( quark_elms.Tin[0][0] == 1 ? 2:1,1,quark_elms.Tin[0][2]);
	   quark_elms.Tin[0][2] = 1;
	 }
       presults.push_back(quark_elms);
       iresults.push_back(presults);
       presults.clear();
       }
     }

    //add disconnected parts
    if(ngluon > 1) {
    iresultsg = gluon_groups( ngluon , gluon_labels );
    iresultsD = qqng_groups(0, nquark, naquark, "", true);
    iresultsD = fold(iresultsD,iresultsg);
    iresults.insert(iresults.end(), iresultsD.begin(), iresultsD.end());
    }
    
    results.insert(results.end(), iresults.begin(), iresults.end());
    }
    
    //Number of gluon connected to a single qqb 
    
    if(ngluon > 0 && nquark+naquark <= 2){
	for(unsigned h = 0; h <= ngluon; h++){
    
    //since the single Tr[T] = 0
    if(h == 1) h++;
    
    //number of connected and disconnected
    std::vector<int> numglu;
    con_to_q = ngluon - h ; 
    ncon_to_q = h;   
    if(con_to_q > 0) numglu.push_back(con_to_q);
    if(ncon_to_q > 0) numglu.push_back(ncon_to_q);

    gluon_labels = "123456789";
    std::vector< std::string > strt;
    
    strt = arrange_groups(numglu , gluon_labels , 1 );
     
       for(unsigned k = 0; k < strt.size(); k++){ 
 
       //cut the indices
       gluon_labels = strt[k].substr(con_to_q,ngluon);
       //The disconnected string
       iresultsD = gluon_groups( ncon_to_q , gluon_labels );
    
       //For connected flows build the straight permutations
       //Return as another ABasis
       iresultsC = TraceConGluons(strt[k].substr(0,con_to_q), 2, quark_labels, aquark_labels);
    
       //Fold together the ABasis
       iresults = fold(iresultsC,iresultsD);
       results.insert(results.end(), iresults.begin(), iresults.end());

        }
      }
    }
    return results;
}



ABasis get_Cbasis(int number_of_gluons, int number_of_quarks, int number_of_aquarks){
  ABasis results;
  ABasis build;

  if(number_of_quarks + number_of_aquarks > 0){
    std::string labels = "123456789";
    build = qqng_groups(number_of_gluons, number_of_quarks, number_of_aquarks, labels, true);
    results.insert(results.end(), build.begin(), build.end());
    return results;
  }
    else if(number_of_quarks + number_of_aquarks == 0){
      std::string gluon_labels = "123456789";
      //std::string gluon_labels = "165432";
    build = gluon_groups(number_of_gluons,gluon_labels);
    results.insert(results.end(), build.begin(), build.end());
    return results;
  }
  else{
    return results;}
  
}

