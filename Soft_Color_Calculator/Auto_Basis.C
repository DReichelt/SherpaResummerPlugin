#include "Auto_Basis.H"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include "math_helpers.C"

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
//


//------------------------------------------------------------------------
//Inlines
//------------------------------------------------------------------------

inline std::vector< CBasis <int> > combineA(std::vector< CBasis <int> > rhs1 , std::vector< CBasis <int> > rhs2 )
{
    std::vector< CBasis <int> > results;
    CBasis <int> plan(0,0);

    for(unsigned i = 0; i < rhs1.size(); i++)
    {
        for(unsigned j = 0; j < rhs2.size(); j++)
        {
            plan = rhs1[i].cadd(rhs2[j]);
            results.push_back(plan);
        }
    }
    return results;
}


//Fold 2 abasis elements into 1
inline abasis fold(abasis rhs1, abasis rhs2)
{
    abasis results;
    if(rhs1.size() == 0) return rhs2;
    if(rhs2.size() == 0) return rhs1;

    int ind_hold = 0;
    results.resize(rhs1.size()*rhs2.size());

    for(unsigned i1 = 0; i1 < rhs1.size(); i1++)
    {
        for(unsigned i2 = 0; i2 < rhs2.size(); i2++)
        {
            results[ind_hold]=combineA(rhs1[i1],rhs2[i2]);
            ind_hold++;
        }
    }
    return results;
}


//Returns an adjoint delta function
inline CBasis<int>   AdjointDelta(int a, int b)
{
    CBasis<int> results(2,0);
    results.Tadd(a,a,b);
    results.Tadd(b,b,a);
    return results;
}


//Get symmetry factor of a partition
inline void  get_partition_symmetry_factor(std::vector< int > partition, size_t &sym_factor, size_t &sym_s)
{
    int last = 0;
    for(unsigned it = 0; it < partition.size(); ++it )
    {
        if(last == partition[it])
        {
            sym_factor++;
            sym_s = partition[it];
        }
        last = partition[it];
    }
}


//get groupings
inline std::vector< std::string > arrange_groups( std::vector<int> parts, std::string labels, bool is_quarks)
{
    std::vector< std::string > results;
    size_t sum_of_c = 0, last , sym_s = 0;
    size_t sym_fact=1;
    std::vector<int> cumsum_of_c;

    //get all index groupings
    if(!is_quarks) get_partition_symmetry_factor( parts , sym_fact, sym_s );

    //determine the symmetry factor for the partition
    for(unsigned it = 0; it < parts.size(); ++it )
    {
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

CBasis<int> ToAdjoint(std::string rhs)
{
    CBasis <int>  results(1,0);
    int int_start=501;
    int nsize = rhs.size();
    results.Fadd(rhs[0]-'0',rhs[1]-'0',int_start);
    results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());

    //loop over internal f's (n > 4)
    for(unsigned m = 0; m < nsize-4; m++)
    {
        results.Fadd(int_start,rhs[m+2]-'0',int_start+1);
        results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());
        int_start++;
    }

    results.Fadd(int_start,rhs[nsize-2]-'0',rhs[nsize-1]-'0');
    results.app(2/sqrt(2.)*results.get_constant(),results.get_pownc());

    return results;
}



std::vector < CBasis<int> > ToFundamental(std::string rhs, int nquark_pairs, std::string qs)
{
    std::vector< CBasis <int> > results;
    int ssize = rhs.size();
    int index;
    int sign = pow(-1,ssize); //Relative sign in Reflection term
    CBasis<int> flow(2*sqrt(2),0),flowR(2*sign*sqrt(2),0);
    flow.Tin.resize(ssize);
    flowR.Tin.resize(ssize);

    //Gluon indices
    for(unsigned i = 0; i<ssize; i++)
    {
        flow.Tin[i].resize(3);
        flow.Tin[i][1] = rhs[i]-'0';
        flow.Tin[i][0] = rhs[i]-'0';
        flow.Tin[(ssize-1+i)%ssize].resize(3);
        flow.Tin[(ssize-1+i)%ssize][2] = (rhs[i]-'0');
        index = (i == 0) ? 0 : ssize-i;
        flowR.Tin[index].resize(3);
        flowR.Tin[index][1] =  rhs[i]-'0';
        flowR.Tin[index][0] = rhs[i]-'0';
        flowR.Tin[(ssize-1+index)%ssize].resize(3);
        flowR.Tin[(ssize-1+index)%ssize][2] = (rhs[i]-'0');
    }

    //Add quark indices
    if(nquark_pairs == 1)
    {
        int qstart = 50;
        flow.Tin[0][1] = qstart+(qs[0]-'0');
        flow.Tin[ssize-1][2] = qstart+(qs[1]-'0');
    }

    if(nquark_pairs == 2 && ssize < 3)
    {
        int qstart = 50;
        flow.Tin[0][1] = qstart+(qs[0]-'0');
        flow.Tin[ssize-1][2] = qstart+(qs[2]-'0');
        flow.Tin[0][2] = qstart+(qs[3]-'0');
        flow.Tin[ssize-1][1] = qstart+(qs[1]-'0');
    }

    if(nquark_pairs == 2 && ssize >= 3)
    {
        int qstart = 50;
        flow.Tin[0][1] = qstart+(qs[0]-'0');
        flow.Tin[ssize-1][2] = qstart+(qs[2]-'0');
        flow.Tin[1][2] = qstart+(qs[3]-'0');
        flow.Tin[ssize-1][1] = qstart+(qs[1]-'0');
    }
    results.push_back(flow);
    if(nquark_pairs == 0) results.push_back(flowR);
    return results;
}



//Returns f-basis elements (pure gluons)
abasis AdjointConGluons(std::string rhs)
{
    abasis results;
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
    int perm_size = fact[nsize-2];
    for(unsigned j = 0; j < perm_size; j++)
    {
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
abasis TraceConGluons( std::string rhs, int num_quark_pairs, std::string ql)
{
    abasis results;
    std::vector< CBasis<int> > presults;
    std::vector< std::string > Perms;
    std::string base, baseT;
    int nsize = rhs.size();
    baseT = rhs;

    if(rhs.size() == 0 && num_quark_pairs>0)
    {
        CBasis<int> qdis(1,0);
        int qstart = 50;
        presults.push_back(qdis.Dadd(qstart+(ql[0]-'0'),qstart+(ql[1]-'0') ) );
        qstart++;
        results.push_back(presults);
        return results;
    }

    //trace-basis element (Ng >= 4)
    //permute (n-1) string
    if(num_quark_pairs == 0) base = baseT.substr(1,nsize-1);
    else
    {
        base = baseT;
    }
    Perms = perms(base);

    //sum over permutations
    int perm_size = fact[base.size()];
    for(unsigned j = 0; j < perm_size; j++)
    {
        base = Perms[j];
        if(num_quark_pairs == 0) baseT = baseT[0]+base;
        else
        {
            baseT = base;
        }

        //Cbasis elements with indices contracted
        //conditional assures that only terms which
        //are not reflection related are included
        //in the sum.
        if(baseT[1] < baseT[baseT.size()-1] && num_quark_pairs == 0)
        {
            presults = ToFundamental(baseT,num_quark_pairs,"");
            results.push_back(presults);
        }
        else if(num_quark_pairs > 0)
        {

            std::vector< std::string > PermsQ;
            PermsQ = perms( ql.substr(0,num_quark_pairs ) );


            for(size_t j(0); j<PermsQ.size(); j++)
            {
                presults = ToFundamental(baseT,num_quark_pairs,PermsQ[j]+ql.substr(num_quark_pairs,num_quark_pairs));
                results.push_back(presults);
            }
        }

        presults.clear();
    }

    return results;

}

abasis gluon_groups( int ngluon, std::string gluon_labels, const bool &multip )
{
    //Pure gluon component of basis
    //returns set of basis elements permuted and
    //distributed over integer partitions
    //if multip is set then the particle labels are
    //multi-perhepherally distributed
    //and basis is constructed with adjoint tensors

    abasis results, iresults, iresults2, iresultsT;
    std::vector< CBasis<int> > presults;
    size_t place;
    //get partition of gluons integers
    Int_mat numglu = gluon_partitions(ngluon);

    //for each partition get elements
    for(unsigned g = 0; g < numglu.size(); g++)
    {
        int Ngroups = numglu[g].size();
        std::string base, baseT;

        //get all index groupings
        std::vector< std::string > strt;
        strt = arrange_groups(numglu[g] , gluon_labels , false);

        //For each grouping get the corresponding elements
        for(unsigned i = 0; i < strt.size(); i++)
        {
            place = 0;

            //sum over groups in partition
            for(unsigned k = 0; k < Ngroups; k++)
            {
                int nsize = numglu[g][k];
                baseT = strt[i].substr(place,nsize);

                //adjoint delta function for ng = 2 partition
                if(nsize == 2)
                {
                    presults.push_back(AdjointDelta(baseT[0]-'0',baseT[1]-'0'));
                    iresults.push_back(presults);
                }
                else if(nsize == 3)
                {
                    iresults =  TraceConGluons(baseT,0,"");
                }
                else
                {
                    //Trace basis connected gluons
                    if(!multip) iresults = TraceConGluons(baseT,0,"");
                    //adjoint basis gluons (for pure gluons only)
                    if(multip) iresults = AdjointConGluons(baseT);
                }

                //Fold together the partial elements
                if(k >= 1)
                {
                    iresultsT = fold(iresults,iresults2);
                    iresults2 = iresultsT;
                }
                //the first or only element
                else
                {
                    iresults2 = iresults;
                }

                //book-keeping
                place += numglu[g][k];
                presults.clear();
                iresults.clear();
            }

            //insert into the basis
            results.insert(results.end(), iresults2.begin(), iresults2.end());
        }
    }

    return results;
}


abasis quark_groups(int nquark, std::string quark_labels)
{
    abasis results, iresults;
    std::vector< CBasis<int> > presults;
    std::vector< std::string > PermsQ;
    int qstart = 50;
    std::string aquark_labels = quark_labels.substr(nquark,nquark);
    PermsQ = perms( quark_labels.substr(0,nquark) );
    int keep = 0;
    int nq;
    int f1,af1,adj;

    for(unsigned k = 0; k < PermsQ.size(); k++)
    {
        CBasis<int> quark_elms(1,0);
        keep = 0;
        nq = nquark;

        for(unsigned j = 0; j < PermsQ[k].size(); j++)
        {
            f1  =  PermsQ[k][j]-'0';
            af1 =  quark_labels[j] -'0';
            quark_elms.Dadd(qstart+PermsQ[k][j]-'0', qstart+aquark_labels[j] -'0');
        }

        presults.push_back(quark_elms);
        iresults.push_back(presults);
        presults.clear();
    }
    results.insert(results.end(), iresults.begin(), iresults.end());
    return results;
}


abasis connected_qqng_groups(int ngluon,int quark_pairs,
                             int con_qp, int con_gl ,std::string gluon_labels,std::string quark_labels)
{
    abasis results, iresults, iresults_L, iresults_R, iresults_Q;

    //get the gluon groups
    std::vector<int> glu_part;
    if(ngluon-con_gl > 0 ) glu_part.push_back(ngluon-con_gl);
    glu_part.push_back(con_gl);
    std::vector< std::string > strt;
    strt = arrange_groups(glu_part,gluon_labels,true);

    for(unsigned i = 0; i < strt.size(); i++)
    {
        //disconnected gluons
        iresults_L = gluon_groups(ngluon-con_gl, strt[i].substr(0,ngluon-con_gl), false);

        //disconnected quarks
        std::string quarks = quark_labels.substr(0,quark_pairs);
        std::string aquarks = quark_labels.substr(quark_pairs,2*quark_pairs);
        std::vector<int> q_part;
        if(quark_pairs-con_qp > 0 ) q_part.push_back(quark_pairs-con_qp);
        q_part.push_back(con_qp);
        std::vector< std::string > strt_q,strt_aq;
        strt_q = arrange_groups(q_part,aquarks,true);
        strt_aq = arrange_groups(q_part,quarks,true);

        for(unsigned j = 0; j < strt_q.size(); j++)
        {

            for(unsigned k = 0; k < strt_aq.size(); k++)
            {

                std::string disQ = strt_aq[j].substr(0,(quark_pairs-con_qp)) + strt_q[k].substr(0,(quark_pairs-con_qp));
                std::string conQ = strt_aq[j].substr((quark_pairs-con_qp),quark_pairs)
                                   + strt_q[k].substr((quark_pairs-con_qp),quark_pairs);

                iresults_Q = quark_groups(quark_pairs-con_qp, disQ);



                //connected
                iresults_R = TraceConGluons( strt[i].substr(ngluon-con_gl,con_gl),con_qp,conQ);

                iresults = fold(fold(iresults_Q,iresults_R),iresults_L);

                results.insert( results.end(), iresults.begin(), iresults.end() );
            }
        }
    }

    return results;
}

abasis qqng_groups(int ngluon, int quark_pairs,
                   std::string gluon_labels, std::string quark_labels, const bool &adjoint)
{
    abasis results, iresults;

    //number of quarks color connected to gluons
    for(size_t i(1); i < quark_pairs+1; i++)
    {
        //get the quark labels

        //get the gluon labels

        iresults.clear();

        for(size_t j(1); j < ngluon+1; j++)
        {
            if(ngluon - j != 1 && j >= i)
            {

                // i quark pairs
                // connected to
                // j gluons
                std::cout << i << " quark pair connected to " << j << " gluons" << std::endl;

                iresults = connected_qqng_groups(ngluon,quark_pairs,i,j,gluon_labels,quark_labels);

                results.insert( results.end(), iresults.begin(), iresults.end() );
                iresults.clear();
            }
        }
    }

    //disconnected quark string (quark_pairs - i)
    //second conditional eliminates singlet gluon
    if(quark_pairs > 0 && ngluon != 1) iresults = quark_groups(quark_pairs,quark_labels);

    //disconnected gluon string
    if( ngluon > 1) iresults = fold( iresults , gluon_groups(ngluon, gluon_labels, adjoint) );

    results.insert( results.end(), iresults.begin(), iresults.end() );

    cout << results.size() << endl;

    return results;

}


abasis auto_Cbasis( int number_of_gluons, int number_of_quarks_pairs, const bool &adjoint )
{

    std::string gluon_labels = "123456789"; //up to 9 gluons (don't try more)
    std::string quark_labels = "12345678";  //up to 4 qqb pairs

    //Multi-peripheral gluon labelling
    if(adjoint) gluon_labels = "1"+gluon_labels.substr(2,number_of_gluons-2)+"2";
    return qqng_groups(number_of_gluons, number_of_quarks_pairs, gluon_labels, quark_labels, adjoint);

}


void printBasis(abasis rhs)
{
    int bdim = rhs.size();
    for(unsigned i = 0; i < bdim; i++)
    {
        for(unsigned j = 0; j < rhs[i].size(); j++)
        {
            cout << "c(" << i << ":" << j << ")" << endl;
            rhs[i][j].printBasInfo();
        }
        std::cout << "=====================================" << std::endl;
    }

}

