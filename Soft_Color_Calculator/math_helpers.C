#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>


unsigned int fact[] = {1, 1, 2, 6, 24, 120, 720,
                       5040, 40320, 362880, 3628800, 39916800, 479001600
                      };


void ppart(std::vector< std::vector<int> > rhs)
{
    for(unsigned i = 0; i < rhs.size(); i++)
    {
        for(unsigned j = 0; j < rhs[i].size(); j++)
        {
            std::cout << rhs[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

inline std::vector< int > cop(int n, int * a)
{
    std::vector< int > ires;
    ires.clear();
    for (unsigned i = 0; i <= n; i++)
    {
        ires.push_back(a[i]);
    }
    return ires;
}


void integerPartition(int n, int * a, int level, std::vector< std::vector<int> >& rhs)
{
    int first;
    int i;
    if (n < 1) return;
    a[level] = n;
    //print(level, a);
    rhs.push_back(cop(level, a));
    first = (level == 0) ? 2 : a[level-1];
    for(i = first; i <= n / 2; i++)
    {
        a[level] = i;
        integerPartition(n - i, a, level + 1, rhs);
    }
}


inline std::vector< std::string > perms(std::string def)
{
    std::vector< std::string > results;
    int perm=1, digits=def.size();

    for (int i=1; i<=digits; perm*=i++);
    for (int a=0; a<perm; a++)
    {
        std::string avail=def;
        std::string thin;
        for (int b=digits,div=perm; b>0; b--)
        {
            div/=b;
            int index = (a/div)%b;
            thin = thin + avail[index];
            avail.erase(index,1) ;  // lexigraphically correct
        }
        results.push_back(thin);
    }
    for(int i = 0; i< results.size(); i++)
    {
    }


    return results;
}




//Generates the disconnected color permutations (for max 2->5)
void allo_discon_groups_impl(
    std::vector<char>& usage, std::vector<int>& order,
    const std::vector<int>& cumsum, int group, int place,
    int min, int sym_factor, int sym_size,std::vector< std::string >& results,
    std::string thenums)
{
    char digitset[10];
    strcpy(digitset,thenums.c_str());

    std::string thegroups("");

    if (place == cumsum[group])
    {
        group++;
        min = 0;
        if (group == cumsum.size())
        {
            for( std::vector<int>::iterator it = order.begin(); it != order.end(); ++it )
            {
                thegroups += digitset[*it];
            }

            //for identical groups, enforce that elements are strictly increasing
            int count = 0;
            for(size_t i(0); i < sym_factor-1 ; i++)
            {
                if(thegroups[i*sym_size] < thegroups[(i+1)*sym_size]) count++;
            }
            if(count == sym_factor-1) results.push_back(thegroups);

            thegroups.clear();
            return;
        }
    }


    for( int v = min, max = usage.size() + place - cumsum[group]; v <= max; ++v )
    {
        if (!usage[v])
        {
            order[place] = v;
            usage[v] = 1;
            allo_discon_groups_impl(usage, order, cumsum, group, place+1, v+1,sym_factor,sym_size,results,digitset);
            usage[v] = 0;
        }
    }
}


inline std::vector< std::vector<int> > gluon_partitions(int ng)
{
    std::vector< std::vector<int> > results;
    int * a = (int * ) malloc(sizeof(int) * ng);
    integerPartition (ng, a, 0, results);
    return results;
}
