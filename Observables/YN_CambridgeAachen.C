#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include <vector>       
#include <algorithm>       

#define s_ymax std::numeric_limits<double>::max()

using namespace ATOOLS;

namespace RESUM {

  template <unsigned int NJETS>
  class YN_CambridgeAachen: public Observable_Base {
  public:

    YN_CambridgeAachen(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<Vec4D>& p,
       const std::vector<Flavour>& fl,
       const size_t &l) {
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      return [](double Rp) {return 1.;};
    }

    double KT2(const Vec4D &p1, const Vec4D &p2) const
    {
      return 2.0*sqr(Min(p1[0],p2[0]))*(1.0-p1.CosTheta(p2));
    }

    double y12(const Vec4D &p1, const Vec4D &p2) const {
      return (1.0-p1.CosTheta(p2));
    }

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin) {
      return _Value(ip,fl,nin,s_ymax);
    }


  private:
    double _Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin,
                 const double cut)
    {
      Vec4D sum;
      size_t nn = ip.size();
      Vec4D_Vector p(&ip[nin],&ip[nn]);
      if(p.size() < NJETS) return 0;
      for (size_t i(0);i<p.size();++i) {
        sum+=p[i];
      }
      Poincare cms(sum);
      for (size_t i(0);i<p.size();++i) cms.Boost(p[i]);
      double Q2(sum.Abs2());
      std::vector<int> imap(p.size());
      for (int i=0;i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > yij
	(imap.size(),std::vector<double>(imap.size()));
      int ii=0, jj=0, n=p.size();
      double dmin=s_ymax;
      for (int i=0;i<n;++i)
	for (int j=0;j<i;++j) {
	  double dij=yij[i][j]=y12(p[i],p[j]);
	  if (dij<dmin) {
            dmin=dij;
            if(p[i][0]<p[j][0]) {
              ii=i;
              jj=j;
            }
            else {
              ii=j;
              jj=i;
            }
          }
	}
      double kt2_max = 0;
      int jets = 1;
      while (n>1) {
	if (ii!=jj) {
          double kt2 = KT2(p[imap[jj]],p[imap[ii]]);
          if(kt2 < cut) {
            p[imap[jj]]+=p[imap[ii]];
            if(kt2_max < kt2) kt2_max = kt2;
          }
          else {
            jets++;
          }
        }
	else THROW(fatal_error,"Invalid clustering");
	--n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	int jjx=imap[jj];
	for (int j=0;j<jj;++j) yij[jjx][imap[j]]=y12(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) yij[imap[i]][jjx]=y12(p[imap[i]],p[jjx]);
	ii=jj=0; dmin=s_ymax;
	for (int i=0;i<n;++i)
	  for (int j=0;j<i;++j) {
	    double dij=yij[imap[i]][imap[j]];
	    if (dij<dmin) {
              dmin=dij;
              if(p[imap[i]][0]<p[imap[j]][0]) {
                ii=i;
                jj=j;
              }
              else {
                ii=j;
                jj=i;
              }
            }
	  }
      }
      if(jets >= NJETS) {
        return cut/Q2;
      }
      else {
        return _Value(ip,fl,nin,kt2_max);
      }
    }

  };// end of class YN_CambridgeAachen

}// end of namespace RESUM

using namespace RESUM;

typedef YN_CambridgeAachen<3> Y3_CambridgeAachen;
DECLARE_GETTER(Y3_CambridgeAachen,"Y3_CambridgeAachen",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y3_CambridgeAachen>::
operator()(const Parameter_Type &args) const 
{ return new Y3_CambridgeAachen(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y3_CambridgeAachen>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y3_CambridgeAachen"; }

typedef YN_CambridgeAachen<4> Y4_CambridgeAachen;
DECLARE_GETTER(Y4_CambridgeAachen,"Y4_CambridgeAachen",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y4_CambridgeAachen>::
operator()(const Parameter_Type &args) const 
{ return new Y4_CambridgeAachen(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y4_CambridgeAachen>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y4_CambridgeAachen"; }

typedef YN_CambridgeAachen<5> Y5_CambridgeAachen;
DECLARE_GETTER(Y5_CambridgeAachen,"Y5_CambridgeAachen",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y5_CambridgeAachen>::
operator()(const Parameter_Type &args) const 
{ return new Y5_CambridgeAachen(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y5_CambridgeAachen>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y5_CambridgeAachen"; }

typedef YN_CambridgeAachen<6> Y6_CambridgeAachen;
DECLARE_GETTER(Y6_CambridgeAachen,"Y6_CambridgeAachen",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y6_CambridgeAachen>::
operator()(const Parameter_Type &args) const 
{ return new Y6_CambridgeAachen(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y6_CambridgeAachen>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y6_CambridgeAachen"; }
