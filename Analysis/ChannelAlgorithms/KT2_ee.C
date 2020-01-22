#include "Analysis/ChannelAlgorithms/KT2_ee.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include <limits>

using namespace RESUM;
using ATOOLS::Vec4D;
using ATOOLS::Flavour;
using ATOOLS::Poincare;
using ATOOLS::sqr;

KT2_ee::KT2_ee(const ChAlg_Key& parameters)
  : ChannelAlgorithm_Base(parameters) {
  m_nborn = RESUM::to_type<int>(m_params[0]);
  if(m_params.size()>1) {
    if(m_params[1] == "") m_mode = 0;
    else if(m_params[1] == "BLAND") m_mode = 1;
    else if(m_params[1] == "BLAND_Z") m_mode = 2;
    else THROW(fatal_error,
               "Channel mode not knwon: Name = "+m_name+", mode = "+m_params[1]+".");
  }
  if(m_mode==0) m_channelNames.push_back("other");
  for(int j=0; j<=m_nborn; j+=2) {
    if(m_mode==0)
      m_channelNames.push_back(std::string(m_nborn-j,'g')+std::string(j,'q'));
    else if(m_mode==1)
      m_channelNames.push_back(std::string(m_nborn-j,'g')+std::string(j,'q')+std::string("_BLAND"));
    else if(m_mode==2)
      m_channelNames.push_back(std::string(m_nborn-j,'g')+std::string(j,'q')+std::string("_BLAND_Z"));
  }
}


bool KT2_ee::flavd(const std::vector<int>& fls) const {
  for (int f: fls) {
    if(f!=0) return true;
  }
  return false;
}

double KT2_ee::KT2(const Vec4D &p1, const Vec4D &p2, 
                         const std::vector<int>& fl1, const std::vector<int>& fl2,
                         int mode, int num_flavd) const {
  if(mode&2 and num_flavd < 3 and flavd(fl1) and flavd(fl2)) {
    return std::numeric_limits<double>::infinity();
  }
  std::vector<int> nfl(6,0);
  for(int i=0; i<6; i++) nfl[i] = fl1[i]+fl2[i];
  if(p1[0] < p2[0]) {
    if(flavd(fl1)) {
      if(mode&1 and flavd(nfl) and flavd(fl2)) return std::numeric_limits<double>::infinity();
      return 2.0*sqr(p2[0])*(1.0-p1.CosTheta(p2));
    }
    else return 2.0*sqr(p1[0])*(1.0-p1.CosTheta(p2));
  }
  else {
    if(flavd(fl2)) {
      if(mode&1 and flavd(nfl) and flavd(fl1)) return std::numeric_limits<double>::infinity();
      return 2.0*sqr(p1[0])*(1.0-p1.CosTheta(p2));
    }
    else return 2.0*sqr(p2[0])*(1.0-p1.CosTheta(p2));
  }
}

std::string KT2_ee::Channel(const std::vector<Vec4D>& ip,
                            const std::vector<Flavour>& fl,
                            const size_t &nin) {
  Vec4D sum;
  size_t nn = ip.size();
  std::vector<Vec4D> p(&ip[nin],&ip[nn]);
  std::vector<std::vector<int>> f(p.size(),std::vector<int>(6,0));
  for (size_t i(0);i<p.size();++i) {
    if(fl[i+nin] == 1) f[i][0] += 1;
    if(fl[i+nin] == -1) f[i][0] -= 1;
    if(fl[i+nin] == 2) f[i][1] += 1;
    if(fl[i+nin] == -2) f[i][1] -= 1;
    if(fl[i+nin] == 3) f[i][2] += 1;
    if(fl[i+nin] == -3) f[i][2] -= 1;
    if(fl[i+nin] == 4) f[i][3] += 1;
    if(fl[i+nin] == -4) f[i][3] -= 1;
    if(fl[i+nin] == 5) f[i][4] += 1;
    if(fl[i+nin] == -5) f[i][4] -= 1;
    if(fl[i+nin] == 6) f[i][5] += 1;
    if(fl[i+nin] == -6) f[i][5] -= 1;
    sum+=p[i];
  }
  int num_flavd = 0;
  msg_Debugging()<<"Cluster input: mode = "<<m_mode<<"\n";
  for(int i=0; i<f.size(); i++) {
    if(flavd(f[i])) num_flavd++;
    msg_Debugging()<<p[i]<<" "<<f[i]<<"\n";
  }
  Poincare cms(sum);
  for (size_t i(0);i<p.size();++i) cms.Boost(p[i]);
  double Q2(sum.Abs2());
  std::vector<int> imap(p.size());
  for (int i=0;i<imap.size();++i) imap[i]=i;
  std::vector<std::vector<double> > kt2ij
    (imap.size(),std::vector<double>(imap.size()));
  int ii=0, jj=0, n=p.size();
  double dmin=Q2;
  for (int i=0;i<n;++i)
    for (int j=0;j<i;++j) {
      double dij=kt2ij[i][j]=KT2(p[i],p[j],f[i],f[j],m_mode,num_flavd);
      if (dij<dmin) { dmin=dij; ii=i; jj=j; }
    }
  while (n>m_nborn) {
    if (ii!=jj) {
      p[imap[jj]]+=p[imap[ii]];
      for(int i=0; i<6; i++) f[imap[jj]][i] += f[imap[ii]][i];
    }
    else {
      msg_Error()<<"\n\n\nSomething went wrong clustering the following: \n";
      msg_Error()<<"ii = "<<ii<<", jj = "<<jj<<"\n";
      for(size_t i=0; i<ip.size(); i++) msg_Error()<<ip[i]<<" "<<fl[i]<<"\n";
      THROW(fatal_error,"Invalid clustering");
    }
    --n;
    for (int i=ii;i<n;++i) imap[i]=imap[i+1];
    num_flavd = 0;
    for(int i=0; i<n; i++) if(flavd(f[imap[i]])) num_flavd++;
    ii=jj=0; dmin=Q2;
    for (int i=0;i<n;++i)
      for (int j=0;j<i;++j) {
        double dij=kt2ij[imap[i]][imap[j]]=KT2(p[imap[i]],p[imap[j]],f[imap[i]],f[imap[j]],m_mode,num_flavd);
        if (dij<dmin) { dmin=dij; ii=i; jj=j; }
      }
  }
  std::string channel = "";
  msg_Debugging()<<"Clustered:\n";
  for(int i=0; i<n; i++) {
    msg_Debugging()<<p[imap[i]]<<" "<<f[imap[i]]<<"\n";
    if(flavd(f[imap[i]])) {
      if(m_mode&1) channel += "q";
      else {
        bool found = false;
        for(int fla: f[imap[i]]) {
          // msg_Out()<<fla<<" ";
          if(fla != 0) {
            if(found or abs(fla) > 1) {
              // msg_Out()<<"other\n";
              channel = "other";
              break;
            }
            else found = true;
          }
        }
        // msg_Out()<<"\n";
        if(channel != "other") {
          // msg_Out()<<"not other\n";
          channel += "q";
        }
        else break;
      }
    }
    else channel = "g"+channel;
  }
  msg_Debugging()<<"Channel = "<<channel<<"\n\n";
  return channel;
}


DECLARE_GETTER(KT2_ee,"KT2_ee",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_ee>::
operator()(const Parameter_Type &args) const 
{ return new KT2_ee(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_ee>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"KT2_ee"; }

