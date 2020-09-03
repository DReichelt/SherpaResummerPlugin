#include "Analysis/ChannelAlgorithms/KT2_pp_Ordered.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include <limits>
#include <algorithm>

using namespace RESUM;
using ATOOLS::Vec4D;
using ATOOLS::Flavour;
using ATOOLS::Poincare;
using ATOOLS::sqr;
using ATOOLS::Min;
using ATOOLS::Max;
using ATOOLS::IsZero;
using std::vector;

KT2_pp_Ordered::KT2_pp_Ordered(const ChAlg_Key& parameters, bool orderByPT)
  : ChannelAlgorithm_Base(parameters), m_orderByPT(orderByPT) {
  m_nborn = RESUM::to_type<int>(m_params[0]);
  m_orderByPT = RESUM::to_type<bool>(parameters.KwArg("OrderByPT",m_orderByPT?"1":"0"));
  if(m_params.size()>1) {
    if(m_params[1] == "" or
       m_params[1] == "''" or
       m_params[1] == "ALL") m_mode = MODE::ALL;
    else if(m_params[1] == "BLAND") m_mode = MODE::BLAND;
    else if(m_params[1] == "BLAND_Z") m_mode = MODE::BLAND_Z;
    else THROW(fatal_error,
               "Channel mode not knwon: Name = "+m_name+", mode = "+m_params[1]+".");
    m_sumNeutral = parameters.KwArg("SUMNEUTRAL","0") != "0";
    if(m_sumNeutral) msg_Debugging()<<"Will sum momenta of color-neutral particles.\n";
    else msg_Debugging()<<"Will use momenta of color-neutral particles individually.\n";
  }
  if(m_orderByPT) {
    if(m_mode==MODE::ALL) {
      m_channelNames.push_back("other");
      for(int j=0; j<=m_nborn; j+=2) {
        std::string help = std::string(m_nborn-j,'g')+std::string(j,'q');
        do {
          m_channelNames.push_back("qqTO"+help);
        } while ( std::next_permutation(help.begin(), help.end()) );
        help = std::string(m_nborn-j,'g')+std::string(j,'q');
        do {
          m_channelNames.push_back("ggTO"+help);
        } while ( std::next_permutation(help.begin(), help.end()) );
        if(m_nborn-j > 0) {
          help = std::string(m_nborn-j-1,'g')+std::string(j+1,'q');
          do {
            m_channelNames.push_back("gqTO"+help);
          } while ( std::next_permutation(help.begin(), help.end()) );
        }
      }
    }
    else if(m_mode==MODE::BLAND) {
      for(int j=0; j<=m_nborn; j+=2) {
        std::string help = std::string(m_nborn-j,'g')+std::string(j,'q');
        do {
          m_channelNames.push_back("qqTO"+help+std::string("_BLAND"));
        } while ( std::next_permutation(help.begin(), help.end()) );
        help = std::string(m_nborn-j,'g')+std::string(j,'q');
        do {
          m_channelNames.push_back("ggTO"+help+std::string("_BLAND"));
        } while ( std::next_permutation(help.begin(), help.end()) );
        if(m_nborn-j > 0) {
          help = std::string(m_nborn-j-1,'g')+std::string(j+1,'q');
          do {
            m_channelNames.push_back("gqTO"+help+std::string("_BLAND"));
          } while ( std::next_permutation(help.begin(), help.end()) );
        }
      }
    }
    else if(m_mode==MODE::BLAND_Z) {
      for(int j=0; j<=m_nborn; j+=2) {
        std::string help = std::string(m_nborn-j,'g')+std::string(j,'q');
        do {
          m_channelNames.push_back("qqTO"+help+std::string("_BLAND_Z"));
        } while ( std::next_permutation(help.begin(), help.end()) );
        if(j>1) {
          help = std::string(m_nborn-j,'g')+std::string(j,'q');
          do {
            m_channelNames.push_back("ggTO"+help+std::string("_BLAND_Z"));
          } while ( std::next_permutation(help.begin(), help.end()) );
        }
        if(m_nborn-j > 0) {
          help = std::string(m_nborn-j-1,'g')+std::string(j+1,'q');
          do {
            m_channelNames.push_back("gqTO"+help+std::string("_BLAND_Z"));
          } while ( std::next_permutation(help.begin(), help.end()) );
        }
      }
    }
  }
  else {
    if(m_mode==MODE::ALL) {
      m_channelNames.push_back("other");
      for(int j=0; j<=m_nborn; j+=2) {
        m_channelNames.push_back("qqTO"+std::string(m_nborn-j,'g')
                                 +std::string(j,'q'));
        m_channelNames.push_back("ggTO"+std::string(m_nborn-j,'g')
                                 +std::string(j,'q'));
        if(m_nborn-j > 0) {
          m_channelNames.push_back("gqTO"+std::string(m_nborn-j-1,'g')
                                   +std::string(j+1,'q'));
        }
      }
    }
    else if(m_mode==MODE::BLAND) {
      for(int j=0; j<=m_nborn; j+=2) {
        m_channelNames.push_back("qqTO"+std::string(m_nborn-j,'g')
                                 +std::string(j,'q')+std::string("_BLAND"));
        m_channelNames.push_back("ggTO"+std::string(m_nborn-j,'g')
                                 +std::string(j,'q')+std::string("_BLAND"));

        if(m_nborn-j > 0) {
          m_channelNames.push_back("gqTO"+std::string(m_nborn-j-1,'g')
                                   +std::string(j+1,'q')+std::string("_BLAND"));
        }
      }
    }
    else if(m_mode==MODE::BLAND_Z) {
      for(int j=0; j<=m_nborn; j+=2) {
        m_channelNames.push_back("qqTO"+std::string(m_nborn-j,'g')
                                 +std::string(j,'q')+std::string("_BLAND_Z"));
        if(j>1) {
          m_channelNames.push_back("ggTO"+std::string(m_nborn-j,'g')
                                   +std::string(j,'q')+std::string("_BLAND_Z"));
        }
        if(m_nborn-j > 0) {
          m_channelNames.push_back("gqTO"+std::string(m_nborn-j-1,'g')
                                   +std::string(j+1,'q')+std::string("_BLAND_Z"));
        }
      }
    }
  }
}


bool KT2_pp_Ordered::flavd(const vector<int>& fls) const {
  for (int f: fls) {
    if(f!=0) return true;
  }
  return false;
}

double KT2_pp_Ordered::KT2(const Vec4D &p1, const Vec4D &p2, 
                   const vector<int>& fl1,
                   const vector<int>& fl2,
                   int mode, int num_flavd) const {
  if(mode&MODE::ZPROD and num_flavd < 4 and flavd(fl1) and flavd(fl2)) {
    return std::numeric_limits<double>::infinity();
  }
  // determine combined flavour
  vector<int> nfl(6,0);
  if(mode&MODE::BLAND)
    for(int i=0; i<6; i++) nfl[i] = fl1[i]+fl2[i];
  // angular distance
  const double deltaR2 = sqr(p1.Y()-p2.Y())+sqr(p1.DPhi(p2));

  if(p1.PPerp2() < p2.PPerp2()) {
    // softer jet is p1
    if(flavd(fl1)) {
      // softer jet is flavoured, use max kt2
      if(mode&MODE::BLAND and flavd(nfl) and flavd(fl2))
        return std::numeric_limits<double>::infinity();
      return p2.PPerp2()*deltaR2;
    }
    else {
      // softer jet is flavourless, use min kt2
      return p1.PPerp2()*deltaR2;
    }
  }
  else {
    // softer jet is p2
    if(flavd(fl2)) {
      // softer jet is flavoured, use max kt2
      if(mode&MODE::BLAND and flavd(nfl) and flavd(fl1))
        return std::numeric_limits<double>::infinity();
      return p1.PPerp2()*deltaR2;
    }
    
    else {
      // softer jet is flavourless, use min kt2
      return p2.PPerp2()*deltaR2;
    }
  }
}

double KT2_pp_Ordered::KT2(const Vec4D &p1,
                   const double ktB2,
                   const vector<int>& fl1,
                   const vector<int>& fl2,
                   int mode, int num_flavd) const {
  if(mode&MODE::ZPROD and num_flavd < 4 and flavd(fl1) and flavd(fl2)) {
    return std::numeric_limits<double>::infinity();
  }
  if(flavd(fl1)) {
    if(m_mode&MODE::BLAND  and flavd(fl2)) {
      for(int i=0; i<6; i++) {
        if(fl1[i] != fl2[i]) {
          return std::numeric_limits<double>::infinity();
        }
      }
    }
    return Max(p1.PPerp2(),ktB2);
  }
  else {
    return Min(p1.PPerp2(),ktB2);
  }
}

double KT2_pp_Ordered::KT2(const vector<Vec4D> &ps,
                           const double eta,
                           const int beamIdx) const {
  if(beamIdx != 0 and beamIdx != 1)
    THROW(fatal_error, "Beam idx > 1 provided.");

  double ktB = 0;
  for(const Vec4D& p: ps) {
    // TODO: is this correct???
    const double deltaEta = beamIdx==0 ?  p.Y()-eta : -p.Y()-eta;
    if(IsZero(deltaEta)) ktB += sqrt(p.PPerp2()+p.Abs2())/2. * (1. + exp(deltaEta)); 
    else ktB += sqrt(p.PPerp2()+p.Abs2()) * (deltaEta > 0 ? 1. : exp(deltaEta)); 
  }
  return sqr(ktB);
}

double KT2_pp_Ordered::KT2(const vector<Vec4D> &pp,
                           const double eta,
                           const int beamIdx,
                           const vector<int> imap,
                           const size_t n,
                           const vector<Vec4D> &neutralFinal) const {
  if(beamIdx != 0 and beamIdx != 1)
    THROW(fatal_error, "Beam idx > 1 provided.");
  double ktB = 0;
  for(size_t i=0; i<n+neutralFinal.size(); i++) {
    const Vec4D& p = i<n ? pp[imap[i]] : neutralFinal[i-n];
    const double deltaEta = beamIdx==0 ?  p.Y()-eta : -(p.Y()-eta);
    if(IsZero(deltaEta)) ktB += sqrt(p.PPerp2()+p.Abs2())/2. * (1. + exp(deltaEta)); 
    else ktB += sqrt(p.PPerp2()+p.Abs2()) * (deltaEta > 0 ? 1. : exp(deltaEta)); 
  }
  return sqr(ktB);
}


std::string KT2_pp_Ordered::Channel(const vector<Vec4D>& ip,
                                    const vector<Flavour>& fl,
                                    const size_t &nin,
                                    vector<Vec4D>* pout=nullptr) {
  if(!fl[0].Strong() or !fl[1].Strong())
    THROW(fatal_error,"This algorithm is only valid for pp collisions, but one beam had no strong charge.");
  
  size_t nn = ip.size();
  msg_Debugging()<<"Got "<<nn<<" momenta and "<<fl.size()<<" flavours. "
                 <<"From these, "<<nin<<" are initial states.\n";
  // TODO: was this doing something useful?
  // if(nn == nin) {
  //   std::string channel = "";
  //   if(fl[0] == 21) channel += "g";
  //   if(fl[1] == 21) channel += "g";
  //   if(fl[0] != 21) channel += "q";
  //   if(fl[1] != 21) channel += "q";
  // }
  Vec4D sum;
  vector<Vec4D> p;
  vector<Vec4D> neutralFinal;
  vector<vector<int>> f;//(nn,vector<int>(6,0));
  for (size_t i=0; i<nn; i++) {
    msg_Debugging()<<i<<"th particle: "<<fl[i]<<", p = "<<ip[i];
    if(fl[i].Strong()) {
      msg_Debugging()<<" is colour-charged ";
      vector<int>  fi(6,0);
      if(fl[i] == 1) fi[0] += 1;
      if(fl[i] == -1) fi[0] -= 1;
      if(fl[i] == 2) fi[1] += 1;
      if(fl[i] == -2) fi[1] -= 1;
      if(fl[i] == 3) fi[2] += 1;
      if(fl[i] == -3) fi[2] -= 1;
      if(fl[i] == 4) fi[3] += 1;
      if(fl[i] == -4) fi[3] -= 1;
      if(fl[i] == 5) fi[4] += 1;
      if(fl[i] == -5) fi[4] -= 1;
      if(fl[i] == 6) fi[5] += 1;
      if(fl[i] == -6) fi[5] -= 1;
      sum+=p[i];
      f.push_back(fi);
      if(i < nin) {
        msg_Debugging()<<" but not a final state.\n";
      }
      else {
        msg_Debugging()<<" and a final state.\n";
        p.push_back(ip[i]);
      }
    }
    else {
      msg_Debugging()<<" is not colour charged ";
      if(i < nin) {
        msg_Debugging()<<" but not a final state.\n";
      }
      else {
        msg_Debugging()<<" and a final state.\n";
        if(!m_sumNeutral or neutralFinal.size()==0) neutralFinal.push_back(ip[i]);
        else neutralFinal[0] += ip[i];
      }
    }
  }
  int  n = p.size();
  msg_Debugging()<<"Found "<<n<<" colour charged final states to use as cluster input.\n";
  int num_flavd = 0;
  msg_Debugging()<<"Cluster input: mode = "<<m_mode<<"\n";
  for(int i=0; i<f.size(); i++) {
    if(flavd(f[i])) num_flavd++;
    if(i<2) {
      msg_Debugging()<<"("<<i<<") "<<f[i]<<"\n";
    }
    else {
      msg_Debugging()<<"("<<i<<") "<<f[i]<<" "<<p[i-2]<<"\n";
    }
  }
  if(m_mode&MODE::ZPROD and num_flavd < 2) {
    THROW(fatal_error, "Clustering to states compatible with Z-production was requested, but there are no flavours in input.");
  }
  
  vector<int> imap(p.size());
  for (int i=0;i<imap.size();++i) imap[i]=i;
  vector<vector<double> > kt2ij
    (imap.size(),vector<double>(imap.size()));
  int ii=0;
  int jj=0;
  // min distance between final states
  double dmin=std::numeric_limits<double>::infinity();
  // min distance to first beam
  double dBmin=std::numeric_limits<double>::infinity();
  // min distance to second beam
  double dBBmin=std::numeric_limits<double>::infinity(); 
  for (int i=0;i<n;++i) {
    // distance of i to first beam
    const double dB = KT2(p[i],KT2(p,p[i].Y(),0,imap,n,neutralFinal),
                          f[i+2],f[0],m_mode,num_flavd);
    msg_Debugging()<<i+2<<" with \n0: "<<dB<<"\n";
    // distance of i to second beam
    const double dBB = KT2(p[i],KT2(p,p[i].Y(),1,imap,n,neutralFinal),
                           f[i+2],f[1],m_mode,num_flavd);
    msg_Debugging()<<"1: "<<dBB<<"\n";
    const double di=Min(dB,dBB);
    if(di<dmin) { dmin=di; ii=i; jj=i; dBmin=dB; dBBmin=dBB; }
    for (int j=0;j<i;++j) {
      // distance between i and j
      double dij=kt2ij[i][j]=KT2(p[i],p[j],f[i+2],f[j+2],m_mode,num_flavd);
      msg_Debugging()<<j+2<<": "<<dij<<"\n";
      if (dij<dmin) { dmin=dij; ii=i; jj=j; }
    }
  }
  msg_Debugging()<<"Starting with "<<n<<" final states to cluster.\n";
  while (n>m_nborn) {
    msg_Debugging()<<"Need clustering because "<<m_nborn<<" < "<<n<<"\n";
    msg_Debugging()<<"Cluster "<<imap[ii]+2<<" to "<<imap[jj]+2<<"\n";
    if (ii!=jj) {
      // combine particles imap[ii] and imap[jj]
      p[imap[jj]]+=p[imap[ii]];
      for(int i=0; i<6; i++) f[imap[jj]+2][i] += f[imap[ii]+2][i];
    }
    else {
      if(dBmin < dBBmin) {
        if(dBmin!=dmin) THROW(fatal_error,"Something went wrong.");
        // combine particle with first beam
        for(int i=0; i<6; i++) f[0][i] -= f[imap[ii]+2][i];
      }
      else {
        if(dBBmin!=dmin) THROW(fatal_error,"Something went wrong.");
        // combine particle with second beam
        for(int i=0; i<6; i++) f[1][i] -= f[imap[ii]+2][i];
      }
    }
    --n;
    // remove ii from list
    for (int i=ii;i<n;++i) imap[i]=imap[i+1];
    // determine how many flavoured objects are left
    num_flavd = 0;
    if(flavd(f[0])) num_flavd++;
    if(flavd(f[1])) num_flavd++;
    for(int i=0; i<n; i++) if(flavd(f[imap[i]+2])) num_flavd++;
    
    ii=jj=0; dmin=std::numeric_limits<double>::infinity();
    // repeat distance calculation
    for (int i=0;i<n;++i) {
      const double dB = KT2(p[imap[i]],KT2(p,p[imap[i]].Y(),0,imap,n,neutralFinal),
                            f[imap[i]+2],f[0],m_mode,num_flavd);
      const double dBB = KT2(p[imap[i]],KT2(p,p[imap[i]].Y(),1,imap,n,neutralFinal),
                             f[imap[i]+2],f[1],m_mode,num_flavd);
      const double di=Min(dB,dBB);
      if(di<dmin) { dmin=di; ii=i; jj=i; dBmin=dB; dBBmin=dBB; }
      for (int j=0;j<i;++j) {
        double dij=kt2ij[imap[i]][imap[j]]=KT2(p[imap[i]],p[imap[j]],
                                               f[imap[i]+2],f[imap[j]+2],
                                               m_mode,num_flavd);
        if (dij<dmin) { dmin=dij; ii=i; jj=j; }
      }
    }
  }
  if(n != m_nborn) {
    msg_Error()<<"Clustered to "<<n<<" final states, requested were "<<m_nborn<<".\n";
    THROW(fatal_error, "Something went wrong, did not cluster to requested multiplicity.")
  }
  else {
    msg_Debugging()<<"Finished clustering to "<<n<<" = "<<m_nborn<<" final states.\n";
  }
  std::string channel = "";
  msg_Debugging()<<"Clustered:\n";
  msg_Debugging()<<"(0)"<<f[0]<<"\n";
  msg_Debugging()<<"(1)"<<f[1]<<"\n";
  for(int i=0; i<n; i++) {
    msg_Debugging()<<"("<<i+2<<") "<<f[imap[i]+2]<<" "<<p[imap[i]]<<"\n";
  }
  if(!flavd(f[0])) {
    msg_Debugging()<<"Beam 1 not flavoured.\n";
    channel += "g";
  }
  if(!flavd(f[1])) {
    msg_Debugging()<<"Beam 2 not flavoured.\n";
    channel += "g";
  }
  if(flavd(f[0])) {
    msg_Debugging()<<"Beam 1 flavoured.\n";
    if(m_mode&MODE::BLAND) channel += "q";
    else {
        bool found = false;
        for(int fla: f[0]) {
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
        if(channel != "other") {
          channel += "q";
        }
    }
  }
  if(flavd(f[1])) {
    msg_Debugging()<<"Beam 2 flavoured.\n";
    if(m_mode&MODE::BLAND) channel += "q";
    else {
        bool found = false;
        for(int fla: f[1]) {
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
        if(channel != "other") {
          channel += "q";
        }
    }
  }
  if(channel == "other") {
    msg_Debugging()<<"Beams not single flavoured.\n";
    if(pout) {
      const Vec4D p0 = ip[0];
      const Vec4D p1 = ip[1];
      pout->clear();
      pout->push_back(p0);    
      pout->push_back(p1);    
      msg_Debugging()<<"Added beam to pout, pout = "<<*pout<<".\n";
      for(int i=0; i<n; i++) {
        pout->push_back(p[imap[i]]);
      }
      msg_Debugging()<<"pout = "<<*pout<<"\n";
    }
    return channel;
  }
  else {
    channel += "TO";
    msg_Debugging()<<"Channel after beam analysis -> "<<channel<<"\n";
  }

  std::string final_channel = "";
  if(pout) {
    const Vec4D p0 = ip[0];
    const Vec4D p1 = ip[1];
    pout->clear();
    pout->push_back(p0);    
    pout->push_back(p1);    
    msg_Debugging()<<"Added beam to pout, pout = "<<*pout<<".\n";
  }
  while(n>0) {
    int ind = -1;
    if(m_orderByPT) {
      double pTmax = -1;
      for(int i=0; i<n; i++) {
        if(p[imap[i]].PPerp() > pTmax) {
          ind = i;
          pTmax = p[imap[i]].PPerp();
        }
      }
    }
    else {
      ind = m_nborn-n;
    }
    if(pout) {
      pout->push_back(p[imap[ind]]);
      msg_Debugging()<<"pout = "<<*pout<<".\n";
    }
    msg_Debugging()<<p[imap[ind]]<<" "<<f[imap[ind]+2]<<"\n";
    if(flavd(f[imap[ind]+2])) {
      if(m_mode&MODE::BLAND) final_channel += "q";
      else {
        bool found = false;
        for(int fla: f[imap[ind]+2]) {
          // msg_Out()<<fla<<" ";
          if(fla != 0) {
            if(found or abs(fla) > 1) {
              // msg_Out()<<"other\n";
              final_channel = "other";
              break;
            }
            else found = true;
          }
        }
        // msg_Out()<<"\n";
        if(final_channel != "other") {
          // msg_Out()<<"not other\n";
          final_channel += "q";
        }
        else break;
      }
    }
    else {
      if(m_orderByPT) final_channel += "g";
      else final_channel = "g"+final_channel;
    }
    --n;
    for (int i=ind;i<n;++i) imap[i]=imap[i+1];
  }
  if(final_channel == "other") channel = "other";
  else channel += final_channel;
  if(m_mode&BLAND) channel += "_BLAND";
  if(m_mode&ZPROD) channel += "_Z";
  msg_Debugging()<<"Channel = "<<channel<<"\n\n";
  return channel;
}


DECLARE_GETTER(KT2_pp_Ordered,"KT2_pp_byPT",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_pp_Ordered>::
operator()(const Parameter_Type &args) const 
{ return new KT2_pp_Ordered(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_pp_Ordered>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"KT2_pp_byPT"; }

class KT2_pp_Unordered : public KT2_pp_Ordered {
public:
  KT2_pp_Unordered(const ChAlg_Key& parameters) : KT2_pp_Ordered(parameters, false) {}
};

DECLARE_GETTER(KT2_pp_Unordered,"KT2_pp",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_pp_Unordered>::
operator()(const Parameter_Type &args) const 
{ return new KT2_pp_Unordered(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,KT2_pp_Unordered>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"KT2_pp"; }

