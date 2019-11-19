#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"

namespace RESUM {

  class NLO_Analysis: public ANALYSIS::Analysis_Base {
  public:

    NLO_Analysis(const ANALYSIS::Argument_Matrix &params);

    ~NLO_Analysis();

    double sumW = 0.;
    int n;
    void Evaluate(double weight, double ncount,int mode) override;
    void EndEvaluation(double scale) override;
    void Output(const std::string& pname) override;
    void EvaluateNLOevt() override;
    std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                        const std::vector<ATOOLS::Flavour>& fl,
                        const size_t &nin,
                        int njets, int mode);
    double KT2(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2, const std::vector<int>& fl1, const std::vector<int>& fl2, 
               int mode, int num_flavd) const;
    bool flavd(const std::vector<int>& fls) const;
    int num_flavd(const std::vector<std::vector<int>>& flavs) const;
    void Fill(int i, double value, double weight, double ncount, int mode);
    std::map<std::string, NLO_Analysis*> m_channels;
    Primitive_Observable_Base * Copy() const;
    std::string m_addition = "";
    int m_nborn = -1;
    int m_fills = 0;
    int m_fills_tmp = 0;
    enum MODE {
               DEFAULT = 0,
               SKIP_REAL = 1        << 0,  // skip the real correction
               SKIP_SUBT = 1        << 1,  // skip all subtraction terms
               SKIP_SCOL = 1        << 2,  // skip the term related to the soft part of the splitting function
               SKIP_COLL = 1        << 3,  // skip the term related to the collinear part if the splitting function
               SKIP_NORM = 1        << 4,  // skip the term related to the normalization of the observable (dbar - b log(...))
               SKIP_SOFT = 1        << 5,  // skip the soft anomalous dimenstion matrix
               SKIP_PDFR = 1        << 6   // skip the pdf ratio
    };

    const std::map<std::string,MODE> m_ModeToEnum = {{"DEFAULT", MODE::DEFAULT},
                                                     {"SKIP_REAL", MODE::SKIP_REAL},
                                                     {"SKIP_SUBT", MODE::SKIP_SUBT},
                                                     {"SKIP_SCOL", MODE::SKIP_SCOL},
                                                     {"SKIP_COLL", MODE::SKIP_COLL},
                                                     {"SKIP_NORM", MODE::SKIP_NORM},
                                                     {"SKIP_SOFT", MODE::SKIP_SOFT},
                                                     {"SKIP_PDFR", MODE::SKIP_PDFR}};


    
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    std::map<std::string,std::vector<std::vector<double>>> m_sigma;
    
    MODE m_mode = MODE::DEFAULT;

    double EtaBeam(const ATOOLS::Vec4D& p, size_t beamId);
    

  };// end of class NLO_Analysis

}// end of namespace RESUM

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "Tools/StringTools.H"
#include <algorithm>


using namespace RESUM;
using namespace ANALYSIS;
using namespace ATOOLS;

NLO_Analysis::NLO_Analysis(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_Resum";
  Data_Reader reader(",",";","!","=");
  const std::string& mode = reader.GetValue<std::string>("MATCHING_TEST","DEFAULT");
  m_nborn = reader.GetValue<int>("RESUM::NBORN", -1);
  if(m_nborn > 0) {
    if(rpa->gen.Variable("RESUM::ANALYSIS_DO_CHANNELS")!="NO") {
      rpa->gen.SetVariable("RESUM::ANALYSIS_DO_CHANNELS","NO");
      for(int j=0; j<=m_nborn; j+=2) {
        m_channels.emplace(std::string(m_nborn-j,'g')+std::string(j,'q'),
                           dynamic_cast<NLO_Analysis*>(Copy()));
        m_channels.emplace(std::string(m_nborn-j,'g')+std::string(j,'q')+std::string("_BLAND"),
                           dynamic_cast<NLO_Analysis*>(Copy()));
        if(j>0)
          m_channels.emplace(std::string(m_nborn-j,'g')+std::string(j,'q')+std::string("_BLAND_Z"),
                             dynamic_cast<NLO_Analysis*>(Copy()));
      }
      m_channels.emplace("other",dynamic_cast<NLO_Analysis*>(Copy()));
      rpa->gen.SetVariable("RESUM::ANALYSIS_DO_CHANNELS","YES"); 
    }
  }
  if(is_int(mode)) {
    m_mode = static_cast<MODE>(to_type<int>(mode));
  }
  else {
    for(const std::string& m: split(mode,"\\|")) {
      m_mode = static_cast<MODE>(m_mode | m_ModeToEnum.at(m));
    }
  }
  Algebra_Interpreter *ip=reader.Interpreter();
  for (size_t i(1);i<params.size();++i) {
    if (params[i].size()<5) continue;
    Observable_Base *obs(Observable_Getter::GetObject
			 (params[i][0],Observable_Key(params[i][0])));
    if (obs == nullptr) {
      msg_Error()<<METHOD<<"(): Observable not found '"<<params[i][0]<<"'.\n";
      continue;
    }
    double xmin = ToType<double>(ip->Interprete(params[i][1]));
    double xmax = ToType<double>(ip->Interprete(params[i][2]));
    size_t nbin = ToType<size_t>(ip->Interprete(params[i][3]));
    int tp = HistogramType(params[i][4]);
    msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp
		   <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    m_histos.push_back(new Histogram(tp,xmin,xmax,nbin,params[i][0]));
    m_histos.push_back(new Histogram(tp,0,1,2,params[i][0]+"_Sigma"));
    for(int j=0; j<nbin; j++)
      m_histos.push_back(new Histogram(tp,0,1,2,params[i][0]+"_Sigma"));
    m_obss.push_back(obs);
  }
}

NLO_Analysis::~NLO_Analysis()
{
  for (size_t i(0);i<m_obss.size();++i) delete m_obss[i];
  for (auto& ch: m_channels) if(ch.second) delete ch.second;
}

double NLO_Analysis::EtaBeam(const Vec4D& p, size_t beamId) {
  double eta = 1./(2.0*rpa->gen.PBeam(beamId)[0]);
  if(beamId == 0) {
    eta *= p.PPlus();
  }
  else if(beamId == 1) {
    eta *= p.PMinus();
  }
  else THROW(fatal_error, "Something went wrong, beamId > 1.");
  return eta;
}


bool NLO_Analysis::flavd(const std::vector<int>& fls) const {
  for (int f: fls) {
    if(f!=0) return true;
  }
  return false;
}


double NLO_Analysis::KT2(const Vec4D &p1, const Vec4D &p2, 
                         const std::vector<int>& fl1, const std::vector<int>& fl2,
                         int mode, int num_flavd) const {
  // msg_Out()<<mode<<" "<<num_flavd;
  // msg_Out()<<" "<<fl1<<" "<<fl2<<" ";
  if(mode&2 and num_flavd < 3 and flavd(fl1) and flavd(fl2)) {
    // msg_Out()<<" "<<"lksdjflakjsf\n";
    return std::numeric_limits<double>::infinity();
  }
  // msg_Out()<<"\n";
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


std::string NLO_Analysis::Channel(const std::vector<Vec4D>& ip,
                                  const std::vector<Flavour>& fl,
                                  const size_t &nin,
                                  int njets, int mode) {
  Vec4D sum;
  size_t nn = ip.size();
  Vec4D_Vector p(&ip[nin],&ip[nn]);
  if(p.size() < njets) return "None";
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
  msg_Debugging()<<"Cluster input: mode = "<<mode<<"\n";
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
      double dij=kt2ij[i][j]=KT2(p[i],p[j],f[i],f[j],mode,num_flavd);
      // msg_Out()<<f[i]<<" "<<f[j]<<" "<<dij<<"\n";
      if (dij<dmin) { dmin=dij; ii=i; jj=j; }
    }
  // msg_Out()<<"\n\n";
  while (n>njets) {
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
    // int jjx=imap[jj];
    // for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=KT2(p[jjx],p[imap[j]],f[jjx],f[imap[j]],mode,num_flavd);
    // for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=KT2(p[imap[i]],p[jjx],f[imap[i]],f[jjx],mode,num_flavd);
    ii=jj=0; dmin=Q2;
    for (int i=0;i<n;++i)
      for (int j=0;j<i;++j) {
        double dij=kt2ij[imap[i]][imap[j]]=KT2(p[imap[i]],p[imap[j]],f[imap[i]],f[imap[j]],mode,num_flavd);
        // msg_Out()<<dij<<" "<<dmin<<" "<<f[imap[i]]<<" "<<f[imap[j]]<<"\n";
        if (dij<dmin) { dmin=dij; ii=i; jj=j; }
      }
  }
  std::string channel = "";
  msg_Debugging()<<"Clustered:\n";
  for(int i=0; i<n; i++) {
      msg_Debugging()<<p[imap[i]]<<" "<<f[imap[i]]<<"\n";
    if(flavd(f[imap[i]])) {
      if(mode&1) channel += "q";
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


void NLO_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  // sumW += weight;
  // n += ncount;
  // msg_Out()<<"weight: "<<sumW/n<<"\n";
  // calculate observable for event and fill into histo
  Particle_List all(*p_ana->GetParticleList(m_listname));
  Vec4D_Vector mom(2+all.size());
  Flavour_Vector fl(2+all.size());
  for (size_t i=0; i<all.size(); i++) {
    mom[2+i]=all[i]->Momentum();
    fl[2+i]=all[i]->Flav();
  }
  fl[0]=rpa->gen.Beam1();
  fl[1]=rpa->gen.Beam2();
  std::string ch = "";
  std::string ch_bland = "";
  std::string ch_bland_z = "";
  if(m_nborn > 0) {
    ch = Channel(mom,fl,2,m_nborn,0); 
    if(m_channels.find(ch) == m_channels.end()) {
      msg_Error()<<"Channel not found: "<<ch<<"\n";
      msg_Error()<<"Available:\n";
      for(auto& c: m_channels) msg_Error()<<c.first<<"\n";
    } 
    ch_bland = Channel(mom,fl,2,m_nborn,0|1)+"_BLAND"; 
    if(m_channels.find(ch_bland) == m_channels.end()) {
      msg_Error()<<"Channel not found: "<<ch_bland<<"\n";
      msg_Error()<<"Available:\n";
      for(auto& c: m_channels) msg_Error()<<c.first<<"\n";
    }
    ch_bland_z = Channel(mom,fl,2,m_nborn,0|1|2)+"_BLAND_Z"; 
    if(m_channels.find(ch_bland_z) == m_channels.end()) {
      msg_Error()<<"Channel not found: "<<ch_bland_z<<"\n";
      msg_Error()<<"Available:\n";
      for(auto& c: m_channels) msg_Error()<<c.first<<"\n";
    }
    
    for(auto& c: m_channels) {
      if(c.first != ch and c.first != ch_bland and c.first != ch_bland_z)  {
        c.second->AddZeroPoint(ncount,mode);
      }
    }

  }
  int obsId = 0;
  for (size_t i=0; i<m_histos.size(); i+=2+m_histos[i]->Nbin()) {
    if(mode == 0) m_fills += ncount;
    if(mode == 1) m_fills_tmp = ncount;
    const double value = m_obss[obsId]->Value(mom, fl);
    // msg_Out()<<ch<<"\n";
    // if(ch=="qqqq") exit(1);
    msg_Debugging()<<"value["<<i<<"] = "<<value
                   <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
    Fill(i,value,weight,ncount,mode);
    if(m_nborn > 0) {
      m_channels.at(ch)->Fill(i,value,weight,ncount,mode);
      m_channels.at(ch_bland)->Fill(i,value,weight,ncount,mode);
      m_channels.at(ch_bland_z)->Fill(i,value,weight,ncount,mode);
    }
    obsId++;
    // if(value > 0) exit(1);
  }
  return;
}

void NLO_Analysis::Fill(int i, double value, double weight, double ncount, int mode) {
  FillHisto(i,value,weight,ncount,mode);
  int fill = 1;
  if(value < m_histos[i]->LowEdge(0))
    fill = 0;
  FillHisto(i+1,fill,weight,ncount,mode);
  for(int j=0; j<m_histos[i]->Nbin(); j++) {
    int fill = 1;
    if(value < m_histos[i]->HighEdge(j))
      fill = 0;
    FillHisto(i+j+2,fill,weight,ncount,mode);
  }
}

void NLO_Analysis::EvaluateNLOevt() {
  m_fills += m_fills_tmp;
  m_fills_tmp = 0;
  for(auto& ch: m_channels) ch.second->EvaluateNLOevt();
  Analysis_Base::EvaluateNLOevt();  
}

void NLO_Analysis::EndEvaluation(double scale)
{
  for(auto& ch: m_channels) {
    ch.second->EndEvaluation(scale);
  }
  for(auto& h: m_histos) {
    h->MPISync();
  }
  auto histos = m_histos;
  m_histos.clear();
  for (size_t i=0; i<histos.size(); i+=2+histos[i]->Nbin()) {
    m_histos.push_back(histos[i]);
    m_sigma[histos[i]->Name()] = std::vector<std::vector<double>>(histos[i]->Nbin()+1);


    auto& h = histos[i+1];
    double v = histos[i]->LowEdge(0);
    double N = h->Fills();
    double sigW = h->Bin(0)/N;
    double sigW2 = h->Bin2(0)/N;
    double sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
    // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
    double BsigW = h->Bin(1)/N;
    double BsigW2 = h->Bin2(1)/N;
    double BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
    // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
    m_sigma[histos[i]->Name()][0] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};

    for(int j=0; j<histos[i]->Nbin(); j++) {
      auto& h = histos[i+j+2];
      double v = histos[i]->HighEdge(j);
      double N = h->Fills();
      double sigW = h->Bin(0)/N;
      double sigW2 = h->Bin2(0)/N;
      double sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
      // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
      double BsigW = h->Bin(1)/N;
      double BsigW2 = h->Bin2(1)/N;
      double BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
      // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
      m_sigma[histos[i]->Name()][j+1] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};
      delete histos[i+j+1];
    }
  }
  Analysis_Base::EndEvaluation(scale);
}

void NLO_Analysis::Output(const std::string& pname) {
// #ifdef USING__MPI
//   if (MPI::COMM_WORLD.Get_rank()) return;
// #endif
  for(auto& ch: m_channels) {
    ch.second->m_addition = "Channel_"+ch.first+"_";
    ch.second->Output(pname);
  }
  for(auto& sig: m_sigma) {
    My_Out_File ofile(pname+"/"+m_addition+sig.first+"_Sigma.dat");
    ofile.Open();
    ofile->precision(ToType<int>(rpa->gen.Variable("HISTOGRAM_OUTPUT_PRECISION")));
    *ofile<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
    for(auto& vals: sig.second) {
      for(double val: vals) *ofile<<val<<" ";
      *ofile<<"\n";
    }                         
  }
  Analysis_Base::Output(pname);
}

Primitive_Observable_Base *NLO_Analysis::Copy() const 
{
  return new NLO_Analysis(m_params);
}

DECLARE_GETTER(NLO_Analysis,"NLOResum",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLO_Analysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return nullptr;
  return new NLO_Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLO_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
