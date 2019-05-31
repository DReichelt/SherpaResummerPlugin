#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"

namespace RESUM {

  class NLO_Analysis: public ANALYSIS::Analysis_Base {
  public:

    NLO_Analysis(const ANALYSIS::Argument_Matrix &params);

    ~NLO_Analysis();

    
    void Evaluate(double weight, double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

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
#include "Tools/StringTools.H"
#include <algorithm>

#ifdef USING__ROOT
#undef USING__ROOT
#include "TH2D.h"
#include "AddOns/Analysis/Tools/My_Root.H"
#endif

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
    m_obss.push_back(obs);
  }
  
#ifdef USING__ROOT
  if (MYROOT::myroot==NULL) {
    MYROOT::myroot = new MYROOT::My_Root();
    MYROOT::myroot->InitFile();
    ATOOLS::exh->AddTerminatorObject(MYROOT::myroot);
  }
  (*MYROOT::myroot)(new TH2D("yzplot","yzplot",100,-16.,0.,100,-16.,0.),"yzplot");
#endif

}

NLO_Analysis::~NLO_Analysis()
{
#ifdef USING__ROOT
  if (MYROOT::myroot) {
    exh->RemoveTerminatorObject(MYROOT::myroot);
    delete MYROOT::myroot;
    MYROOT::myroot=NULL;
  }
#endif
  for (size_t i(0);i<m_obss.size();++i) delete m_obss[i];
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

void NLO_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
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
  for (size_t i=0; i<m_obss.size(); i++) {
    const double value = m_obss[i]->Value(mom, fl);
    msg_Debugging()<<"value["<<i<<"] = "<<value
                   <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
    FillHisto(i,value,weight,ncount,mode);
    }
  return;
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
