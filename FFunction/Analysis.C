#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "Analysis/Observable_Base.H"
#include "FFunction/Observables.H"
#include "FFunction/FFunction.H"

#include "MultiplePrecision.H"
#include "Analysis/Observable_Base.H"




namespace RESUM {

  class FFunctionAnalysis: public ANALYSIS::Analysis_Base {
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base_MP::Ptr> m_obss;

    FFunction* p_ffunc;

  public:

    FFunctionAnalysis(const ANALYSIS::Argument_Matrix &params);

    void Evaluate(const ATOOLS::Blob_List & blobs,
		  double weight,double ncount);
    void Evaluate(double weight, double ncount,int mode) {}
    Primitive_Observable_Base * Copy() const;

  };// end of class FFunctionAnalysis

}// end of namespace RESUM

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"


#include <algorithm>

using namespace RESUM;
using namespace ANALYSIS;
using namespace ATOOLS;

FFunctionAnalysis::FFunctionAnalysis(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_FFunction";
  Data_Reader reader(",",";","!","=");
  Algebra_Interpreter *ip=reader.Interpreter();
  msg_Debugging()<<"Look for SHOWER_GENERATOR: "<<rpa->gen.Variable("SHOWER_GENERATOR")<<".\n";
  p_ffunc=(FFunction*)ToType<void*>(rpa->gen.Variable("SHOWER_GENERATOR"));
  if (dynamic_cast<FFunction*>(p_ffunc)==NULL)
    THROW(fatal_error,"FFunction plugin not loaded");
  for (size_t i(1);i<params.size();++i) {
    if (params[i].size()<5) continue;
    double xmin(ToType<double>(ip->Interprete(params[i][1])));
    double xmax(ToType<double>(ip->Interprete(params[i][2])));
    size_t nbin(ToType<size_t>(ip->Interprete(params[i][3])));
    int tp(HistogramType(params[i][4]));
    msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp
		   <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    std::string name = params[i][0];
    m_obss.push_back(std::shared_ptr<Observable_Base_MP>(Observable_Getter_MP::GetObject
                                (name,Observable_Key(name))));
    if (m_obss.back()==NULL)
      THROW(not_implemented,"No such observable: "+params[i][0]);
    m_histos.push_back(new Histogram(tp,xmin,xmax,nbin,params[i][0]));
    msg_Debugging()<<"Add "<<m_obss.back()->Name()<<".\n";
    p_ffunc->AddObservable(m_obss.back(),m_histos.back());
  }
}

void FFunctionAnalysis::Evaluate(const ATOOLS::Blob_List &blobs,
			double weight,double ncount)
{
  DEBUG_FUNC("");
  Particle_List *all(p_ana->GetParticleList(m_listname));
  if (all==NULL) AddZero(ncount,0);
  Blob *sb(blobs.FindFirst(btp::Shower));
  if (sb==NULL || sb->TypeSpec()!="RESUM") {
    static int mess(true);
    if (mess) msg_Error()<<METHOD<<"(): Invalid shower. Skip.\n";
    mess=false;
    AddZero(ncount,0);
  }
  for (size_t n(0);n<m_obss.size();++n) {
    const std::vector<double> &res(p_ffunc->Result(n));
    msg_Debugging()<<"Fill '"<<m_obss[n]<<"' -> "<<res<<"\n";
    FillHisto(n,int(res[0]),weight*res[1],ncount,0);
  }
}

Primitive_Observable_Base *FFunctionAnalysis::Copy() const 
{
  return new FFunctionAnalysis(m_params);
}

DECLARE_GETTER(FFunctionAnalysis,"FFunction",Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,FFunctionAnalysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new FFunctionAnalysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,FFunctionAnalysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
