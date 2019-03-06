#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "Analysis/Observable_Base.H"
#include "Main/Resum.H"

namespace RESUM {

  class Analysis: public ANALYSIS::Analysis_Base {
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    Resum *p_resum;

  public:

    Analysis(const ANALYSIS::Argument_Matrix &params);

    void Evaluate(const ATOOLS::Blob_List & blobs,
		  double weight,double ncount);
    void Evaluate(double weight, double ncount,int mode) {}
    Primitive_Observable_Base * Copy() const;

  };// end of class Analysis

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

Analysis::Analysis(const Argument_Matrix& params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_Resum";
  Data_Reader reader(",",";","!","=");
  Algebra_Interpreter* ip = reader.Interpreter();
  p_resum=(Resum*)ToType<void*>(rpa->gen.Variable("SHOWER_GENERATOR"));

  if (dynamic_cast<Resum*>(p_resum) == nullptr) {
    THROW(fatal_error,"Resummer plugin not loaded");
  }
  for (size_t i=1; i<params.size(); i++) {
    if (params[i].size()<5) continue;
    double xmin = ToType<double>(ip->Interprete(params[i][1]));
    double xmax = ToType<double>(ip->Interprete(params[i][2]));
    size_t nbin = ToType<size_t>(ip->Interprete(params[i][3]));
    int tp = HistogramType(params[i][4]);
    msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp
		   <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    m_obss.push_back(Observable_Getter::GetObject
		     (params[i][0],Observable_Key(params[i][0])));
    if (m_obss.back() == nullptr)
      THROW(not_implemented,"No such observable: "+params[i][0]);
    m_histos.push_back(new Histogram(tp,xmin,xmax,nbin,params[i][0]));
    p_resum->AddObservable(m_obss.back(),m_histos.back());
  }
}

void Analysis::Evaluate(const ATOOLS::Blob_List& blobs,
			double weight, double ncount)
{
  DEBUG_FUNC("");
  Particle_List* all = p_ana->GetParticleList(m_listname);
  if (all == nullptr) AddZero(ncount,0);
  Blob* sb = blobs.FindFirst(btp::Shower);
  if (sb == nullptr || sb->TypeSpec()!="RESUM") {
    /// @TODO: what is that good for?
    static int mess(true);
    if (mess) msg_Error()<<METHOD<<"(): Invalid shower. Skip.\n";
    mess=false;
    AddZero(ncount,0);
  }
  for (size_t n=0; n<m_obss.size(); n++) {
    const std::pair<int, double>& res = p_resum->Result(n);
    msg_Debugging()<<"Fill '"<<m_obss[n]<<"' -> "<<res.first<<" "<<res.second<<"\n";
    FillHisto(n, res.first, weight*res.second, ncount, 0);
  }
}

Primitive_Observable_Base *Analysis::Copy() const 
{
  return new Analysis(m_params);
}

DECLARE_GETTER(Analysis,"Resum",Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,Analysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
