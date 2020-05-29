#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PDF/Main/ISR_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "Main/Resum.H"
#include "Tools/Key_Base.H"

// bad hack :(
#define private public
#include "AddOns/Python/MEProcess.H"

int main(int argc,char* argv[])
{
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  SHERPA::Sherpa *Generator(new SHERPA::Sherpa());
  // initialize the framework
  try {
    Generator->InitializeTheRun(argc,argv);

    // create a MEProcess instance
    MEProcess Process(Generator);
    Process.Initialize();

    auto showers = Generator->GetInitHandler()->GetShowerHandlers().at(PDF::isr::id::hard_process);
    RESUM::Resum* shower = static_cast<RESUM::Resum*>(showers->GetShower());

    // shower->AddObservable({"JetAngularities",{"tag:JA05","alpha:0.5","R:0.8"}},{0.08});
    // shower->AddObservable({"JetAngularities",{"tag:JA1","alpha:1","R:0.8"}},{0.08});
    // shower->AddObservable({"JetAngularities",{"tag:JA2","alpha:2","R:0.8"}},{0.08});

    shower->AddObservable({"SD_JetAngularities",{"tag:JA05","alpha:0.5","R:0.8","beta:0","zcut:0.1","ep:1/3"}},{0.08});
    shower->AddObservable({"SD_JetAngularities",{"tag:JA1","alpha:1","R:0.8","beta:0","zcut:0.1","ep:1/3"}},{0.08});
    shower->AddObservable({"SD_JetAngularities",{"tag:JA2","alpha:2","R:0.8","beta:0","zcut:0.1","ep:1/3"}},{0.08});


    for (size_t n=1;n<=Process.NumberOfPoints();++n) {
      // set momenta from file
      Process.SetMomenta(n);

      auto ampl = Process.GetAmp();
      ampl->SetMuR2(pow(0.8*500,2)); 
      ampl->SetMuQ2(pow(0.8*500,2));
      ampl->SetMuF2(pow(0.8*500,2));
      ampl->SetProc(Process.p_proc);
      ampl->SetProcs(Process.p_proc->AllProcs());
      msg_Out()<<*ampl<<"\n";
      shower->PrepareShower(ampl);
      shower->PerformShowers();
      for(auto res: shower->m_resNLL) {
        for(double r: res) msg_Out()<<r<<"\n";
        msg_Out()<<"\n";
      } 
      msg_Out()<<"\n";
    }
  }
  catch(ATOOLS::Exception exception) {
    std::terminate();
  }
  delete Generator;
#ifdef USING__MPI
  MPI_Finalize();
#endif
  return 0;
}

