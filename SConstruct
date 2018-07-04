import os, distutils.spawn
vars = Variables('.SConstruct')
vars.Add(PathVariable('sherpa','path to sherpa',
    os.popen('Sherpa-config --prefix').read().rstrip() if
    distutils.spawn.find_executable('Sherpa-config') else
    '/path/to/sherpa',PathVariable.PathIsDir))
env = Environment(variables=vars,
    CPPPATH=['${sherpa}/include/SHERPA-MC',os.getcwd()])
vars.Add('CXX','The C++ Compiler',
    os.popen(env['sherpa']+'/bin/Sherpa-config --cxx').read().rstrip())
vars.Add('CXXFLAGS','The C++ Flags',['-g','-O2','-std=c++11'])
vars.Update(env)
Help(vars.GenerateHelpText(env))
vars.Save('.SConstruct',env)
env['ENV']=os.environ
if env['PLATFORM']=='darwin':
   env.Append(LINKFLAGS=['-Wl,-undefined','-Wl,dynamic_lookup'])

resumlib = env.SharedLibrary('SherpaResum',
	['Math/r8lib.cpp',
	'Math/c8lib.cpp',
	'Math/matexp.cpp',
	'Math/asa007.cpp',
	'Tools/CBasis.C',
	'Tools/CMetric_Base.C',
	'Tools/Hard_Matrix.C',
	'Bases/QCD_Generic.C',
	'Main/Comix_Interface.C',
	'Main/Resum.C',
	'Main/Cluster_Definitions.C'])

analysislib = env.SharedLibrary('ResumAnalysis',
	['Analysis/Analysis.C',
	'Analysis/Observable_Base.C',
	'Analysis/Matching_Analysis.C',
	'Analysis/Observable_Selector.C',
	'Observables/Y3_FF.C',
	'Observables/Y2_IF.C',
	'Observables/Y1_II.C',
	'Observables/Thrust_FF.C',
	'Observables/Thrust_IF.C',
	'Observables/Thrust_II.C',
	'Observables/TThrust_FF.C',
	'Observables/Thrust.C'],
	LIBPATH=['${sherpa}/lib/SHERPA-MC'],
	RPATH=['${sherpa}/lib/SHERPA-MC'],
	LIBS=['SherpaAnalyses','SherpaAnalysis'])

rratiolib = env.SharedLibrary('SherpaRRatios',
	['Math/r8lib.cpp',
	'Math/c8lib.cpp',
	'Math/matexp.cpp',
	'Math/asa007.cpp',
	'Tools/CBasis.C',
	'Tools/CMetric_Base.C',
	'Tools/Hard_Matrix.C',
	'Bases/QCD_Generic.C',
	'Main/Comix_Interface.C',
        'RRatios/RRatios.C',
	'Main/Cluster_Definitions.C'])

env.Install('${sherpa}/lib/SHERPA-MC', [resumlib,analysislib,rratiolib])
env.Alias('install', '${sherpa}/lib/SHERPA-MC')
