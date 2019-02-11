#-*- mode: python-*- 
import os, distutils.spawn, subprocess
from functools import partial
vars = Variables('.SConstruct')
vars.Add(PathVariable('sherpa','path to sherpa',
    os.popen('Sherpa-config --prefix').read().rstrip() if
    distutils.spawn.find_executable('Sherpa-config') else
    '/path/to/sherpa',PathVariable.PathIsDir))
vars.Add(PathVariable('include','directories that will be passed to the compiler with the -I option.','.'))
env = Environment(variables=vars,
                  CPPPATH=['${sherpa}/include/SHERPA-MC',os.getcwd(),'${include}'])
vars.Add('CXX','The C++ Compiler',
    os.popen(env['sherpa']+'/bin/Sherpa-config --cxx').read().rstrip())
vars.Add('CXXFLAGS','The C++ Flags',['-g','-O2','-std=c++11'])
vars.Update(env)
Help(vars.GenerateHelpText(env))
vars.Save('.SConstruct',env)
env['ENV']=os.environ
if env['PLATFORM']=='darwin':
   env.Append(LINKFLAGS=['-Wl,-undefined','-Wl,dynamic_lookup','-L/opt/local/lib', '-lgmp', '-lgmpxx',
                         '-lmpfr'])

resumlib = env.SharedLibrary('SherpaResum',
	['Math/r8lib.cpp',
	'Math/c8lib.cpp',
	'Math/matexp.cpp',
	'Math/asa007.cpp',
	'Tools/CBasis.C',
	'Tools/CMetric_Base.C',
	'Tools/Hard_Matrix.C',
        'Tools/StringTools.C',
        'Tools/Reader.C',
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
	'Observables/Thrust.C',
	'Observables/HeavyJetMass.C',
	'Observables/Thrust_F_table.C',
	'Observables/Durham_3Jet_res.C'],
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
        'Tools/StringTools.C',
        'Tools/Reader.C',
        'Bases/QCD_Generic.C',  
	'Main/Comix_Interface.C',
        'RRatios/RRatios.C',
	'Main/Cluster_Definitions.C'])

def replace(target, source, env, old, new):
    with open(str(source[0]), "rt") as fin:
         with open(str(target[0]), "wt") as fout:
              for line in fin:
                  fout.write(line.replace(old, new))

def copy(target, source, env):
   subprocess.check_output(['cp', str(source[0]), str(target[0])])
                  
def make_exe(target, source, env, cp=copy):
   cp(target, source, env)
   subprocess.check_output(['chmod', 'ug+x', str(target[0])])

env.Command(target="Tools/Files.H", source="Tools/Files.H.in",
            action=partial(replace, old='@share_dir',
	    			    new=os.path.join(env.subst('${sherpa}'),
							'share/RESUM')))
env.Command(target='${sherpa}/bin/dat2yoda', source="Scripts/dat2yoda",
	    action=partial(make_exe,
                           cp=partial(replace,
                                      old="-*- mode: python-*-",
		                      new="!"+subprocess.check_output(['which',
                                                                       'python']))))
env.Install('${sherpa}/lib/SHERPA-MC', [resumlib,analysislib,rratiolib])
env.Install('${sherpa}/share/RESUM',['share/pre_calc','share/FFunctions'])
env.Alias('install', ['Tools/Files.H',
		      '${sherpa}/bin/dat2yoda',
                      '${sherpa}/share/RESUM',
                      '${sherpa}/lib/SHERPA-MC'])

