#-*- mode: python-*- 
import os, distutils.spawn, subprocess
from functools import partial

AddOption('--enable-fastjetcontrib',
          dest='fjc',
          action='store_true',
          default=False,
          help=('Enable Observables that are implemented using fastjet-contrib.'
                'Needs to be installed as libfastjetcontribfragile.'
                'Use fjcpath to set path if needed.'))

AddOption('--enable-yoda',
          dest='yoda',
          action='store_true',
          default=False,
          help=('Enable use of yoda.',
                'Use yodapath to set path if needed.'))


vars = Variables('.SConstruct')
vars.Add(PathVariable('sherpa','path to sherpa',
    os.popen('Sherpa-config --prefix').read().rstrip() if
    distutils.spawn.find_executable('Sherpa-config') else
    '/path/to/sherpa',PathVariable.PathIsDir))
vars.Add(PathVariable('include','directories that will be passed to the compiler with the -I option.','.'))
vars.Add(BoolVariable('fjc','Whether to enable fjcontrib.', GetOption('fjc')))
vars.Add(PathVariable('fjcpath','path to fastjetcontrib','${sherpa}'))
vars.Add(BoolVariable('yoda','Whether to enable yoda.', GetOption('yoda')))
vars.Add(PathVariable('yodapath','path to yoda','${sherpa}'))

env = Environment(variables=vars,
                  CPPPATH=['${sherpa}/include/SHERPA-MC',
                           os.getcwd(),
                           '${include}',
                           '${yodapath}/include',
                           '${fjcpath}/include'])

vars.Add('CXX','The C++ Compiler',
    os.popen(env['sherpa']+'/bin/Sherpa-config --cxx').read().rstrip())
vars.Add('CXXFLAGS','The C++ Flags',['-g','-O2','-std=c++11'])
vars.Update(env)
Help(vars.GenerateHelpText(env))
vars.Save('.SConstruct',env)
env['ENV']=os.environ
if env['PLATFORM']=='darwin':
   env.Append(LINKFLAGS=['-Wl,-undefined','-Wl,dynamic_lookup'])

fjc = env['fjc']
yoda = env['yoda']

env.Append(CCFLAGS=[] +
           (['-D USING_YODA'] if yoda else []) +
           (['-D USING_FJCONTRIB'] if fjc else []))

resumcommon = env.SharedLibrary('ResumCommon',
                                 ['Math/r8lib.cpp',
	                          'Math/c8lib.cpp',
	                          'Math/matexp.cpp',
                                  'Math/asa007.cpp',
                                  'Math/spline.cpp',
                                  'Math/linpack_d.cpp',
                                  'Math/blas0.cpp',
                                  'Math/blas1_d.cpp',
                                  'Tools/StringTools.C',
                                  'Tools/Key_Base.C'])
   
resumlib = env.SharedLibrary('SherpaResum',
	                     ['Math/InterpolationWrapper.C',
	                      'Tools/CBasis.C',
	                      'Tools/CMetric_Base.C',
	                      'Tools/Hard_Matrix.C',
                              'Tools/Reader.C',
	                      'Bases/QCD_Generic.C',
	                      'Main/Comix_Interface.C',
	                      'Main/Resum.C',
                              'Main/Cluster_Definitions.C'
                             ],
	                     LIBPATH=['${sherpa}/lib/SHERPA-MC'],
	                     RPATH=['${sherpa}/lib/SHERPA-MC'],
	                     LIBS=['ResumCommon'])

observables = ['Observables/Y2_IF.C',
	       'Observables/Y1_II.C',
	       'Observables/Thrust_FF.C',
               'Observables/ColorSinglet/Thrust.C',
	       'Observables/Thrust_IF.C',
	       'Observables/Thrust_II.C',
	       'Observables/TThrust_FF.C',
	       'Observables/Thrust.C',
	       'Observables/HeavyJetMass.C',
               'Observables/YN_CambridgeAachen.C',
               'Observables/YN_Durham.C',
               'Observables/CParameter.C',
               'Observables/EParameter.C',
               'Observables/ThrustMajor.C',
               'Observables/ThrustMinor.C',
               'Observables/WideBroad.C',
               'Observables/TotBroad.C',
               'Observables/HeavyMinusLightHemMass.C',
               'Observables/Oblateness.C']

obsFjcontrib = ['Observables/JetAngularities.C',
                'Observables/Algorithms/FastjetAlg.C']

analysislib = env.SharedLibrary('ResumAnalysis',
                                ['Analysis/Observable_Base.C',
	                         'Analysis/Observable_Selector.C',
                                 'Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.C',
                                 'Analysis/ChannelAlgorithms/KT2_ee.C',
                                 'Analysis/ChannelAlgorithms/KT2_pp.C',
                                 'Analysis/NLO_Analysis.C',
                                 # 'Analysis/Matching_Analysis.C',
                                 'Analysis/NLL_Analysis.C',
                                 'Analysis/Resum_Enhance_Observable.C',
                                 'Scales/Resum_Scale_Setter.C',
                                 'Scales/Resum_Scale_Setter_Durham.C',
                                 'FFunction/FFunctions.C'
                                ] + (observables +
                                     (obsFjcontrib if fjc else [])),
	                        LIBPATH=(['${sherpa}/lib/SHERPA-MC']
                                         + (['${fjcpath}/lib']
                                            if fjc else [])
                                         + (['${yodapath}/lib']
                                            if yoda else [])),
	                        RPATH=(['${sherpa}/lib/SHERPA-MC']
                                       + (['${fjcpath}/lib']
                                            if fjc else [])
                                       + (['${yodapath}/lib']
                                          if yoda else [])),
	                        LIBS=(['SherpaAnalyses',
                                      'SherpaAnalysis',
                                      'ResumCommon']
                                      + (['fastjetcontribfragile']
                                         if fjc else [])
                                      + (['YODA']
                                         if yoda else [])))


rratiolib = env.SharedLibrary('SherpaRRatios',
	                      ['Tools/CBasis.C',
	                       'Tools/CMetric_Base.C',
	                       'Tools/Hard_Matrix.C',
                               'Tools/Reader.C',
                               'Bases/QCD_Generic.C',  
	                       'Main/Comix_Interface.C',
                               'RRatios/RRatios.C',
	                       'Main/Cluster_Definitions.C'],
	                      LIBPATH=['${sherpa}/lib/SHERPA-MC'],
	                      RPATH=['${sherpa}/lib/SHERPA-MC'],
	                      LIBS=['ResumCommon'])

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

print env['sherpa']

env.Command(target="Tools/Files.H", source="Tools/Files.H.in",
            action=partial(replace, old='@share_dir',
	    			    new=os.path.abspath(os.path.join(env.subst('${sherpa}'),
							             'share/RESUM'))))
env.Command(target='${sherpa}/bin/dat2yoda', source="Scripts/dat2yoda",
	    action=partial(make_exe,
                           cp=partial(replace,
                                      old="-*- mode: python-*-",
		                      new="!"+subprocess.check_output(['which',
                                                                       'python']))))

env.Command(target='${sherpa}/bin/resum-combine', source="Scripts/resum-combine",
	    action=partial(make_exe,
                           cp=partial(replace,
                                      old="-*- mode: python-*-",
		                      new="!"+subprocess.check_output(['which',
                                                                       'python']))))

env.Command(target='${sherpa}/bin/resum-match', source="Scripts/resum-match",
	    action=partial(make_exe,
                           cp=partial(replace,
                                      old="-*- mode: python-*-",
		                      new="!"+subprocess.check_output(['which',
                                                                       'python']))))

env.Install('${sherpa}/lib/SHERPA-MC', [resumcommon,resumlib,analysislib] +(
   [rratiolib] if yoda else []))
env.Install('${sherpa}/share/RESUM',['share/pre_calc','share/FFunction'])
env.Alias('install', ['Tools/Files.H',
		      '${sherpa}/bin/dat2yoda',
                      '${sherpa}/bin/resum-combine',
                      '${sherpa}/bin/resum-match',
                      '${sherpa}/share/RESUM',
                      '${sherpa}/lib/SHERPA-MC'])
