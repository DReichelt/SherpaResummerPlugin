import math as m


value = 0
sumW_Sig = 1
sumW2_Sig = 2
err_Sig = 3
sumW_bSig = 4
sumW2_bSig = 5
err_bSig = 6
nEntries = 7


class Sigma(object):

    
    def __init__(self,values,AbsValues=False):
        self.AbsValues = AbsValues
        self.values = values


    @classmethod
    def from_file(cls, filename, toAbsValues=False):
        c = cls([])
        c.AbsValues = toAbsValues
        with open(filename,'r') as f:
            for line in f.readlines():
                if line.strip().startswith("#"): continue
                vals = line.split()
                if len(vals) != 8: continue
                c.values += [[float(val) for val in vals]]
        if toAbsValues:
            for i in range(len(c.values)):
                N = c.values[i][nEntries]
                c.values[i][sumW_Sig] *= N
                c.values[i][sumW_bSig] *= N
                c.values[i][sumW2_Sig] *= N
                c.values[i][sumW2_bSig] *= N
        return c

    @classmethod
    def from_yoda(cls, scatter):
        n = len(scatter.points)
        c = cls([[0]*8]*n)
        c.AbsValues = False
        sigma = scatter.points[-1].y
        sigma_err = (abs(scatter.points[-1].yErrs[0])
                     +abs(scatter.points[-1].yErrs[1]))/2.
        for i in range(n):
            p = scatter.points[i]
            c.values[i][value] = p.x
            c.values[i][sumW_Sig] = p.y
            c.values[i][err_Sig] = (abs(p.yErrs[0])+abs(p.yErrs[1]))/2.
            c.values[i][sumW_bSig] = sigma-p.y
            c.values[i][err_bSig] = m.sqrt(pow(sigma_err,2)+pow(c.values[i][err_Sig],2))
        return c
    
    def Sigma(self, v):
        for val in self.values:
            if abs(val[value]-v) < 1e-6:
                if self.AbsValues:
                    N = val[nEntries]
                    return (val[sumW_Sig]/N, val[err_Sig])
                else:
                    return (val[sumW_Sig], val[err_Sig])
        print("Value {} not found.".format(v))

    def barSigma(self, v):
        for val in self.values:
            if abs(val[value]-v) < 1e-6:
                if self.AbsValues:
                    N = val[nEntries]
                    return (val[sumW_bSig]/N, val[err_bSig])
                else:
                    return (val[sumW_bSig], val[err_bSig])
        print("Value {} not found.".format(v))

    def __add__(self,other):
        n = len(self.values)
        c = Sigma([[0]*8 for _ in range(n)])
        if self.AbsValues:
            assert other.AbsValues
            for i in range(n):
                vo = other.values[i]
                v = self.values[i]
                assert abs(v[value]-vo[value]) < 1e-6
                c.values[i][value] = v[value]
                c.values[i][sumW_Sig] = (v[sumW_Sig]+vo[sumW_Sig])
                c.values[i][sumW_bSig] = (v[sumW_bSig]+vo[sumW_bSig])
                c.values[i][sumW2_Sig] = (v[sumW2_Sig]+vo[sumW2_Sig])
                c.values[i][sumW2_bSig] = (v[sumW2_bSig]+vo[sumW2_bSig])
                c.values[i][nEntries] = (v[nEntries]+vo[nEntries])
                N = c.values[i][nEntries]
                c.values[i][err_Sig] = m.sqrt((c.values[i][sumW2_Sig]/N-pow(c.values[i][sumW_Sig]/N,2))/(N-1))
                c.values[i][err_bSig] = m.sqrt((c.values[i][sumW2_bSig]/N-pow(c.values[i][sumW_bSig]/N,2))/(N-1))
                c.AbsValues = True
        else:
            for i in range(n):
                vo = other.values[i]
                v = self.values[i]
                assert abs(v[value]-vo[value]) < 1e-6
                c.values[i][value] = v[value]
                c.values[i][sumW_Sig] = (v[sumW_Sig]+vo[sumW_Sig])
                c.values[i][err_Sig] = m.sqrt(pow(v[err_Sig],2)+pow(vo[err_Sig],2))
                c.values[i][sumW_bSig] = (v[sumW_bSig]+vo[sumW_bSig])
                c.values[i][err_bSig] = m.sqrt(pow(v[err_bSig],2)+pow(vo[err_bSig],2))
        return c

    def __getitem__(self,i):
        return self.values[i]
    
    def toYoda(self, path, title=None, mode="ln"):
        import yoda
        if title is None: title = path
        scatter1 = yoda.Scatter2D(path=path+"_Sigma",title=title+"_Sigma")
        scatter2 = yoda.Scatter2D(path=path+"_barSigma", title=title+"_Sigma")
        n = len(self.values)
        for i in range(n):
            v = self.values[i]
            if mode == "lin": x = v[value]
            elif mode == "log": x = m.log10(v[value])
            elif mode == "ln":  x = m.log(v[value])
            if self.AbsValues:
                N = v[nEntries]
                y = v[sumW_Sig]/N
            else:
                y = v[sumW_Sig]
            yerr = v[err_Sig]
            scatter1.addPoint(x,y,yerrs=(yerr,yerr))
            if self.AbsValues:
                y = v[sumW_bSig]/N
            else:
                y = v[sumW_bSig]
            yerr = v[err_bSig]
            scatter2.addPoint(x,y,yerrs=(yerr,yerr))
        return [scatter1,scatter2]


    def toFile(self,filename):
        with open(filename,'w') as f:
            f.write("# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n")
            for v in self.values:
                if self.AbsValues:
                    sig = v[sumW_Sig]/v[nEntries]
                    sig2 = v[sumW2_Sig]/v[nEntries]
                    sigb = v[sumW_bSig]/v[nEntries]
                    sigb2 = v[sumW2_bSig]/v[nEntries]
                else:
                    sig = v[sumW_Sig]
                    sig2 = v[sumW2_Sig]
                    sigb = v[sumW_bSig]
                    sigb2 = v[sumW2_bSig]
                f.write("{VAL} {SIG} {SIGW2} {SIGERR} {BSIG} {BSIGW2} {BSIGERR} {N}\n".format(VAL=v[value],
                                                                                              SIG=sig,
                                                                                              SIGW2=sig2,
                                                                                              SIGERR=v[err_Sig],
                                                                                              BSIG=sigb,
                                                                                              BSIGW2=sigb2,
                                                                                              BSIGERR=v[err_bSig],
                                                                                              N=v[nEntries]))

if __name__=='__main__':
    import argparse
    import sys
    import os
    import subprocess
    import yoda
    parser = argparse.ArgumentParser()
    parser.add_argument('--nll-file', '-nll',default='NLL.dat')
    parser.add_argument('--lo-file', '-lo',default='LO.dat')
    parser.add_argument('--nlo-rs-file', '-nloRS',default='NLO_RS.dat')
    parser.add_argument('--nlo-vi-file', '-nloVI',default='NLO_VI.dat')
    parser.add_argument('--nlo-file', '-nlo', default='NLO_VI.dat')
    parser.add_argument('--observable', '-obs', required=True)
    parser.add_argument('--rebin', default=1, type=int)
    parser.add_argument('--output-nll', '-NLL', default="NLL.yoda")
    parser.add_argument('--output-lo', '-LO', default="LO.yoda")
    parser.add_argument('--output-nlo', '-NLO', default="NLO.yoda")
    parser.add_argument('--channels', '-ch', action='append')
    parser.add_argument('--channel-fo-add', default=None)
    parser.add_argument('--match', action='store_true', default=False)
    
    args = parser.parse_args(sys.argv[1:])

    obs = args.observable

    # NLL files
    if os.path.exists(args.nll_file):
        inpath = "{NLL}/{ACC}/{OBS}_Sigma.dat"
        outpath = "/Resum/{OBS}_{ACC}"
        yodas = Sigma.from_file(
            inpath.format(NLL=args.nll_file,
                          ACC="NLL",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="NLL"))
        yodas += Sigma.from_file(
            inpath.format(NLL=args.nll_file,
                          ACC="LO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="LO"))
        yodas += Sigma.from_file(
            inpath.format(NLL=args.nll_file,
                          ACC="NLO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="NLO"))
        for ch in args.channels:
            inpath = "{NLL}/{ACC}/Channel_{CH}_{OBS}_Sigma.dat"
            outpath = "/Resum/{OBS}_{ACC}_Channel_{CH}"
        
            yodas += Sigma.from_file(
                inpath.format(NLL=args.nll_file,
                              ACC="NLL",
                              CH=ch,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              ACC="NLL",
                                                              CH=ch))
            yodas += Sigma.from_file(
                inpath.format(NLL=args.nll_file,
                            ACC="LO",
                              CH=ch,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              ACC="LO",
                                                              CH=ch))
            yodas += Sigma.from_file(
                inpath.format(NLL=args.nll_file,
                            ACC="NLO",
                              CH=ch,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              ACC="NLO",
                                                              CH=ch))
        yoda.writeYODA(yodas,args.output_nll)
    elif not os.path.exists(args.output_nll):
        args.output_nll = "NO"
        
    # LO files
    if os.path.exists(args.lo_file):
        yodas = []
        inpath = "{ANALYSIS}/{ACC}/{OBS}_Sigma.dat"
        outpath = "/Resum/{OBS}_{ACC}"

        yodas = Sigma.from_file(
            inpath.format(ANALYSIS=args.lo_file,
                          ACC="LO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="LO"))
        for ch in args.channels:
            if args.channel_fo_add is None:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{OBS}_Sigma.dat"
            else:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{ADD}_{OBS}_Sigma.dat"
            outpath = "/Resum/{OBS}_{ACC}_Channel_{CH}"
        
            yodas += Sigma.from_file(
                inpath.format(ANALYSIS=args.lo_file,
                              ACC="LO",
                              CH=ch,
                              ADD=args.channel_fo_add,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              ACC="LO",
                                                              CH=ch))
        yoda.writeYODA(yodas,args.output_lo)
    elif not os.path.exists(args.output_lo):
        args.output_lo = "NO"
        
    # NLO files
    if os.path.exists(args.nlo_rs_file) and os.path.exists(args.nlo_vi_file):
        yodas = []
        inpath = "{ANALYSIS}/{ACC}/{OBS}_Sigma.dat"
        outpath = "/Resum/{OBS}_{ACC}"

        rs = Sigma.from_file(
            inpath.format(ANALYSIS=args.nlo_rs_file,
                          ACC="NLO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="NLO"))
        vi = Sigma.from_file(
            inpath.format(ANALYSIS=args.nlo_rs_file,
                          ACC="NLO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="NLO"))
        yodas = (rs+vi).toYoda(outpath.format(OBS=obs,ACC="LO"))
        for ch in args.channels:
            if args.channel_fo_add is None:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{OBS}_Sigma.dat"
            else:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{ADD}_{OBS}_Sigma.dat"
            outpath = "/Resum/{OBS}_{ACC}_Channel_{CH}"

            rs = Sigma.from_file(
                inpath.format(ANALYSIS=args.nlo_rs_file,
                              ACC="NLO",
                              CH=ch,
                              ADD=args.channel_fo_add,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              CH=ch,
                                                              ACC="NLO"))
            vi = Sigma.from_file(
                inpath.format(ANALYSIS=args.nlo_rs_file,
                              ACC="NLO",
                              CH=ch,
                              ADD=args.channel_fo_add,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              CH=ch,
                                                              ACC="NLO"))
            yodas += (rs+vi).toYoda(outpath.format(OBS=obs,ACC="NLO"))
            
        yoda.writeYODA(yodas,args.output_nlo)
        
    elif os.path.exists(args.nlo_file):
        yodas = []
        inpath = "{ANALYSIS}/{ACC}/{OBS}_Sigma.dat"
        outpath = "/Resum/{OBS}_{ACC}"

        yodas = Sigma.from_file(
            inpath.format(ANALYSIS=args.nlo_file,
                          ACC="LO",
                          OBS=obs)).toYoda(outpath.format(OBS=obs,ACC="LO"))
        for ch in args.channels:
            if args.channel_fo_add is None:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{OBS}_Sigma.dat"
            else:
                inpath = "{ANALYSIS}/{ACC}/Channel_{CH}_{ADD}_{OBS}_Sigma.dat"
            outpath = "/Resum/{OBS}_{ACC}_Channel_{CH}"
        
            yodas += Sigma.from_file(
                inpath.format(ANALYSIS=args.lo_file,
                              ACC="NLO",
                              CH=ch,
                              ADD=args.channel_fo_add,
                              OBS=obs)).toYoda(outpath.format(OBS=obs,
                                                              ACC="NLO",
                                                              CH=ch))
        yoda.writeYODA(yodas,args.output_nlo)
    elif not os.path.exists(args.output_nlo):
        args.output_nlo = 'NO'
        
    if args.match:
        command = ['python','match.py',
                   '-lo',args.output_lo,
                   '-nlo',args.output_nlo,
                   '-nll',args.output_nll,
                   '-obs',obs,
                   '--rebin',str(args.rebin),
                   '-o', 'resum-plots/{OBS}'.format(OBS=obs),
                   '-ch']+args.channels
        print " ".join(command)
        subprocess.call(command)
