import yoda
import math as m
import glob
import argparse
import numpy as np
import sys


class Sigma(object):

    def __init__(self,args):
        self.args = args
        lo_file = ["LO.yoda"] if self.args.lo_file is None else args.lo_file 
        nlo_file = ["NLO.yoda"] if self.args.nlo_file is None else args.nlo_file 
        nll_file = ["NLL.yoda"] if self.args.nll_file is None else args.nll_file 

        obs = self.args.observable
        
        f_nll = yoda.readYODA(nll_file)
        Sigma_NLL = self.get_Sigma("/Resum/{OBS}_NLL_Sigma".format(OBS=obs),f_nll)
        self.xvals = np.array([p.x for p in Sigma_NLL.points])
        self.nll = np.array([p.y for p in Sigma_NLL.points])
        self.ch_nll = {}
        for ch in args.channels:
            try:
                Sigma_NLL_ch = self.get_Sigma("/Resum/{OBS}_NLL_Channel_{CH}_Sigma".format(OBS=obs,
                                                                                           CH=ch),f_nll)
            except KeyError:
                print "No channel {} in NLL file.".format(ch)
                continue
            self.match_bins(Sigma_NLL,Sigma_NLL_ch)
            self.ch_nll[ch] = np.array([p.y for p in Sigma_NLL_ch.points])

        Sigma_expLO = self.get_Sigma("/Resum/{OBS}_LO_Sigma".format(OBS=obs),f_nll)
        self.match_bins(Sigma_NLL,Sigma_expLO)
        self.expLO = np.array([p.y for p in Sigma_expLO.points])
        self.ch_expLO = {}
        for ch in self.ch_nll:
            Sigma_expLO = self.get_Sigma("/Resum/{OBS}_LO_Channel_{CH}_Sigma".format(OBS=obs,
                                                                                  CH=ch),f_nll)
            self.match_bins(Sigma_NLL,Sigma_expLO)
            self.ch_expLO[ch] = np.array([p.y for p in Sigma_expLO.points])
            
        Sigma_expNLO = self.get_Sigma("/Resum/{OBS}_NLO_Sigma".format(OBS=obs),f_nll)
        self.match_bins(Sigma_NLL,Sigma_expNLO)
        self.expNLO = np.array([p.y for p in Sigma_expNLO.points])
        self.ch_expNLO = {}
        for ch in self.ch_nll:
            Sigma_expNLO = self.get_Sigma("/Resum/{OBS}_NLO_Channel_{CH}_Sigma".format(OBS=obs,
                                                                                       CH=ch),f_nll)
            self.match_bins(Sigma_NLL,Sigma_expNLO)
            self.ch_expNLO[ch] = np.array([p.y for p in Sigma_expNLO.points])

        self.sigma_0 = self.nll[-1]
        self.ch_sigma_0 = {}
        for ch in self.args.channels:
            self.ch_sigma_0[ch] = self.ch_nll[ch][-1] if ch in self.ch_nll else 0

        self.ch_lo = {}
        try:
            f_lo = yoda.readYODA(lo_file)
        except:
            print "No LO file ({FILE}).".format(FILE=lo_file)
            self.lo = self.expLO
            for ch in args.channels:
                self.ch_lo[ch] = self.ch_expLO[ch]
        else:
            Sigma_LO = self.get_Sigma("/Resum/{OBS}_LO_Sigma".format(OBS=obs),f_lo)
            self.match_bins(Sigma_NLL,Sigma_LO)
            self.lo = np.array([p.y-self.sigma_0 for p in Sigma_LO.points])
            for ch in args.channels:
                Sigma_LO = self.get_Sigma("/Resum/{OBS}_LO_Channel_{CH}_Sigma".format(OBS=obs,CH=ch),f_lo)
                self.match_bins(Sigma_NLL,Sigma_LO)
                self.ch_lo[ch] = np.array([p.y-self.ch_sigma_0[ch] for p in Sigma_LO.points])

        self.ch_nlo = {}
        try:
            f_nlo = yoda.readYODA(nlo_file)
        except:
            print "No NLO file ({FILE}).".format(FILE=nlo_file)
            self.nlo = self.expNLO
            for ch in args.channels:
                self.ch_nlo[ch] = self.ch_expNLO[ch]
        else:
            Sigma_NLO = self.get_Sigma("/Resum/{OBS}_NLO_barSigma".format(OBS=obs),f_nlo)
            self.match_bins(Sigma_NLL,Sigma_NLO)
            self.nlo = np.array([p.y for p in Sigma_NLO.points])

            for ch in args.channels:
                Sigma_NLO = self.get_Sigma("/Resum/{OBS}_NLO_Channel_{CH}_barSigma".format(OBS=obs,
                                                                                           CH=ch),f_nlo)
                self.match_bins(Sigma_NLL,Sigma_NLO)
                self.ch_nlo[ch] = np.array([p.y for p in Sigma_NLO.points])

        self.aSC = (self.lo-self.expLO)/self.sigma_0
        self.ch_aSC = {ch : ((self.ch_lo[ch]
                              -self.ch_expLO[ch])/self.ch_sigma_0[ch]
                             if ch in self.ch_nll else np.zeros_like(self.aSC))
                       for ch in self.args.channels}


        # Additive
        self._nll_add = self.match_add(self.nll,
                                       self.lo,
                                       self.expLO,
                                       self.sigma_0,
                                       aSC=self.aSC[0],
                                       nlo=self.nlo,
                                       expnlo=self.expNLO)

        self._nll_add_noC = self.match_add(self.nll,
                                           self.lo,
                                           self.expLO,
                                           self.sigma_0,
                                           aSC=0.,
                                           nlo=self.nlo,
                                           expnlo=self.expNLO)

        
        self.ch_nll_add = {ch : (self.match_add_ch(ch) if ch in self.ch_nll
                                 else np.zeros_like(self._nll_add))
                           for ch in self.args.channels}

        self.ch_nll_add_noC = {ch : (self.match_add(self.ch_nll[ch],
                                                    self.ch_lo[ch],
                                                    self.ch_expLO[ch],
                                                    self.ch_sigma_0[ch],
                                                    aSC=0.,
                                                    nlo=self.ch_nlo[ch],
                                                    expnlo=self.ch_expNLO[ch])
                                     if ch in self.ch_nll
                                     else np.zeros_like(self._nll_add))
                               for ch in self.args.channels}
        
        # Multiplicative
        self._nll_mult = self.match_mult(self.nll,
                                         self.lo,
                                         self.expLO,
                                         self.sigma_0,
                                         nlo=self.nlo,
                                         expnlo=self.expNLO)
        self.ch_nll_mult = {ch : (self.match_mult_ch(ch) if ch in self.ch_nll
                                  else np.zeros_like(self._nll_mult))
                            for ch in self.args.channels}

        # LogR
        self._nll_logr = self.match_logr(self.nll,
                                         self.lo,
                                         self.expLO,
                                         self.sigma_0,
                                         nlo=self.nlo,
                                         expnlo=self.expNLO)
        self.ch_nll_logr = {ch : (self.match_logr_ch(ch) if ch in self.ch_nll
                                  else np.zeros_like(self._nll_logr))
                            for ch in self.args.channels}

        self.ch_other = {ch: np.array([self.ch_lo[ch], self.ch_lo[ch]+self.ch_nlo[ch]])
                         if (not ch in self.ch_nll) else 0. for ch in self.args.channels}

        if len(self.args.channels) > 0:
            self.nll_add_noC = sum(self.ch_nll_add_noC.values()+self.ch_other.values())
            self.nll_add = sum(self.ch_nll_add.values()+self.ch_other.values())
            self.nll_mult = sum(self.ch_nll_mult.values()+self.ch_other.values())
            self.nll_logr = sum(self.ch_nll_logr.values()+self.ch_other.values())
        
            # Add missing channels to nll, with 0 contribution
            for ch in self.args.channels:
                if not ch in self.ch_nll:
                    self.ch_nll[ch] = np.zeros_like(self.nll)
                    self.ch_expLO[ch] = np.zeros_like(self.nll)
                    self.ch_expNLO[ch] = np.zeros_like(self.nll)
        else:
            self.nll_add = self._nll_add
            self.nll_add_noC = self._nll_add_noC
            self.nll_mult = self._nll_mult
            self.nll_logr = self._nll_logr
            
    def toYODA(self,histo,normalize=True, scale=1.,variations=[],rebin=None):
        path = "/ResumAnalysis/{OBS}".format(OBS=self.args.observable)
        Sigma = yoda.Scatter2D(path+"_Sigma",path+"_Sigma")
        dSigma = yoda.Scatter2D(path,path)
        if rebin is None:
            rebin = self.args.rebin
        xerr = self.xvals[::rebin][1:]-self.xvals[::rebin][:-1]
        deriv = (histo[1:]-histo[:-1])/xerr
        varder = np.empty((len(variations)*len(variations),len(deriv)))

        for i in range(len(variations)):
            for j in range(len(variations)):
                varder[i*len(variations)+j] =  ([variations[i][1:]-variations[j][:-1]])/xerr
        if len(variations) > 0:
            vardown = np.array([min(varder[:,i]) for i in range(len(deriv))])
            varup = np.array([max(varder[:,i]) for i in range(len(deriv))])
        else:
            vardown = np.array([d for d in deriv])
            varup = np.array([d for d in deriv])
        if normalize:
            if abs(histo[-1]) < 1e-10:
                print "Can not normalize."
                hist = np.array([h for h in histo])
            else:
                deriv /= histo[-1]
                vardown /= histo[-1]
                varup /= histo[-1]
                hist = histo/histo[-1]
            deriv *= scale
            vardown *= scale
            varup *= scale
            hist *= scale
        else:
            hist = scale*histo
            deriv *= scale
            vardown *= scale
            varup *= scale
        yup = np.abs(varup-deriv)
        ydown = np.abs(vardown-deriv)
        for i in range(len(self.xvals[::rebin])):
            x = self.xvals[::rebin][i]
            Sigma.addPoint(x,hist[i])
            if i < len(deriv):
                xe = xerr[i]/2.
                dSigma.addPoint(x+xe,deriv[i],xerrs=(xe,xe),
                                yerrs=(ydown[i],yup[i]))
        return [dSigma,Sigma]

    def match_add_ch(self,ch):
        return self.match_add(self.ch_nll[ch],
                              self.ch_lo[ch],
                              self.ch_expLO[ch],
                              self.ch_sigma_0[ch],
                              aSC=self.ch_aSC[ch][0],
                              nlo=self.ch_nlo[ch],
                              expnlo=self.ch_expNLO[ch])
    
    def match_add(self,nll,lo,explo,s0,aSC=0.,nlo=None,expnlo=None):
        LO = (1.+aSC)*nll+lo-explo-aSC*s0
        if nlo is None: return LO
        NLO = -nlo - expnlo - aSC*explo
        return np.array([LO, LO+NLO])

    def match_mult_ch(self,ch):
        return self.match_mult(self.ch_nll[ch],
                               self.ch_lo[ch],
                               self.ch_expLO[ch],
                               self.ch_sigma_0[ch],
                               nlo=self.ch_nlo[ch],
                               expnlo=self.ch_expNLO[ch])

    
    def match_mult(self,nll,lo,explo,s0,aSC=0.,nlo=None,expnlo=None):
        LO = nll*(1.+(lo-explo)/s0)
        if nlo is None: return LO
        NLO = nll*(-nlo-expnlo-explo*(lo-explo)/s0)/s0
        return np.array([LO, LO+NLO])

    def match_logr_ch(self,ch):
        return self.match_logr(self.ch_nll[ch],
                               self.ch_lo[ch],
                               self.ch_expLO[ch],
                               self.ch_sigma_0[ch],
                               nlo=self.ch_nlo[ch],
                               expnlo=self.ch_expNLO[ch])

    
    def match_logr(self,nll,lo,explo,s0,asC=0.,nlo=None,expnlo=None):
        LO = nll*np.exp((lo-explo)/s0)
        if nlo is None: return LO
        NLO = np.exp((-nlo-expnlo-
                     1./2.*(np.power(lo,2)-np.power(explo,2))/s0)/s0)
        return np.array([LO, LO*NLO])
        
            
        
    def get_Sigma(self,name, yodafile):
        return yodafile[name]

    def match_bins(self,nom, other):
        assert len(nom.points) == len(other.points)
        for i in range(len(nom.points)):
            if(abs(nom.points[i].x-other.points[i].x)<1e-4):
                "Warning: matched unaligned bins {} vs. {}".format(nom.points[i].x,
                                                                   other.points[i].x)
            other.points[i].x = nom.points[i].x



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lo-file', '-lo')
    parser.add_argument('--nlo-file', '-nlo')
    parser.add_argument('--nll-file', '-nll')
    parser.add_argument('--observable', '-obs', required=True)
    parser.add_argument('--rebin', default=1, type=int)
    parser.add_argument('--output-dir', '-o', default="resum-plots")
    parser.add_argument('--file-prefix', '-O', default="Resum")
    parser.add_argument('--channels', '-ch', default=None, nargs='+')
    parser.add_argument('--scale', '-s', default=None, type=float)
    
    args = parser.parse_args(sys.argv[1:])
    if args.channels is None: args.channels = []
    nll_file = np.array([args.nll_file.split(":")])
    print nll_file
    args.nll_file = nll_file[:,0]
    sigma = Sigma(args)
    if args.scale is None: args.scale = 1./sigma.sigma_0
    variations = []
    for i in range(len(nll_file[0])):
        args.nll_file = nll_file[:,i]
        print args.nll_file
        variations += [Sigma(args)]
    print variations

    
    import os

    file_prefix = args.file_prefix
    dirName = args.output_dir
    if not os.path.exists(dirName): os.makedirs(dirName)
    yoda.writeYODA(sigma.toYODA(sigma.nll_mult[0][::args.rebin],variations=[s.nll_mult[0][::args.rebin] for
                                                                            s in variations]),
                   "{DIR}/{FILE_PREF}_mult_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll_logr[0][::args.rebin], variations=[s.nll_logr[0][::args.rebin] for
                                                               s in variations]),
                   "{DIR}/{FILE_PREF}_LogR_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll_add[0][::args.rebin], variations=[s.nll_add[0][::args.rebin] for
                                                              s in variations]),
                   "{DIR}/{FILE_PREF}_add_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll_add_noC[0][::args.rebin],
                                variations=[s.nll_add_noC[0][::args.rebin]
                                            for s in variations]), 
                   "{DIR}/{FILE_PREF}_add_LO_noC.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    yoda.writeYODA(sigma.toYODA(sigma.nll_mult[1][::args.rebin], variations=[s.nll_mult[1][::args.rebin] for
                                                               s in variations]), 
                   "{DIR}/{FILE_PREF}_mult_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll_logr[1][::args.rebin], variations=[s.nll_logr[1][::args.rebin] for
                                                               s in
                                                               variations]),
                   "{DIR}/{FILE_PREF}_LogR_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll_add[1][::args.rebin], variations=[s.nll_mult[1][::args.rebin] for
                                                              s in variations]),
                   "{DIR}/{FILE_PREF}_add_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA(sigma.nll[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_NLL.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+sigma.lo)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+sigma.lo-sigma.nlo)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+sigma.expLO)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_expLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+sigma.expLO+sigma.expNLO)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_expNLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+(1.+sigma.aSC)*sigma.expLO+sigma.expNLO)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_expNLO_Capprox.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+(1.+sigma.aSC[0])*sigma.expLO+sigma.expNLO)[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_expNLO_C.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    expNLO = np.zeros_like(sigma.expLO)
    for ch in args.channels:
        expNLO += sigma.ch_sigma_0[ch]+(1.+sigma.ch_aSC[ch])*sigma.ch_expLO[ch]+sigma.ch_expNLO[ch]
    yoda.writeYODA(sigma.toYODA(expNLO[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_expNLO_Channels.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
    
    yoda.writeYODA(sigma.toYODA((sigma.sigma_0+(1.+sigma.aSC[0])*sigma.expLO+sigma.expNLO-(sigma.sigma_0+sigma.lo-sigma.nlo))[::args.rebin],normalize=False,scale=args.scale),
                   "{DIR}/{FILE_PREF}_NLOdiff.yoda".format(DIR=dirName,FILE_PREF=file_prefix))

    # print sigma.nll_logr[0][-1]
    # for p in sigma.toYODA(sigma.nll_logr[0], variations=[s.nll_logr[0] for
    #                                                      s in variations],normalize=True)[0].points:
    #     print p.x, p.y
    for ch in args.channels:
        # print sigma.ch_nll_logr[ch][0][-1]
        # for p in sigma.toYODA(sigma.ch_nll_logr[ch][0],normalize=True)[0].points:
        #     print p.x, p.y
        dirName = args.output_dir+"-{CH}".format(CH=ch)
        if not os.path.exists(dirName): os.mkdir(dirName)
        print ch, sigma.ch_sigma_0[ch]
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_mult[ch][0][::args.rebin],variations=[s.ch_nll_mult[ch][0][::args.rebin]
                                                                         for s in
                                                                         variations]),
                       "{DIR}/{FILE_PREF}_mult_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_logr[ch][0][::args.rebin],variations=[s.ch_nll_logr[ch][0][::args.rebin]
                                                                         for s in
                                                                         variations]),
                       "{DIR}/{FILE_PREF}_LogR_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_add[ch][0][::args.rebin],variations=[s.ch_nll_add[ch][0][::args.rebin]
                                                                        for s in
                                                                        variations]),
                       "{DIR}/{FILE_PREF}_add_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_mult[ch][1][::args.rebin],variations=[s.ch_nll_mult[ch][1][::args.rebin]
                                                                         for s in
                                                                         variations]),
                       "{DIR}/{FILE_PREF}_mult_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_logr[ch][1][::args.rebin],variations=[s.ch_nll_logr[ch][1][::args.rebin]
                                                                         for s in
                                                                         variations]),
                       "{DIR}/{FILE_PREF}_LogR_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll_add[ch][1][::args.rebin],variations=[s.ch_nll_add[ch][1][::args.rebin]
                                                                        for s in
                                                                        variations]),
                       "{DIR}/{FILE_PREF}_add_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        scale = 1./sigma.ch_sigma_0[ch] if sigma.ch_sigma_0[ch] > 0 else 1.
        yoda.writeYODA(sigma.toYODA(sigma.ch_nll[ch][::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_NLL.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+sigma.ch_lo[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_LO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+sigma.ch_lo[ch]-sigma.ch_nlo[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_NLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+sigma.ch_expLO[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_expLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+sigma.ch_expLO[ch]+sigma.ch_expNLO[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_expNLO.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+(1.+sigma.ch_aSC[ch])*sigma.ch_expLO[ch]+sigma.ch_expNLO[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_expNLO_Capprox.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
        yoda.writeYODA(sigma.toYODA((sigma.ch_sigma_0[ch]+(1.+sigma.ch_aSC[ch][0])*sigma.ch_expLO[ch]+sigma.ch_expNLO[ch])[::args.rebin],normalize=False,scale=scale),
                       "{DIR}/{FILE_PREF}_expNLO_C.yoda".format(DIR=dirName,FILE_PREF=file_prefix))
                    
