import optparse
import numpy
import ROOT
import os.path
import math

from LatinoAnalysis.Gardener.gardening import TreeCloner

#
#  _ _|      |         /     _ _|                    ___|                |            ____|             |                           
#    |    _` |        /        |    __|   _ \      \___ \    __|   _` |  |   _ \      |     _` |   __|  __|   _ \    __|  __| 
#    |   (   |       /         |  \__ \  (   |           |  (     (   |  |   __/      __|  (   |  (     |    (   |  |   \__ \ 
#  ___| \__,_|     _/        ___| ____/ \___/      _____/  \___| \__,_| _| \___|     _|   \__,_| \___| \__| \___/  _|   ____/ 
#                                                                                                                             
#

class IdIsoSFFiller(TreeCloner):

    def __init__(self):
        pass

    def __del__(self):
        pass

    def help(self):
        return '''Add a new lepton scale factor weight based on id/isolation scale factors data/MC.'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)

        group.add_option('-c', '--cmssw', dest='cmssw', help='cmssw version (naming convention may change)', default='763', type='string')

        group.add_option( '--idmu',       dest='idScaleFactorsFileMu' ,       help='file with scale factors for id for muons', default=None)
        group.add_option( '--isoTightmu', dest='isoTightScaleFactorsFileMu' , help='file with scale factors for isolation ,tight definition, for muons', default=None)
        group.add_option( '--isoLoosemu', dest='isoLooseScaleFactorsFileMu' , help='file with scale factors for isolation ,loose definition, for muons', default=None)
        
        group.add_option( '--isoidele'   , dest='idIsoScaleFactorsFileElectron',            help='file with scale factors for isolation and id for electrons',                  default=None)
        group.add_option( '--tkSCele'    , dest='tkSCFileElectron',                         help='file with scale factors for track-SC efficiency for electrons',               default=None)
        group.add_option( '--isoideleAlt', dest='idIsoScaleFactorsFileElectronAlternative', help='file with scale factors for isolation and id for electrons, alternative',     default=None)
        group.add_option( '--isoideleAltLumiRatio', dest='idIsoScaleFactorsFileElectronAlternativeLumiRatio', help='Luminosity ratio between first period and the whole', type='float'  ,    default=-1.0)

        parser.add_option_group(group)
        return group



    def checkOptions(self,opts):
       
        # ~~~~
        idIsoScaleFactors = {}

        cmssw_base = os.getenv('CMSSW_BASE')
        if opts.idScaleFactorsFileMu == None :
          if opts.cmssw == "ICHEP2016" :  opts.idScaleFactorsFileMu =        cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/muons.txt'  
          else :                          opts.idScaleFactorsFileMu =        cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/muons_Moriond76x.txt'
        if opts.isoTightScaleFactorsFileMu == None :
          if opts.cmssw == "ICHEP2016" :  opts.isoTightScaleFactorsFileMu = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/muons_iso_tight.txt'  
          else :                          opts.isoTightScaleFactorsFileMu = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/muons_iso_tight_Moriond76x.txt'
        if opts.isoLooseScaleFactorsFileMu == None :
          if opts.cmssw == "ICHEP2016" :  opts.isoLooseScaleFactorsFileMu = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/muons_iso_loose.txt'  
          else :                          opts.isoLooseScaleFactorsFileMu = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/muons_iso_loose_Moriond76x.txt'

        if opts.idIsoScaleFactorsFileElectron == None :
          if opts.cmssw == "ICHEP2016" :  opts.idIsoScaleFactorsFileElectron = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/electrons.txt'   
          else :                          opts.idIsoScaleFactorsFileElectron = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/electrons_Moriond76x.txt' 
        
        if opts.tkSCFileElectron == None :
          opts.tkSCFileElectron = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/eleRECO.txt.egamma_SF2D.root'
 
        if opts.cmssw == "ICHEP2016" :  opts.idIsoScaleFactorsFileElectronAlternative = cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/electrons_firstPart.txt'   
 
           
        file_idScaleFactorsFileMu  = open (opts.idScaleFactorsFileMu)
        file_isoTightScaleFactorsFileMu  = open (opts.isoTightScaleFactorsFileMu)
        file_isoLooseScaleFactorsFileMu  = open (opts.isoLooseScaleFactorsFileMu)

        file_idIsoScaleFactorsFileElectron = open (opts.idIsoScaleFactorsFileElectron)
        if opts.cmssw == "ICHEP2016" :  file_idIsoScaleFactorsFileElectronAlternative = open (opts.idIsoScaleFactorsFileElectronAlternative)


        self.idIsoScaleFactors = {}
        #                                       create the list               from the line                                if there is no "#"
        self.idIsoScaleFactors['ele']    =    [line.rstrip().split()    for line in file_idIsoScaleFactorsFileElectron     if '#' not in line]
        if opts.cmssw == "ICHEP2016" : 
          self.idIsoScaleFactors['eleAlt'] =    [line.rstrip().split()    for line in file_idIsoScaleFactorsFileElectronAlternative     if '#' not in line]
        
        self.idIsoScaleFactors['mu']    =    [line.rstrip().split()    for line in file_idScaleFactorsFileMu              if '#' not in line]

        self.isoScaleFactors = {}
        self.isoScaleFactors['muTight']   =    [line.rstrip().split()    for line in file_isoTightScaleFactorsFileMu        if '#' not in line]
        self.isoScaleFactors['muLoose']   =    [line.rstrip().split()    for line in file_isoLooseScaleFactorsFileMu        if '#' not in line]

        self.tkSCElectronRootFile = self._openRootFile(opts.tkSCFileElectron)
        self.tkSCElectronHisto = self._getRootObj(self.tkSCElectronRootFile, 'EGamma_SF2D')



        self.cmssw = opts.cmssw
        self.idIsoScaleFactorsFileElectronAlternativeLumiRatio = opts.idIsoScaleFactorsFileElectronAlternativeLumiRatio


        
        self.minpt_mu = 10
        self.maxpt_mu = 200
        self.mineta_mu = -2.4
        self.maxeta_mu = 2.4
        
        self.minpt_ele = 10
        self.maxpt_ele = 200
        self.mineta_ele = -2.5
        self.maxeta_ele = 2.5


    def _getHistoValue(self, h2, pt, eta):

        nbins = h2.GetNbinsY()
        ptmax = -1
        if (ptmax <= 0.) : 
          ptmax = h2.GetYaxis().GetBinCenter(nbins)
        
        # eta on x-axis, pt on y-axis
        value = h2.GetBinContent(h2.FindBin(eta, min(pt, ptmax)))
        error = h2.GetBinError  (h2.FindBin(eta, min(pt, ptmax)))
        
        #print ' x,y(max),z,err = ', eta, ' - ', min(pt, ptmax), '(', ptmax, ') - ', value, ' - ', error
        return value, error



    #                                              wantOnlyRecoEff = 0 ( idiso scale factors only ), 1 (reco scale factors only), 2 ( idiso * reco scale factors )
    def _getWeight (self, kindLep, pt, eta, tight, wantOnlyRecoEff):

        # fix underflow and overflow

        # print " kindLep = ", kindLep
        
        if kindLep == 'ele' :          
          if pt < self.minpt_ele:
            pt = self.minpt_ele
          if pt > self.maxpt_ele:
            pt = self.maxpt_ele
          
          if eta < self.mineta_ele:
            eta = self.mineta_ele
          if eta > self.maxeta_ele:
            eta = self.maxeta_ele

        if kindLep == 'mu' :          
          if pt < self.minpt_mu:
            pt = self.minpt_mu
          if pt > self.maxpt_mu:
            pt = self.maxpt_mu
          
          if eta < self.mineta_mu:
            eta = self.mineta_mu
          if eta > self.maxeta_mu:
            eta = self.maxeta_mu
 
 
        #print " self.idIsoScaleFactors = ", self.idIsoScaleFactors
        
        # decide if to use the first period of 2016 electron data
        # or the second period
        toss_a_coin = 1.
        if self.cmssw == "ICHEP2016" : 
          toss_a_coin = ROOT.gRandom.Rndm()
          if kindLep == 'ele' :
            if toss_a_coin < self.idIsoScaleFactorsFileElectronAlternativeLumiRatio: 
              kindLep == 'eleAlt'
        
        # idiso * reco scale factors
        if wantOnlyRecoEff == 2 or wantOnlyRecoEff == 0:
          if kindLep in self.idIsoScaleFactors.keys() : 
            #print " self.idIsoScaleFactors = ", self.idIsoScaleFactors
            #print " eta,pt = ",eta, ", ", pt
            # get the efficiency
            if kindLep == 'ele' :
              #print " self.idIsoScaleFactors[", kindLep, "] = ", self.idIsoScaleFactors[kindLep]
              
              tkSC, tkSC_err = self._getHistoValue(self.tkSCElectronHisto, pt, eta)
              #print ' pt, eta, tkSC, tkSC_err = ', pt, ' ', eta, ' ', tkSC, ' ', tkSC_err
              
              for point in self.idIsoScaleFactors[kindLep] : 
          
               #            eta       |      pt     | eff_data   stat  |  eff_mc   stat |      other nuisances
               #       -2.500  -2.000  10.000  20.000  0.358   0.009     0.286   0.002       0.094   0.048   0.071   0.127   -1      -1
          
                #
                # Procedure required by EGamma:
                # - electrons scale factors are provided in absolute eta bins
                #
                eta = abs(eta)
          
                if ( eta >= float(point[0]) and eta <= float(point[1]) and         # the "=" in both directions is only used by the overflow bin
                     pt  >= float(point[2]) and pt  <= float(point[3]) ) :         # in other cases the set is (min, max]
                    
                    data = float(point[4])
                    mc   = float(point[6])
          
                    sigma_data = float(point[5])
                    sigma_mc   = float(point[7])
                    
                    scaleFactor = data / mc
                    error_scaleFactor = math.sqrt((sigma_data / mc) * (sigma_data / mc) + (data / mc / mc * sigma_mc)*(data / mc / mc * sigma_mc))
                    
                    # systematic uncertainty
                    error_syst_scaleFactor = math.sqrt( float(point[8]) * float(point[8])   + 
                                                        float(point[9]) * float(point[9])   +
                                                        float(point[10]) * float(point[10]) +
                                                        float(point[11]) * float(point[11])  )
                    
                    error_syst_scaleFactor = error_syst_scaleFactor / mc
                    
                    #  idiso * reco scale factors 
                    if tkSC != 0 and wantOnlyRecoEff == 2:
                      # sum in quadrature the relative uncertainty
                      error_scaleFactor = scaleFactor * math.sqrt(error_scaleFactor/scaleFactor*error_scaleFactor/scaleFactor + tkSC_err/tkSC*tkSC_err/tkSC )
                      # now scale by the correction factor
                      scaleFactor *= tkSC
                      error_scaleFactor *= tkSC 
                      error_syst_scaleFactor *= tkSC 
          
                    return scaleFactor, error_scaleFactor, error_scaleFactor, error_syst_scaleFactor
          
              # default ... it should never happen!
              #print " default ele ???"
              return 1.0, 0.0, 0.0, 0.0
          
          
            elif kindLep == 'mu' :
              kindTight = ""
              if tight == 1 :
                kindTight = "muTight"
              else :
                kindTight = "muLoose"
               
               
              for point in self.isoScaleFactors[kindTight] : 
                iso_scaleFactor = 1
                iso_error_scaleFactor_up = 0
                iso_error_scaleFactor_do = 0
                
                if ( eta >= float(point[0]) and eta <= float(point[1]) and         # the "=" in both directions is only used by the overflow bin
                     pt  >= float(point[2]) and pt  <= float(point[3]) ) :         # in other cases the set is (min, max]
                    data = float(point[4])
                    mc   = float(point[7])
          
                    sigma_up_data = float(point[5])
                    sigma_up_mc   = float(point[8])
          
                    sigma_do_data = float(point[6])
                    sigma_do_mc   = float(point[9])
                    
                    iso_scaleFactor = data / mc
                    iso_error_scaleFactor_up = (data + sigma_up_data) / (mc - sigma_do_mc)  - iso_scaleFactor
                    iso_error_scaleFactor_do = iso_scaleFactor -   (data - sigma_do_data) / (mc + sigma_up_mc)  
             
                    break
              
             
              for point in self.idIsoScaleFactors[kindLep] : 
               #            eta       |      pt     | eff_data   stat up   stat down |  eff_mc   stat up   stat down  |      other nuisances
               #       -2.500  -2.000  10.000  20.000  0.358   0.009        0.009       0.286   0.002       0.009          0.094   0.048   0.071   0.127   -1      -1
          
                if ( eta >= float(point[0]) and eta <= float(point[1]) and         # the "=" in both directions is only used by the overflow bin
                     pt  >= float(point[2]) and pt  <= float(point[3]) ) :         # in other cases the set is (min, max]
                    data = float(point[4])
                    mc   = float(point[7])
          
                    sigma_up_data = float(point[5])
                    sigma_up_mc   = float(point[8])
          
                    sigma_do_data = float(point[6])
                    sigma_do_mc   = float(point[9])
                    
                    scaleFactor = data / mc
                    error_scaleFactor_up = (data + sigma_up_data) / (mc - sigma_do_mc)  - scaleFactor
                    error_scaleFactor_do = scaleFactor -   (data - sigma_do_data) / (mc + sigma_up_mc)  
                    
                    # multiply for isolation scale factor
                    #  -> sum in quadrature the relative uncertainties
                    error_scaleFactor_up = scaleFactor * iso_scaleFactor * math.sqrt(error_scaleFactor_up*error_scaleFactor_up/scaleFactor/scaleFactor +  iso_error_scaleFactor_up*iso_error_scaleFactor_up/iso_scaleFactor/iso_scaleFactor)
                    error_scaleFactor_do = scaleFactor * iso_scaleFactor * math.sqrt(error_scaleFactor_do*error_scaleFactor_do/scaleFactor/scaleFactor +  iso_error_scaleFactor_do*iso_error_scaleFactor_do/iso_scaleFactor/iso_scaleFactor)
                    scaleFactor *= iso_scaleFactor
                    
                    #                                                             no systematic uncertainty for the time being
                    return scaleFactor, error_scaleFactor_do, error_scaleFactor_up, 0.0
          
              # default ... it should never happen!
              #print " default mu ???"
              return 1.0, 0.0, 0.0, 0.0
          
          
            # not a lepton ... like some default value: and what can it be if not a lepton? ah ah 
            # --> it happens for default values -9999.
            return 1.0, 0.0, 0.0, 0.0
          
          # not a lepton ... like some default value: and what can it be if not a lepton? ah ah 
          # --> it happens for default values -9999.
          return 1.0, 0.0, 0.0, 0.0
        
        # reco scale factors only
        elif wantOnlyRecoEff == 1:
         
          scaleFactor = 1 
          error_scaleFactor = 0. 
          if kindLep == 'ele' :
            tkSC, tkSC_err = self._getHistoValue(self.tkSCElectronHisto, pt, eta)
            scaleFactor *= tkSC
            error_scaleFactor = tkSC_err 
          
          return scaleFactor, error_scaleFactor, error_scaleFactor, 0
             


   
   
    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)

        self.namesOldBranchesToBeModifiedVector = [
           'std_vector_lepton_recoW',
           'std_vector_lepton_recoW_Up',
           'std_vector_lepton_recoW_Down',                              
           'std_vector_lepton_idisoW',
           'std_vector_lepton_idisoW_Up',
           'std_vector_lepton_idisoW_Down',                              
           'std_vector_lepton_idisoW_Syst',                               
           'std_vector_lepton_idisoLooseW',
           'std_vector_lepton_idisoLooseW_Up',
           'std_vector_lepton_idisoLooseW_Down',                              
           'std_vector_lepton_idisoLooseW_Syst'                              
           ]
        
        self.clone(output,self.namesOldBranchesToBeModifiedVector)


        bvector_reco =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW',bvector_reco)
        bvector_reco_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW_Up',bvector_reco_Up)
        bvector_reco_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW_Down',bvector_reco_Down)

        bvector_idiso =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW',bvector_idiso)
        bvector_idiso_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Up',bvector_idiso_Up)
        bvector_idiso_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Down',bvector_idiso_Down)
        bvector_idiso_Syst =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Syst',bvector_idiso_Syst)
            
        bvector_idisoLoose =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoLooseW',bvector_idisoLoose)
        bvector_idisoLoose_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoLooseW_Up',bvector_idisoLoose_Up)
        bvector_idisoLoose_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoLooseW_Down',bvector_idisoLoose_Down)
        bvector_idisoLoose_Syst =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoLooseW_Syst',bvector_idisoLoose_Syst)
            
            
        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 
        savedentries = 0
                
        # avoid dots to go faster
        itree     = self.itree
        otree     = self.otree

        print '- Starting eventloop'
        step = 5000
        for i in xrange(nentries):
            itree.GetEntry(i)

            ## print event count
            if i > 0 and i%step == 0.:
              print i,'events processed.'

            bvector_reco.clear()
            bvector_reco_Up.clear()
            bvector_reco_Down.clear()

            bvector_idiso.clear()
            bvector_idiso_Up.clear()
            bvector_idiso_Down.clear()
            bvector_idiso_Syst.clear()

            bvector_idisoLoose.clear()
            bvector_idisoLoose_Up.clear()
            bvector_idisoLoose_Down.clear()
            bvector_idisoLoose_Syst.clear()

            for iLep in xrange(len(itree.std_vector_lepton_pt)) :
             
              pt = itree.std_vector_lepton_pt [iLep]
              eta = itree.std_vector_lepton_eta [iLep]
              flavour = itree.std_vector_lepton_flavour [iLep]
              
              kindLep = 'nonlep' # ele or mu
              if abs (flavour) == 11 : 
                kindLep = 'ele'
              elif abs (flavour) == 13 :
                kindLep = 'mu'
 
              #                                                              is tight lepton? 1=tight, 0=loose
              w, error_w_lo, error_w_up, error_w_syst = self._getWeight (kindLep, pt, eta, 1,                   0)
             
              bvector_idiso.push_back(w)
              bvector_idiso_Up.push_back(w+error_w_up)
              bvector_idiso_Down.push_back(w-error_w_lo)             
              bvector_idiso_Syst.push_back(w+error_w_syst)             

              w, error_w_lo, error_w_up, error_w_syst = self._getWeight (kindLep, pt, eta, 1,                   1)
             
              bvector_reco.push_back(w)
              bvector_reco_Up.push_back(w+error_w_up)
              bvector_reco_Down.push_back(w-error_w_lo)             


              loose_w, error_loose_w_lo, error_loose_w_up, error_loose_w_syst = self._getWeight (kindLep, pt, eta, 0,       0)
             
              bvector_idisoLoose.push_back(loose_w)
              bvector_idisoLoose_Up.push_back(loose_w+error_loose_w_up)
              bvector_idisoLoose_Down.push_back(loose_w-error_loose_w_lo)             
              bvector_idisoLoose_Syst.push_back(loose_w+error_loose_w_syst)             


            otree.Fill()
            savedentries+=1

        self.disconnect()
        print '- Eventloop completed'
        print '   Saved: ', savedentries, ' events'


