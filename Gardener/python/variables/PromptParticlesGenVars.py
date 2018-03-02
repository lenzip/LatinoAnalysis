#
#
#
#
#  Stable particles Gen Vars 
#                                
#
#


from LatinoAnalysis.Gardener.gardening import TreeCloner
import numpy
import ROOT
import sys
import optparse
import re
import warnings
import os.path
from array import array;
import math
import copy

class PromptParticlesGenVars(TreeCloner):
    def __init__(self):
       pass

    def help(self):
        return ''' Calculate gen variables '''

    def addOptions(self,parser):
        pass

    def checkOptions(self,opts):
        pass

    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        particles=["leptonGen", "neutrinoGen", "photonGen"]
        intQuantities   = ["MotherPID", "MotherStatus"] 
        floatQuantities = ["eta", "phi", "pt", "pid"]
        boolQuantities  = ["fromHardProcess", "isDirectHadronDecayProduct", "isDirectPromptTauDecayProduct", "isPrompt", "isTauDecayProduct"]
        allQuantities = copy.deepcopy(intQuantities)
        allQuantities.extend(floatQuantities)
        allQuantities.extend(boolQuantities)

        self.newbranches = []
        for particle in particles:
          self.newbranches.append("n"+particle)
          for quantity in allQuantities:
            self.newbranches.append(particle+"_"+quantity)

        # does that work so easily and give new variable itree and otree?
        self.connect(tree,input)
 

        self.clone(output,self.newbranches)

                
        newbranchesVecotor = {}
        for particle in particles:
          bvariable = numpy.zeros(1, dtype=numpy.int32) 
          newbranchesVecotor["n"+] = bvariable
          for quantity in intQuantities:
            bvector =  ROOT.std.vector(int) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector
          for quantity in floatQuantities:
            bvector =  ROOT.std.vector(float) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector
          for quantity in boolQuantities:
            bvector =  ROOT.std.vector(bool) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector

        for bname, bvariable in newbranchesVecotor.iteritems():
            print " bname   = ", bname
            print " bvariable = ", bvariable
            self.otree.Branch(bname,bvariable,bname+'/F')


  
        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 

        # input tree and output tree
        itree     = self.itree
        otree     = self.otree
        #----------------------------------------------------------------------------------------------------
        print '- Starting eventloop'
        step = 5000

        #for i in xrange(5000):
        for i in xrange(nentries):

          itree.GetEntry(i)

          if i > 0 and i%step == 0.:
              print i,'events processed.'

          nGen = itree.nGenPart
          for ipart in range(nGen):
            isDecayedAndLastCopy=int("0100000000001000", 2)
            particle = ""
            if ( (abs(itree.GenPart_pdgId[ipart])==11 or abs(itree.GenPart_pdgId[ipart])==13) and itree.GenPart_status == 1 ) or \
                 (abs(itree.GenPart_pdgId[ipart])==15 and (itree.GenPart_statusFlags[ipart] & isDecayedAndLastCopy) == isDecayedAndLastCopy) :
              particle = "leptonGen"
            if itree.GenPart_pdgId[ipart]) == 22 and itree.GenPart_status == 1:
              particle = "photonGen"
            if (abs(itree.GenPart_pdgId[ipart])=12 or abs(itree.GenPart_pdgId[ipart])==14 or abs(itree.GenPart_pdgId[ipart])==16) and itree.GenPart_status == 1:
              particle = "neutrinoGen"
            if particle != "":
              newbranchesVecotor["n"+particle][0] += 1
              newbranchesVecotor[particle+"_pid"].push_back(itree.GenPart_pdgId[ipart])
              newbranchesVecotor[particle+"_pt"].push_back(itree.GenPart_pt[ipart])
              newbranchesVecotor[particle+"_eta"].push_back(itree.GenPart_eta[ipart])
              newbranchesVecotor[particle+"_phi"].push_back(itree.GenPart_phi[ipart])
              newbranchesVecotor[particle+"_MotherPID"].push_back(itree.GenPart_pdgId[itree.GenPart_genPartIdxMother[ipart]])
              newbranchesVecotor[particle+"_MotherStatus"].push_back(itree.GenPart_status[itree.GenPart_genPartIdxMother[ipart]])
              isFromHardProcess = int("000000100000000", 2)
              newbranchesVecotor[particle+"_fromHardProcess"].push_back(bool(itree.GenPart_statusFlags[ipart] & isFromHardProcess))
              isDirectHadronDecayProduct = int("000000010000000", 2)  
              newbranchesVecotor[particle+"_isDirectHadronDecayProduct"].push_back(bool(itree.GenPart_statusFlags[ipart] & isDirectHadronDecayProduct))
              isDirectPromptTauDecayProduct = int("000000001000000", 2)
              newbranchesVecotor[particle+"_isDirectPromptTauDecayProduct"].push_back(bool(itree.GenPart_statusFlags[ipart] & isDirectPromptTauDecayProduct))
              isPrompt = int("000000000000001", 2)
              newbranchesVecotor[particle+"_isPrompt"].push_back(bool(itree.GenPart_statusFlags[ipart] & isPrompt))
              isTauDecayProduct = int("000000000001000", 2)
              newbranchesVecotor[particle+"_isTauDecayProduct"].push_back(bool(itree.GenPart_statusFlags[ipart] & isTauDecayProduct))
                

          otree.Fill()
          for bname, bvector in newbranchesVecotor.iteritems():
            bvector.clear()
          for particle in particles:
            newbranchesVecotor["n"+particle][0]=0
        self.disconnect()
        print '- Eventloop completed'

