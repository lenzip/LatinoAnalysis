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

    
    def sortByPt(self, vectors):
      zipped=zip(*vectors)
      sortzipped=sorted(zipped, reverse=True)
      sortedvectors = zip(*sortzipped)

      for iv,sortedvector in enumerate(sortedvectors):
        #print "vector", iv
        size = len(sortedvector)
        for i in range(size):
          #print "old\t",(vectors[iv])[i],"\tnew\t",sortedvector[i] 
          (vectors[iv])[i] = sortedvector[i]  
        
      

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

        newbranchesSimpleVars = {}        
        newbranchesVecotor = {}
        for particle in particles:
          bvariable = numpy.zeros(1, dtype=numpy.int32) 
          newbranchesSimpleVars["n"+particle] = bvariable
          for quantity in intQuantities:
            bvector =  ROOT.std.vector(int) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector
          for quantity in floatQuantities:
            bvector =  ROOT.std.vector(float) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector
          for quantity in boolQuantities:
            bvector =  ROOT.std.vector(int) ()
            newbranchesVecotor[particle+"_"+quantity] = bvector

        for bname, bvariable in newbranchesVecotor.iteritems():
            print " bname   = ", bname
            print " bvariable = ", bvariable
            self.otree.Branch(bname,bvariable)
        for bname, bvariable in newbranchesSimpleVars.iteritems():
            print " bname   = ", bname
            print " bvariable = ", bvariable
            self.otree.Branch(bname,bvariable,bname+"/I")


  
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
            particle = ""
            if ( (abs(itree.GenPart_pdgId[ipart])==11 or abs(itree.GenPart_pdgId[ipart])==13) and itree.GenPart_status[ipart] == 1 ) or \
                 (abs(itree.GenPart_pdgId[ipart])==15 and itree.GenPart_statusFlags[ipart] >> 1 & 1 and itree.GenPart_statusFlags[ipart] >> 13 & 1) : # isDecayed and LastCOpy
              particle = "leptonGen"
            if itree.GenPart_pdgId[ipart] == 22 and itree.GenPart_status[ipart] == 1:
              particle = "photonGen"
            if (abs(itree.GenPart_pdgId[ipart])==12 or abs(itree.GenPart_pdgId[ipart])==14 or abs(itree.GenPart_pdgId[ipart])==16) and itree.GenPart_status[ipart] == 1:
              particle = "neutrinoGen"
            if particle != "":
              #print particle  
              newbranchesSimpleVars["n"+particle][0] += 1
              newbranchesVecotor[particle+"_pid"].push_back(itree.GenPart_pdgId[ipart])
              newbranchesVecotor[particle+"_pt"].push_back(itree.GenPart_pt[ipart])
              newbranchesVecotor[particle+"_eta"].push_back(itree.GenPart_eta[ipart])
              newbranchesVecotor[particle+"_phi"].push_back(itree.GenPart_phi[ipart])
              #print itree.GenPart_pdgId[ipart], itree.GenPart_genPartIdxMother[ipart], nGen
              motherid = -9999
              motherstatus = -9999
              if itree.GenPart_genPartIdxMother[ipart] > 0:
                motherid = itree.GenPart_pdgId[itree.GenPart_genPartIdxMother[ipart]] 
                motherstatus = itree.GenPart_status[itree.GenPart_genPartIdxMother[ipart]]
              newbranchesVecotor[particle+"_MotherPID"].push_back(motherid)
              newbranchesVecotor[particle+"_MotherStatus"].push_back(motherstatus)
              newbranchesVecotor[particle+"_fromHardProcess"].push_back(int(itree.GenPart_statusFlags[ipart] >> 8 & 1))
              newbranchesVecotor[particle+"_isDirectHadronDecayProduct"].push_back(int(itree.GenPart_statusFlags[ipart] >> 6 & 1))
              newbranchesVecotor[particle+"_isDirectPromptTauDecayProduct"].push_back(int(itree.GenPart_statusFlags[ipart] >> 5 & 1))
              newbranchesVecotor[particle+"_isPrompt"].push_back(int(itree.GenPart_statusFlags[ipart] & 1))
              newbranchesVecotor[particle+"_isTauDecayProduct"].push_back(int(itree.GenPart_statusFlags[ipart] >> 3 & 1))
          for particle in particles:
            vectorsToOrder = []
            vectorsToOrder.append(newbranchesVecotor[particle+"_pt"])
            #print "old",newbranchesVecotor[particle+"_pt"][0]
            for quantity in allQuantities:
              if quantity != "pt":
                vectorsToOrder.append(newbranchesVecotor[particle+"_"+quantity])
            self.sortByPt(vectorsToOrder)
            #print "new",newbranchesVecotor[particle+"_pt"][0]
          otree.Fill()
          for particle in particles:
            for quantity in allQuantities:
              newbranchesVecotor[particle+"_"+quantity].clear()
          for particle in particles:
            newbranchesSimpleVars["n"+particle][0]=0
        self.disconnect()
        print '- Eventloop completed'

