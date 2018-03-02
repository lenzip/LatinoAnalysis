#
#
#
#
#   Quark TOP pT
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

class HiggsGenVar(TreeCloner):
    def __init__(self):
       pass

    def help(self):
        return ''' Calculate higgs gen variables '''

    def addOptions(self,parser):
        pass

    def checkOptions(self,opts):
        pass

    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.newbranches = ['higgsGenPt',
                            'higgsGenEta',    
                            'higgsGenPhi',
                            'higgsGenMass',
                            ]

        # does that work so easily and give new variable itree and otree?
        self.connect(tree,input)
 

        self.clone(output,self.newbranches)

                
        newbranchesVecotor = {}
        for bname in self.newbranches:
          #bvector =  ROOT.std.vector(float) ()
          bvariable = numpy.array([-9999.], dtype=numpy.float32)
          newbranchesVecotor[bname] = bvariable


        for bname, bvariable in newbranchesVecotor.iteritems():
            print " bname   = ", bname
            print " bvariable = ", bvariable
            self.otree.Branch(bname,bvariable,bname+'/F')


  
        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 

        # input tree and output tree
        itree     = self.itree
        otree     = self.otree
        numTOP=0.
        numAntiTOP=0.
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
            isLastCopy=int("010000000000000", 2)
            if itree.GenPart_pdgId[ipart] == 25 and (itree.GenPart_statusFlags[ipart] & isLastCopy): #isLastCopy
              newbranchesVecotor["higgsGenPt"][0]   = itree.GenPart_pt[ipart]
              newbranchesVecotor["higgsGenEta"][0]  = itree.GenPart_eta[ipart]
              newbranchesVecotor["higgsGenPhi"][0]  = itree.GenPart_phi[ipart]
              newbranchesVecotor["higgsGenMass"][0] = itree.GenPart_mass[ipart]

          #for bname, bvector in newbranchesVecotor.iteritems():
          #    bvector.clear()
                
          #print "----------------"
          
          otree.Fill()

        self.disconnect()
        print '- Eventloop completed'

