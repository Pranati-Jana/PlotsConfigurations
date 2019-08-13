import os
import copy
from LatinoAnalysis.Tools.commonTools import getSampleFiles, getBaseW, addSampleWeight

def nanoGetSampleFiles(inputDir, sample):
    try:
        if _samples_noload:
            return []
    except NameError:
        pass

    return getSampleFiles(inputDir, sample, True, 'nanoLatino_')

# samples

try:
    len(samples)
except NameError:
    import collections
    samples = collections.OrderedDict()

################################################
############### Fiducial bins ##################
################################################

#pthBinning = ['0', '10', '15', '20', '30', '45', '60', '80', '100', '120', '155', '200', '260', '350', 'inf']
pthBinning = ['0', '20', '45', '80', '120', '200', '350', 'inf']
pthBins = []
for ibin in range(len(pthBinning) - 1):
    low, high = pthBinning[ibin:ibin + 2]
    if high == 'inf':
      pthBins.append('GT%s' % low)
    else:
      pthBins.append('%s_%s' % (low, high))

#yhBinning = [0., 0.15, 0.3, 0.6, 0.9, 1.2, 2.5, 10.]
njetBinning = ['0', '1', '2', '3', 'GE4']

################################################
################# SKIMS ########################
################################################

mcProduction = 'Autumn18_102X_nAODv5_Full2018v5'

dataReco = 'Run2018_102X_nAODv5_Full2018v5'

mcSteps = 'MCl1loose2018v5__MCCorr2018v5__l2loose__l2tightOR2018v5{var}' #FIXME, no wwSel skim yet

fakeSteps = 'DATAl1loose2016__l2loose__fakeW__wwSel'

dataSteps = 'DATAl1loose2018v5__l2loose__l2tightOR2018v5' #FIXME, no wwSel skim yet

##############################################
###### Tree base directory for the site ######
##############################################

SITE=os.uname()[1]
if    'iihe' in SITE:
  treeBaseDir = '/pnfs/iihe/cms/store/user/xjanssen/HWW2015'
elif  'cern' in SITE:
  treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano'

def makeMCDirectory(var=''):
    if var:
        return os.path.join(treeBaseDir, mcProduction, mcSteps.format(var='__' + var))
    else:
        return os.path.join(treeBaseDir, mcProduction, mcSteps.format(var=''))

mcDirectory = makeMCDirectory()
fakeDirectory = os.path.join(treeBaseDir, dataReco, fakeSteps)
dataDirectory = os.path.join(treeBaseDir, dataReco, dataSteps)

################################################
############ DATA DECLARATION ##################
################################################

DataRun = [
            ['A','Run2018A-Nano1June2019-v1'] ,
            ['B','Run2018B-Nano1June2019-v1'] ,
            ['C','Run2018C-Nano1June2019-v1'] ,
            ['D','Run2018D-Nano1June2019_ver2-v1'] ,
          ]

DataSets = ['MuonEG','DoubleMuon','SingleMuon','EGamma']

DataTrig = {
            'MuonEG'         : 'Trigger_ElMu' ,
            'DoubleMuon'     : '!Trigger_ElMu && Trigger_dblMu' ,
            'SingleMuon'     : '!Trigger_ElMu && !Trigger_dblMu && Trigger_sngMu' ,
            'EGamma'         : '!Trigger_ElMu && !Trigger_dblMu && !Trigger_sngMu && (Trigger_sngEl || Trigger_dblEl)' ,
           }


#########################################
############ MC COMMON ##################
#########################################

# SFweight does not include btag weights
mcCommonWeightNoMatch = 'XSWeight*SFweight*METFilter_MC'
mcCommonWeight = 'XSWeight*SFweight*1*METFilter_MC' #'XSWeight*SFweight*PromptGenLepMatch2l*METFilter_MC' #FIXME include genMatch when we have data fakes

###########################################
#############  BACKGROUNDS  ###############
###########################################

###### DY #######

ptllDYW_NLO = '(((0.623108 + 0.0722934*gen_ptll - 0.00364918*gen_ptll*gen_ptll + 6.97227e-05*gen_ptll*gen_ptll*gen_ptll - 4.52903e-07*gen_ptll*gen_ptll*gen_ptll*gen_ptll)*(gen_ptll<45)*(gen_ptll>0) + 1*(gen_ptll>=45))*(abs(gen_mll-90)<3) + (abs(gen_mll-90)>3))'
ptllDYW_LO = '((0.632927+0.0456956*gen_ptll-0.00154485*gen_ptll*gen_ptll+2.64397e-05*gen_ptll*gen_ptll*gen_ptll-2.19374e-07*gen_ptll*gen_ptll*gen_ptll*gen_ptll+6.99751e-10*gen_ptll*gen_ptll*gen_ptll*gen_ptll*gen_ptll)*(gen_ptll>0)*(gen_ptll<100)+(1.41713-0.00165342*gen_ptll)*(gen_ptll>=100)*(gen_ptll<300)+1*(gen_ptll>=300))'

files = nanoGetSampleFiles(mcDirectory, 'DYJetsToLL_M-50-LO') + \
    nanoGetSampleFiles(mcDirectory, 'DYJetsToLL_M-10to50-LO') #FIXME swithch to NLO M50 sample when abailable

samples['DY'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 4,
}
addSampleWeight(samples,'DY','DYJetsToTT_MuEle_M-50',ptllDYW_LO)  #FIXME... and remember to change this as well to NLO
addSampleWeight(samples,'DY','DYJetsToLL_M-10to50-LO',ptllDYW_LO)

###### Top #######

files = nanoGetSampleFiles(mcDirectory, 'TTTo2L2Nu') + \
    nanoGetSampleFiles(mcDirectory, 'ST_s-channel_ext1') + \
    nanoGetSampleFiles(mcDirectory, 'ST_t-channel_antitop') + \
    nanoGetSampleFiles(mcDirectory, 'ST_t-channel_top') + \
    nanoGetSampleFiles(mcDirectory, 'ST_tW_antitop_ext1') + \
    nanoGetSampleFiles(mcDirectory, 'ST_tW_top_ext1')

samples['top'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 1,
    #'EventsPerJob': 100000
}

Top_pTrw = "(TMath::Sqrt( TMath::Exp(0.0615-0.0005*topGenPt) * TMath::Exp(0.0615-0.0005*antitopGenPt) ) )"

addSampleWeight(samples,'top','TTTo2L2Nu',Top_pTrw)

###### WW ########

samples['WW'] = {
    'name': nanoGetSampleFiles(mcDirectory, 'WWTo2L2Nu'),
    #'weight': mcCommonWeight + '*nllW', # temporary
    'weight': mcCommonWeight + '*nllWOTF', # temporary
    'FilesPerJob': 8
}

samples['WWewk'] = {
    'name': nanoGetSampleFiles(mcDirectory, 'WpWmJJ_EWK'),
    'weight': mcCommonWeight + '*(Sum$(abs(GenPart_pdgId)==6 || GenPart_pdgId==25)==0)', #filter tops and Higgs
    'FilesPerJob': 4
}

samples['ggWW'] = {
    'name': nanoGetSampleFiles(mcDirectory, 'GluGluToWWToENEN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToENMN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToENTN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToMNEN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToMNMN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToMNTN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToTNEN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToTNMN') \
          + nanoGetSampleFiles(mcDirectory, 'GluGluToWWToTNTN'),
    'weight': mcCommonWeight + '*1.53/1.4', # updating k-factor
    'FilesPerJob': 4
}

######## Vg ########

files = nanoGetSampleFiles(mcDirectory, 'Wg_MADGRAPHMLM') + \
    nanoGetSampleFiles(mcDirectory, 'Zg')

samples['Vg'] = {
    'name': files,
    'weight': mcCommonWeightNoMatch + '*!(Gen_ZGstar_mass > 0 && Gen_ZGstar_MomId == 22)',
    'FilesPerJob': 4
}

######## VgS ########

files = nanoGetSampleFiles(mcDirectory, 'Wg_MADGRAPHMLM') + \
    nanoGetSampleFiles(mcDirectory, 'Zg') + \
    nanoGetSampleFiles(mcDirectory, 'WZTo3LNu_mllmin01')

samples['VgS'] = {
    'name': files,
    'weight': mcCommonWeight,# + ' * (gstarLow * 0.94 + gstarHigh * 1.14)',
    'FilesPerJob': 4,
    'subsamples': {
      'L': 'gstarLow',
      'H': 'gstarHigh'
    }
}
addSampleWeight(samples, 'VgS', 'Wg_MADGRAPHMLM', '(Gen_ZGstar_mass > 0 && Gen_ZGstar_MomId == 22 && Gen_ZGstar_mass < 0.1)')
addSampleWeight(samples, 'VgS', 'Zg', '(Gen_ZGstar_mass > 0 && Gen_ZGstar_MomId == 22)')
addSampleWeight(samples, 'VgS', 'WZTo3LNu_mllmin01', '(Gen_ZGstar_mass > 0.1 || Gen_ZGstar_mass < 0)')

############ VZ ############

files = nanoGetSampleFiles(mcDirectory, 'ZZTo2L2Nu_ext1') + \
    nanoGetSampleFiles(mcDirectory, 'ZZTo2L2Q') + \
    nanoGetSampleFiles(mcDirectory, 'ZZTo4L_ext1') + \
    nanoGetSampleFiles(mcDirectory, 'WZTo2L2Q')

samples['VZ'] = {
    'name': files,
    'weight': mcCommonWeight + '*1.11',
    'FilesPerJob': 4
}

########## VVV #########

files = nanoGetSampleFiles(mcDirectory, 'ZZZ') + \
    nanoGetSampleFiles(mcDirectory, 'WZZ') + \
    nanoGetSampleFiles(mcDirectory, 'WWZ') + \
    nanoGetSampleFiles(mcDirectory, 'WWW')
#+ nanoGetSampleFiles(mcDirectory, 'WWG'), #should this be included? or is it already taken into account in the WW sample?

samples['VVV'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 4
}

###########################################
#############   SIGNALS  ##################
###########################################

signals = []

#### ggH -> WW

samples['ggH_hww'] = {
    #'name': nanoGetSampleFiles(mcDirectory, 'GluGluHToWWTo2L2NuPowheg_M125_PrivateNano'),
    'name': nanoGetSampleFiles(mcDirectory, 'GluGluHToWWTo2L2NuPowheg_M125'),
    'weight': [mcCommonWeight, {'class': 'Weight2MINLO', 'args': '%s/src/LatinoAnalysis/Gardener/python/data/powheg2minlo/NNLOPS_reweight.root' % os.getenv('CMSSW_BASE')}],
    'FilesPerJob': 4,
    'linesToAdd': ['.L %s/work/PlotsConfigurations/Configurations/Differential/weight2MINLO.cc+' % os.getenv('CMSSW_BASE')]
}

signals.append('ggH_hww')

############ VBF H->WW ############
samples['qqH_hww'] = {
    'name': nanoGetSampleFiles(mcDirectory, 'VBFHToWWTo2L2NuPowheg_M125'),
    'weight': mcCommonWeight,
    'FilesPerJob': 4
}

signals.append('qqH_hww')

############ ZH H->WW ############

samples['ZH_hww'] = {
    'name':   nanoGetSampleFiles(mcDirectory, 'HZJ_HToWW_M125'),
    'weight': mcCommonWeight,
    'FilesPerJob': 4
}

signals.append('ZH_hww')

samples['ggZH_hww'] = {
    'name':   nanoGetSampleFiles(mcDirectory, 'GluGluZH_HToWW_M125'),
    'weight': mcCommonWeight,
    'FilesPerJob': 4
}

signals.append('ggZH_hww')

############ WH H->WW ############

samples['WH_hww'] = {
    'name':   nanoGetSampleFiles(mcDirectory, 'HWplusJ_HToWW_M125') + nanoGetSampleFiles(mcDirectory, 'HWminusJ_HToWW_M125'),
    'weight': mcCommonWeight,
    'FilesPerJob': 4
}

signals.append('WH_hww')

############ ttH ############

samples['ttH_hww'] = {
    'name':   nanoGetSampleFiles(mcDirectory, 'ttHToNonbb_M125'),
    'weight': mcCommonWeight,
    'FilesPerJob': 1
}

signals.append('ttH_hww')

############ bbH ############
#FIXME Missing samples

############ H->TauTau ############

#samples['ggH_htt'] = {
#    'name': nanoGetSampleFiles(mcDirectory, 'GluGluHToTauTau_M125'),
#    'weight': mcCommonWeight,
#    'FilesPerJob': 4
#}
#
#signals.append('ggH_htt')
#
#samples['qqH_htt'] = {
#    'name': nanoGetSampleFiles(mcDirectory, 'VBFHToTauTau_M125'),
#    'weight': mcCommonWeight,
#    'FilesPerJob': 4
#}
#
#signals.append('qqH_htt')
#
#samples['ZH_htt'] = {
#    'name': nanoGetSampleFiles(mcDirectory, 'HZJ_HToTauTau_M125'),
#    'weight': mcCommonWeight,
#    'FilesPerJob': 4
#}
#
#signals.append('ZH_htt')
#
#samples['WH_htt'] = {
#    'name':  nanoGetSampleFiles(mcDirectory, 'HWplusJ_HToTauTau_M125') + nanoGetSampleFiles(mcDirectory, 'HWminusJ_HToTauTau_M125'),
#    'weight': mcCommonWeight,
#    'FilesPerJob': 4
#}
#
#signals.append('WH_htt')

for sname in signals:
  sample = samples[sname]
  sample['subsamples'] = {}

  for pth in pthBins:
    binName = 'PTH_%s' % pth
    if pth.startswith('GT'):
      cut = 'higgsGenPt > %s' % pth[2:]
    else:
      cut = 'higgsGenPt > %s && higgsGenPt < %s' % tuple(pth.split('_'))

    sample['subsamples'][binName] = cut

  for nj in njetBinning:
    binName = 'NJ_%s' % nj
    if nj.startswith('GE'):
      cut = 'nCleanGenJet >= %s' % nj[2:]
    else:
      cut = 'nCleanGenJet == %s' % nj

    sample['subsamples'][binName] = cut

#  for flabel, fidcut in [('fid', 'fiducial'), ('nonfid', '!fiducial')]:
#    for pth in pthBins:
#      binName = '%s_PTH_%s' % (flabel, pth)
#      if pth.startswith('GT'):
#        cut = '%s && genPth > %s' % (fidcut, pth[2:])
#      else:
#        cut = '%s && genPth > %s && genPth < %s' % ((fidcut,) + tuple(pth.split('_')))
#  
#      sample['subsamples'][binName] = cut
#  
#    for nj in njetBinning:
#      binName = '%s_NJ_%s' % (flabel, nj)
#      if nj.startswith('GE'):
#        cut = '%s && nCleanGenJet >= %s' % (fidcut, nj[2:])
#      else:
#        cut = '%s && nCleanGenJet == %s' % (fidcut, nj)
#  
#      sample['subsamples'][binName] = cut

###########################################
################## FAKE ###################
###########################################

#samples['Fake'] = {
#  'name': [],
#  'weight': 'METFilter_DATA*fakeW',
#  'weights': [],
#  'isData': ['all'],
#  'FilesPerJob': 50
#}
#
#for _, sd in DataRun:
#  for pd in DataSets:
#    files = nanoGetSampleFiles(fakeDirectory, pd + '_' + sd)
#    samples['Fake']['name'].extend(files)
#    samples['Fake']['weights'].extend([DataTrig[pd]] * len(files))
#
#samples['Fake']['subsamples'] = {
#  'em': 'abs(Lepton_pdgId[0]) == 11',
#  'me': 'abs(Lepton_pdgId[0]) == 13'
#}

samples['MCFake']  = {  'name'   :   nanoGetSampleFiles(mcDirectory, 'WJetsToLNu-LO')
                                  + nanoGetSampleFiles(mcDirectory,'TTToSemiLeptonic'),
                       'weight' : mcCommonWeight,
                       'FilesPerJob': 3,
                    }

samples['MCFake']['subsamples'] = {
  'em': 'abs(Lepton_pdgId[0]) == 11',
  'me': 'abs(Lepton_pdgId[0]) == 13'
}



###########################################
################## DATA ###################
###########################################

samples['DATA'] = {
  'name': [],
  'weight': 'METFilter_DATA*LepWPCut',
  'weights': [],
  'isData': ['all'],
  'FilesPerJob': 100
}

for _, sd in DataRun:
  for pd in DataSets:
    files = nanoGetSampleFiles(dataDirectory, pd + '_' + sd)
    samples['DATA']['name'].extend(files)
    samples['DATA']['weights'].extend([DataTrig[pd]] * len(files))
