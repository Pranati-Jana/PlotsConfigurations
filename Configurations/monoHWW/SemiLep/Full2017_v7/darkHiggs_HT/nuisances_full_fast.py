# nuisances

#nuisances = {}

# name of samples here must match keys in samples.py

# imported from samples.py:
# samples, treeBaseDir, mcProduction, mcSteps
# imported from cuts.py
# cuts

from LatinoAnalysis.Tools.commonTools import getSampleFiles, getBaseW, addSampleWeight

def nanoGetSampleFiles(inputDir, Sample):
    return getSampleFiles(inputDir, Sample, False, 'nanoLatino_')

try:
    mc = [skey for skey in samples if skey != 'DATA' and not skey.startswith('FAKE')]
except NameError:
    mc = []
    cuts = {}
    nuisances = {}
    def makeMCDirectory(x=''):
        return ''

from LatinoAnalysis.Tools.HiggsXSection import HiggsXSection
HiggsXS = HiggsXSection()

#mu_cuts = [ c for c in cuts if 'MuCh' in c] 
#el_cuts = [ c for c in cuts if 'ElCh' in c] 

mc_VBS = ['WWewk', 'WZqcd', 'WZewk', 'ZZ']
#for skey in mc:
#    print(skey, samples[skey])
#    if '_smp' in samples[skey]['name'][0]: mc_VBS.append(skey)
#mc_VBS   = [skey for skey in mc if '_smp' in samples[skey]['name'][0]]
mc_noVBS = [skey for skey in mc if not skey in mc_VBS]
mc_noTop = [skey for skey in mc if not 'top' in skey]
mc_top   = [skey for skey in mc if 'top' in skey]
#print('VBS samples:' + str(mc_VBS))

def makeSuffixVar(name, fix, samples_l, tag='', cuts=None, folder_fix=None, as_lnN=False):
    nuis = {}
    nuis['name']    = name
    nuis['kind']    = 'suffix'
    nuis['type']    = 'shape'
    nuis['mapUp']   = fix+'up'
    nuis['mapDown'] = fix+'do'
    if as_lnN: nuis['AsLnN'] = '1'
    if not cuts is None: nuis['cuts'] = cuts
    f_fix = fix
    if not folder_fix is None: f_fix = folder_fix
    
    #if not 'JES' in fix:    
    #    nuis_VBS = copy.deepcopy(nuis)
    #    nuis_VBS['samples']    = dict((skey, ['1.', '1.']) for skey in mc_VBS)
    #    nuis_VBS['folderUp']   = makeMCDirectory(var=f_fix+'up', base='/eos/cms/store/group/phys_smp/VJets_NLO_VBSanalyses')    
    #    nuis_VBS['folderDown'] = makeMCDirectory(var=f_fix+'do', base='/eos/cms/store/group/phys_smp/VJets_NLO_VBSanalyses')    
    #    nuisances[name+'_VBS'] = copy.deepcopy(nuis_VBS)

    detect_VBS = False
    if samples_l[0] in mc_VBS: detect_VBS = True

    nuis_noVBS = copy.deepcopy(nuis)
    nuis_noVBS['samples']    = dict((skey, ['1.', '1.']) for skey in samples_l)
    if not detect_VBS:
        nuis_noVBS['folderUp']   = makeMCDirectory(var=f_fix+'up')    
        nuis_noVBS['folderDown'] = makeMCDirectory(var=f_fix+'do')
    else:
        if 'iihe' in SITE:
            base = '/pnfs/iihe/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano_smp'
        else:
            base='/eos/cms/store/group/phys_smp/VJets_NLO_VBSanalyses'
        nuis_noVBS['folderUp']   = makeMCDirectory(var=f_fix+'up', base=base)    
        nuis_noVBS['folderDown'] = makeMCDirectory(var=f_fix+'do', base=base)
    nuisances[name+tag] = copy.deepcopy(nuis_noVBS)
    

################################ EXPERIMENTAL UNCERTAINTIES  #################################

#### Luminosity

nuisances['lumi'] = {
   'name': 'lumi_13TeV_2017',
   'type': 'lnN',
   'samples': dict((skey, '1.023') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
}

#nuisances['lumi_Uncorrelated'] = {
#    'name': 'lumi_13TeV_2017',
#    'type': 'lnN',
#    'samples': dict((skey, '1.02') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_XYFact'] = {
#    'name': 'lumi_13TeV_XYFact',
#    'type': 'lnN',
#    'samples': dict((skey, '1.008') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_LScale'] = {
#    'name': 'lumi_13TeV_LSCale',
#    'type': 'lnN',
#    'samples': dict((skey, '1.003') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_BBDefl'] = {
#    'name': 'lumi_13TeV_BBDefl',
#    'type': 'lnN',
#    'samples': dict((skey, '1.004') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_DynBeta'] = {
#    'name': 'lumi_13TeV_DynBeta',
#    'type': 'lnN',
#    'samples': dict((skey, '1.005') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_CurrCalib'] = {
#    'name': 'lumi_13TeV_CurrCalib',
#    'type': 'lnN',
#    'samples': dict((skey, '1.003') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}
#
#nuisances['lumi_Ghosts'] = {
#    'name': 'lumi_13TeV_Ghosts',
#    'type': 'lnN',
#    'samples': dict((skey, '1.001') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top'])
#}

#### FAKES E_T, stat up/do var and 30% syst

nuisances['fake_syst'] = {
    'name': 'CMS_fake_syst_2017_18',
    'type': 'lnN',
    'samples': {
        'FAKE': '1.3',
    },
}

el_et = El_jetEt
mu_et = Mu_jetEt
for syst in ['El', 'statEl', 'Mu', 'statMu']:
    name_tag = ''
    #if 'stat' in syst: name_tag += '_stat_2017'
    #if 'El' in syst: name_tag += '_e_2017_18'
    #else: name_tag += '_m_2017_18'

    if 'El' in syst: name_tag += '_e'
    else: name_tag += '_m'
    if 'stat' in syst: name_tag += '_stat_2017'
    else: name_tag += '_recoil_2017_18'

    nuisances['fakeW_'+syst] = {
        'name': 'CMS_fake'+name_tag,
        'kind': 'weight',
        'type': 'shape',
        'samples': {
            'FAKE': ['FW_mu'+str(mu_et)+'_el'+str(el_et)+'_'+syst+'Up[0]', 'FW_mu'+str(mu_et)+'_el'+str(el_et)+'_'+syst+'Down[0]'],
            #'FAKE': ['FW_mu20_el35_'+syst+'Up[0]', 'FW_mu20_el35_'+syst+'Down[0]'],
        },
    }

##### B-tagger

for shift in ['jes', 'lf', 'hf', 'hfstats1', 'hfstats2', 'lfstats1', 'lfstats2', 'cferr1', 'cferr2']:
    btag_syst = ['(btagSF%sup)/(btagSF)' % shift, '(btagSF%sdown)/(btagSF)' % shift]

    name = 'CMS_btag_%s' % shift
    if 'stats' in shift:
        name += '_2017'

    nuisances['btag_shape_%s' % shift] = {
        'name': name,
        'kind': 'weight',
        'type': 'shape',
        'samples': dict((skey, btag_syst) for skey in mc),
    }

##### Trigger Efficiency

#trig_syst = ['((TriggerEffWeight_1l_u)/(TriggerEffWeight_1l))*(TriggerEffWeight_1l>0.02) + (TriggerEffWeight_1l<=0.02)', '(TriggerEffWeight_1l_d)/(TriggerEffWeight_1l)']
trig_syst = ['((TriggerEffWeight_1l_fixed_u)/(TriggerEffWeight_1l_fixed))*(TriggerEffWeight_1l_fixed>0.02) + (TriggerEffWeight_1l_fixed<=0.02)', '(TriggerEffWeight_1l_fixed_d)/(TriggerEffWeight_1l_fixed)']

nuisances['trigg'] = {
    'name': 'CMS_eff_hwwtrigger_2017',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, trig_syst) for skey in mc)
}

prefire_syst = ['PrefireWeight_Up/PrefireWeight', 'PrefireWeight_Down/PrefireWeight']

nuisances['prefire'] = {
    'name': 'CMS_eff_prefiring_2017',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, prefire_syst) for skey in mc),
}

##### Electron Efficiency and energy scale

nuisances['eff_e'] = {
    'name': 'CMS_eff_e_2017',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightEleUp[0]', 'SFweightEleDown[0]']) for skey in mc)
}

#makeSuffixVar('CMS_scale_e_2017', 'ElepT')
makeSuffixVar('CMS_scale_e_2017', 'ElepT', mc_noVBS, as_lnN=True)

##### Muon Efficiency and energy scale

nuisances['eff_m'] = {
    'name': 'CMS_eff_m_2017',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightMuUp[0]', 'SFweightMuDown[0]']) for skey in mc)
}

#makeSuffixVar('CMS_scale_m_2017', 'MupT')
makeSuffixVar('CMS_scale_m_2017', 'MupT', mc_noVBS, as_lnN=True)

##### Jet energy scale
#jes_systs = ['JESAbsolute','JESAbsolute_2017','JESBBEC1','JESBBEC1_2017','JESEC2','JESEC2_2017','JESFlavorQCD','JESHF','JESHF_2017','JESRelativeBal','JESRelativeSample_2017']
#jes_systs = [
#    'JESAbsolute','JESAbsolute_2017','JESBBEC1','JESBBEC1_2017','JESEC2',
#    'JESEC2_2017','JESFlavorQCD','JESHF','JESHF_2017','JESRelativeBal',
#    'JESRelativeSample_2017'
#]
#
#for syst in jes_systs:
#    makeSuffixVar('CMS_scale_'+syst, syst, folder_fix='JES')

## DYlow no longer exists and VgS and VBF are small anyway
#nuisances['CMS_scale_JES_bogus_2017'] = {
#    'name': 'CMS_scale_JES_b_2017',
#    'type': 'lnN',
#    'samples': {
#        'DYlow': '1.06',
#        'VgS_L': '1.04',
#        'VBF-V': '1.02',
#    },
#}

mc_forJES = [skey for skey in mc_noVBS if not skey in ['']]

# Unified
makeSuffixVar('CMS_scale_JES_2017',     'JES', mc_forJES, as_lnN=True)
#makeSuffixVar('CMS_scale_JES_2017',     'JES', mc_forJES, as_lnN=False)

## top/rest split
#makeSuffixVar('CMS_scale_JES_2017',     'JES', [skey for skey in mc_forJES if skey in mc_noTop], as_lnN=True)
#makeSuffixVar('CMS_scale_JES_top_2017', 'JES', [skey for skey in mc_forJES if skey in mc_top] )


##### Jet energy resolution

#makeSuffixVar('CMS_res_j', 'JER')
mc_forJER = [skey for skey in mc_noVBS if not skey in ['WjetsNLO', 'Wjets', 'top']]
makeSuffixVar('CMS_res_j_2017', 'JER', mc_forJER, as_lnN=True)

##### MET energy scale

mc_forMET = [skey for skey in mc_noVBS if not skey in ['WjetsNLO', 'Wjets', 'top']]
#makeSuffixVar('CMS_scale_met_2017', 'MET')
makeSuffixVar('CMS_scale_met_2017', 'MET', mc_forMET, as_lnN=True)

#nuisances['CMS_scale_met_2017']  = {
#    'name'  : 'CMS_scale_met_2017',
#    'type'  : 'lnN',
#    'samples'  : {
#        #'Wjets'  : '1.04',
#        #'WZewk'  : '1.0416406210649',
#        'WW'     : '1.03',
#        'DY'     : '1.1',
#        'Higgs'  : '1.06',
#        #'ZZ'     : '1.115497219168',
#        'VVV'    : '1.01',
#        'ggWW'   : '1.03',
#        'Vg'     : '1.04',
#        'VBF-V'  : '1.03',
#        #'top'    : '1.02',
#        #'WWewk'  : '1.0734941788458',
#        #'WZqcd'  : '1.0679605960038',
#        'VgS_H'  : '1.01',
#        'VgS_L'  : '1.12',
#    },
#}
#for samp in mc:
#    if 'darkHiggs' in samp:
#        nuisances['CMS_scale_met_2017']['samples'][samp] = '1.006'

##### Pileup

pu_syst = '(puWeightUp/puWeight)', '(puWeightDown/puWeight)'

nuisances['PU'] = {
    'name': 'CMS_PU_2017',
    'kind': 'weight',
    'type': 'shape',
    #'samples': dict((skey, pu_syst) for skey in mc),
    'samples': {
        'Wjets': ['0.991238056241*(puWeightUp/puWeight)',  '1.00820697786*(puWeightDown/puWeight)'],
        'WjetsNLO': ['0.991238056241*(puWeightUp/puWeight)',  '1.00820697786*(puWeightDown/puWeight)'],
        'VBF-V': ['0.999925404198*(puWeightUp/puWeight)',  '0.999773846632*(puWeightDown/puWeight)'],
        'top'  : ['0.998758134732*(puWeightUp/puWeight)',  '1.00128739105*(puWeightDown/puWeight)'],
        'DY'   : ['0.985597078151*(puWeightUp/puWeight)',  '1.01399166222*(puWeightDown/puWeight)'],
        'DYlow': ['0.985970121345*(puWeightUp/puWeight)',  '1.012782071*(puWeightDown/puWeight)'],
        'WW'   : ['0.997425366796*(puWeightUp/puWeight)',  '1.00272616036*(puWeightDown/puWeight)'],
        'ggWW' : ['0.998088444607*(puWeightUp/puWeight)',  '1.00190916004*(puWeightDown/puWeight)'],
        'WWewk': ['0.998950556861*(puWeightUp/puWeight)',  '1.00139453776*(puWeightDown/puWeight)'],
        'WZqcd': ['0.997864609854*(puWeightUp/puWeight)',  '1.00212985184*(puWeightDown/puWeight)'],
        'WZewk': ['0.998037316045*(puWeightUp/puWeight)',  '1.00194969263*(puWeightDown/puWeight)'],
        'ZZ'   : ['0.996120500023*(puWeightUp/puWeight)',  '1.00295839391*(puWeightDown/puWeight)'],
        'Vg'   : ['0.988229767437*(puWeightUp/puWeight)',  '1.01187308309*(puWeightDown/puWeight)'],
        'VVV'  : ['0.997103334721*(puWeightUp/puWeight)',  '1.00315476639*(puWeightDown/puWeight)'],
        'Higgs': ['0.994732224977*(puWeightUp/puWeight)',  '1.00505623833*(puWeightDown/puWeight)'],
    },
    'AsLnN': '1',
}
covered_samples = nuisances['PU']['samples'].keys()
for skey in mc:
    if skey not in covered_samples: nuisances['PU']['samples'][skey] = ['(puWeightUp/puWeight)', '(puWeightDown/puWeight)']

#nuisances['JetPUID_sf']  = {
#    'name'  : 'CMS_jetpuid_2017',
#    #'type': 'lnN',
#    #'samples': dict((skey, '0.96/1.001') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top']),
#    'kind'  : 'weight',
#    'type'  : 'shape',
#    'samples'  : dict((skey, ['PUJetIdSF_up/PUJetIdSF','PUJetIdSF_down/PUJetIdSF']) for skey in mc ),
#    'AsLnN': '1',
#}

#handle = open('../WUpDown/JetPUID_2017_cfg.py', 'r')
handle = open('/afs/cern.ch/work/s/svanputt/monoHiggs/PR/CMSSW_10_6_5/src/PlotsConfigurations/Configurations/monoHWW/SemiLep/Full2017_v7/WUpDown/JetPUID_2017_cfg.py', 'r')
exec(handle)
handle.close()

puid_dict = {}
for samp in mc:
    if samp in jetpuid_dict:
        puid_dict[samp] = [jetpuid_dict[samp]['Up'], jetpuid_dict[samp]['Down']]
        #puid_dict[samp] = ['1.', jetpuid_dict[samp]['Down']]
    elif 'darkHiggs' in samp:
        puid_dict[samp] = [jetpuid_dict['darkHiggs']['Up'], jetpuid_dict['darkHiggs']['Down']]
        #puid_dict[samp] = ['1.', jetpuid_dict['darkHiggs']['Down']]
    elif 'VgS' in samp:
        puid_dict[samp] = [jetpuid_dict['VgS_H']['Up'], jetpuid_dict['VgS_H']['Down']]
        #puid_dict[samp] = ['1.', jetpuid_dict['VgS_H']['Down']]

    # Normalize Wjets and top
    elif 'Wjets' in samp:
        # Variations normalized to the W+jets CR
        # We have a rateParam to fix the normalization of the W+jets
        puid_dict[samp] = ['('+jetpuid_dict['Wjets']['Up']+')*1.00215576533', '('+jetpuid_dict['Wjets']['Down']+')*1.01562046258']
    elif 'top' in samp:
        # Variations normalized to the top CR
        # We have a rateParam to fix the normalization of the top
        puid_dict[samp] = ['('+jetpuid_dict[samp]['Up']+')*0.999679663232', '('+jetpuid_dict[samp]['Down']+')*1.00676613138']

nuisances['JetPUID_fake_sf']  = {
    'name'  : 'CMS_jetpuid_2017',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': puid_dict,
    #'samples': dict((skey, [jetpuid_dict[skey]['Up'], jetpuid_dict[skey]['Down']]) for skey in jetpuid_dict)
    #'samples'  : {
    #    'top'  : jetPUID_top,
    #    'Wjets': jetPUID_Wjets,
    #} 
    #'AsLnN': '1',
    #'symmetrize': True
}


##### PS and UE
# Shapes broken in 2017, propagate shapes of 2018
PSw_samples = ['DY','WW','VVV','Vg','VgS','VBF-V','ggWW']
nuisances['PS_ISR']  = {
    'name'  : 'PS_ISR',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': dict((skey, ['1.', '1.']) for skey in PSw_samples)
}
nuisances['PS_FSR']  = {
    'name'  : 'PS_FSR',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': dict((skey, ['1.', '1.']) for skey in PSw_samples)
}

# Wjets
nuisances['PS_ISR_Wjets']  = {
    'name'  : 'PS_ISR_Wjets',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'Wjets' : ['1.', '1.'],
        'WjetsNLO' : ['1.', '1.'],
    }
}
# FSR for W+jets was very noisy, and approx flat => lnN, same size in all regions => removed (absorbed in rateParam)

# top
nuisances['PS_ISR_top']  = {
    'name'  : 'PS_ISR_top',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'top' : ['1.', '1.'],
    }
}
nuisances['PS_FSR_top']  = {
    'name'  : 'PS_FSR_top',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'top' : ['1.', '1.'],
    }
}

## Shapes broken in 2017, propagate shapes of 2018
#samples_PS = ['WjetsNLO', 'Wjets','top','DY','WW','VVV','Vg','VgS','VBF-V','ggWW']
#for samp in samples_PS:
#    nuisances['PS_ISR_'+samp]  = {
#        'name'  : 'PS_ISR_'+samp,
#        'kind'  : 'weight',
#        'type'  : 'shape',
#        'samples'  : {
#            #samp : ['PSWeight[2]', 'PSWeight[0]'],
#            samp : ['1.', '1.'],
#        }
#    }
#    
#    nuisances['PS_FSR_'+samp]  = {
#        'name'  : 'PS_FSR_'+samp,
#        'kind'  : 'weight',
#        'type'  : 'shape',
#        'samples'  : {
#            #samp : ['PSWeight[3]', 'PSWeight[1]'],
#            samp : ['1.', '1.'],
#        }
#    }
#
#    if samp in ['Vg']:
#        nuisances['PS_ISR_'+samp]['AsLnN'] = '1'
#        nuisances['PS_FSR_'+samp]['AsLnN'] = '1'
#    if samp in ['Wjets', 'WjetsNLO']:
#        nuisances['PS_FSR_'+samp]['AsLnN'] = '1'

#nuisances['PS_ISR']  = {
#    'name': 'PS_ISR',
#    'kind': 'weight',
#    'type': 'shape',
#    #'AsLnN': '1',
#    'samples': {
#        'Wjets': ['0.999361844313*(nCleanGenJet==0) + 1.00932187819*(nCleanGenJet==1) + 1.01036957312*(nCleanGenJet==2) + 0.980802963847*(nCleanGenJet>=3)', 
#                  '1.00133729681*(nCleanGenJet==0) + 0.989165456129*(nCleanGenJet==1) + 0.988125359076*(nCleanGenJet==2) + 1.02558360248*(nCleanGenJet>=3)'],
#        #'WZewk': ['*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)', 
#        #          '*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)'],
#        #'WZqcd': ['*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)', 
#        #          '*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)'],
#        'Vg'   : ['1.00227428567253*(nCleanGenJet==0) + 1.00572014989997*(nCleanGenJet==1) + 0.970824885256465*(nCleanGenJet==2) + 0.927346068071086*(nCleanGenJet>=3)', 
#                  '0.996488506572636*(nCleanGenJet==0) + 0.993582795375765*(nCleanGenJet==1) + 1.03643678934568*(nCleanGenJet==2) + 1.09735277266955*(nCleanGenJet>=3)'],
#        'VgS'  : ['1.0000536116408023*(nCleanGenJet==0) + 1.0100100693580492*(nCleanGenJet==1) + 0.959068359375*(nCleanGenJet==2) + 0.9117049260469496*(nCleanGenJet>=3)', 
#                  '0.9999367833485968*(nCleanGenJet==0) + 0.9873682892005163*(nCleanGenJet==1) + 1.0492717737268518*(nCleanGenJet==2) + 1.1176958835210322*(nCleanGenJet>=3)'],
#        'ggWW' : ['1.040233912070831*(nCleanGenJet==0) + 0.9611236379290876*(nCleanGenJet==1) + 0.9014289294088699*(nCleanGenJet==2) + 0.864310738090035*(nCleanGenJet>=3)', 
#                  '0.9510305474211223*(nCleanGenJet==0) + 1.0433432942960381*(nCleanGenJet==1) + 1.1271383507266095*(nCleanGenJet==2) + 1.1885756983901514*(nCleanGenJet>=3)'],
#        'WW'   : ['1.0005237869294796*(nCleanGenJet==0) + 1.0157425373134328*(nCleanGenJet==1) + 0.9644598124510606*(nCleanGenJet==2) + 0.9271488926223369*(nCleanGenJet>=3)', 
#                  '0.9993553300024391*(nCleanGenJet==0) + 0.9806102300995024*(nCleanGenJet==1) + 1.042603303739856*(nCleanGenJet==2) + 1.0950369125887705*(nCleanGenJet>=3)'],
#        'WWewk': ['1.000046239608802*(nCleanGenJet==0) + 1.0133437345113673*(nCleanGenJet==1) + 1.0144565712441801*(nCleanGenJet==2) + 0.9633224709814209*(nCleanGenJet>=3)', 
#                  '0.9972924400977996*(nCleanGenJet==0) + 0.9796025428294364*(nCleanGenJet==1) + 0.9785038122896787*(nCleanGenJet==2) + 1.043894692205412*(nCleanGenJet>=3)'], # From WpWmJJ_EWK_noTop
#        'top'  : ['1.0020618369910668*(nCleanGenJet==0) + 1.0063081530771556*(nCleanGenJet==1) + 1.0094298425968304*(nCleanGenJet==2) + 0.9854207999040726*(nCleanGenJet>=3)', 
#                  '0.9974340279269026*(nCleanGenJet==0) + 0.9920634820709106*(nCleanGenJet==1) + 0.988226385054923*(nCleanGenJet==2) + 1.017968568319235*(nCleanGenJet>=3)'],
#        'DY'   : ['0.9998177685645392*(nCleanGenJet==0) + 1.0080838149428026*(nCleanGenJet==1) + 1.0057948912950987*(nCleanGenJet==2) + 0.9721358221196619*(nCleanGenJet>=3)', 
#                  '1.0003244155266309*(nCleanGenJet==0) + 0.9897992135367016*(nCleanGenJet==1) + 0.9928782069009531*(nCleanGenJet==2) + 1.0348902921423981*(nCleanGenJet>=3)'],
#        'VVV'  : ['1.0270826786253018*(nCleanGenJet==0) + 1.0198703447307862*(nCleanGenJet==1) + 1.0109191915514344*(nCleanGenJet==2) + 0.9838184220287978*(nCleanGenJet>=3)', 
#                  '0.9661665482954546*(nCleanGenJet==0) + 0.9751744967838527*(nCleanGenJet==1) + 0.9859624782745712*(nCleanGenJet==2) + 1.0202995039288625*(nCleanGenJet>=3)'],
#        'VZ'   : ['1.0005649495173758*(nCleanGenJet==0) + 1.0126108996187784*(nCleanGenJet==1) + 0.9655814439488959*(nCleanGenJet==2) + 0.9231853653899337*(nCleanGenJet>=3)', 
#                  '0.9992093891493296*(nCleanGenJet==0) + 0.9843983411237471*(nCleanGenJet==1) + 1.0409271492937164*(nCleanGenJet==2) + 1.1001990664321764*(nCleanGenJet>=3)'], # From ZZTo2L2Nu_ext1
#        'Higgs': ['1.0007510488273352*(nCleanGenJet==0) + 1.0152476349471342*(nCleanGenJet==1) + 0.9645590929269297*(nCleanGenJet==2) + 0.9189171704206691*(nCleanGenJet>=3)', 
#                    '0.9989909143752528*(nCleanGenJet==0) + 0.9814978813068076*(nCleanGenJet==1) + 1.0416554335980368*(nCleanGenJet==2) + 1.1060543963750413*(nCleanGenJet>=3)'],
#    },
#}
#
#nuisances['PS_FSR']  = {
#    'name': 'PS_FSR',
#    'kind': 'weight',
#    'type': 'shape',
#    #'AsLnN': '1',
#    'samples': {
#        'Wjets': ['0.951581197919*(nCleanGenJet==0) + 0.997755474747*(nCleanGenJet==1) + 1.02983517401*(nCleanGenJet==2) + 1.01143623738*(nCleanGenJet>=3)', 
#                   '1.08454892632*(nCleanGenJet==0) + 1.0003688082*(nCleanGenJet==1) + 0.980741215888*(nCleanGenJet==2) + 0.97189538569*(nCleanGenJet>=3)'],
#        #'WZewk': ['*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)', 
#        #          '*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)'],
#        #'WZqcd': ['*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)', 
#        #          '*(nCleanGenJet==0) + *(nCleanGenJet==1) + *(nCleanGenJet==2) + *(nCleanGenJet>=3)'],
#        'Vg'     : ['0.999935529935028*(nCleanGenJet==0) + 0.997948255568351*(nCleanGenJet==1) + 1.00561645493085*(nCleanGenJet==2) + 1.0212896960035*(nCleanGenJet>=3)', 
#                    '1.00757702771109*(nCleanGenJet==0) + 1.00256681166083*(nCleanGenJet==1) + 0.93676371569867*(nCleanGenJet==2) + 0.956448336052435*(nCleanGenJet>=3)'],
#        'VgS'    : ['0.9976593177227735*(nCleanGenJet==0) + 1.0016125187585532*(nCleanGenJet==1) + 1.0049344618055556*(nCleanGenJet==2) + 1.0195631514301164*(nCleanGenJet>=3)', 
#                    '1.0026951855766457*(nCleanGenJet==0) + 1.0008132148661049*(nCleanGenJet==1) + 1.003949291087963*(nCleanGenJet==2) + 0.9708160910230832*(nCleanGenJet>=3)'],
#        'ggWW'   : ['0.9910563426395067*(nCleanGenJet==0) + 1.0069894351287263*(nCleanGenJet==1) + 1.016616376034912*(nCleanGenJet==2) + 1.015902717074592*(nCleanGenJet>=3)', 
#                    '1.0147395976461193*(nCleanGenJet==0) + 0.9860219489006646*(nCleanGenJet==1) + 0.9694680606617647*(nCleanGenJet==2) + 0.9489845115821678*(nCleanGenJet>=3)'],
#        'WW'     : ['0.995462478372054*(nCleanGenJet==0) + 1.0052129975124378*(nCleanGenJet==1) + 1.008836750560578*(nCleanGenJet==2) + 0.9984120564941189*(nCleanGenJet>=3)', 
#                    '1.008927720738437*(nCleanGenJet==0) + 0.995163868159204*(nCleanGenJet==1) + 0.9911024228315418*(nCleanGenJet==2) + 0.9763787172658678*(nCleanGenJet>=3)'],
#        'WWewk'  : ['0.9919173789731052*(nCleanGenJet==0) + 0.994273375713824*(nCleanGenJet==1) + 1.0009783801226202*(nCleanGenJet==2) + 1.0046943184676753*(nCleanGenJet>=3)', 
#                    '1.0081780733496333*(nCleanGenJet==0) + 1.0052172826204072*(nCleanGenJet==1) + 0.9984786981975752*(nCleanGenJet==2) + 0.9881601609947174*(nCleanGenJet>=3)'],
#        'top'    : ['0.9910899786333963*(nCleanGenJet==0) + 0.9990635702054794*(nCleanGenJet==1) + 1.002141744200183*(nCleanGenJet==2) + 1.0129742776372779*(nCleanGenJet>=3)', 
#                    '1.0068843378231833*(nCleanGenJet==0) + 0.998988498438759*(nCleanGenJet==1) + 0.9952696584115224*(nCleanGenJet==2) + 0.9790955840673237*(nCleanGenJet>=3)'],
#        'DY'     : ['0.9958763409773141*(nCleanGenJet==0) + 1.0041335498093422*(nCleanGenJet==1) + 1.0163363150953029*(nCleanGenJet==2) + 1.0296733670670226*(nCleanGenJet>=3)', 
#                    '1.0066775262249232*(nCleanGenJet==0) + 0.9945601465681602*(nCleanGenJet==1) + 0.9662459619335311*(nCleanGenJet==2) + 0.9479423453563661*(nCleanGenJet>=3)'],
#        'VVV'    : ['0.9809047855490748*(nCleanGenJet==0) + 0.9823641498350338*(nCleanGenJet==1) + 0.9976414629808243*(nCleanGenJet==2) + 1.0077953569413387*(nCleanGenJet>=3)', 
#                    '1.035388723727876*(nCleanGenJet==0) + 1.0347339790465233*(nCleanGenJet==1) + 1.0017058788771533*(nCleanGenJet==2) + 0.9829344116371653*(nCleanGenJet>=3)'],
#        'VZ'     : ['0.9969613063607379*(nCleanGenJet==0) + 1.0055822932729375*(nCleanGenJet==1) + 1.0115617769727228*(nCleanGenJet==2) + 1.003930218451846*(nCleanGenJet>=3)', 
#                    '1.0044825167062181*(nCleanGenJet==0) + 0.9934087678296067*(nCleanGenJet==1) + 0.9816200440919792*(nCleanGenJet==2) + 0.9733590995315957*(nCleanGenJet>=3)'],
#        'Higgs'  : ['0.9936588910230489*(nCleanGenJet==0) + 1.0087564198432573*(nCleanGenJet==1) + 1.014636529653396*(nCleanGenJet==2) + 1.00399261707105*(nCleanGenJet>=3)', 
#                    '1.0125063182369591*(nCleanGenJet==0) + 0.9846168672324244*(nCleanGenJet==1) + 0.9778204449152542*(nCleanGenJet==2) + 1.0014057292097962*(nCleanGenJet>=3)'],
#    },
#}

# An overall 1.5% UE uncertainty will cover all the UEup/UEdo variations
# And we don't observe any dependency of UE variations on njet
nuisances['UE']  = {
    'name'  : 'UE_CP5',
    'skipCMS' : 1,
    'type': 'lnN',
    'samples': dict((skey, '1.015') for skey in mc if skey not in ['WjetsNLO', 'Wjets', 'top']), 
}

####### Generic "cross section uncertainties"
lastcopy     = (1 << 13)
topGenPt     = '(Sum$((GenPart_pdgId == 6 && (GenPart_statusFlags & %d)) * GenPart_pt))' % lastcopy
antitopGenPt = '(Sum$((GenPart_pdgId == -6 && (GenPart_statusFlags & %d)) * GenPart_pt))' % lastcopy

nuisances['singleTopToTTbar'] = {
    'name': 'singleTopToTTbar',
    'skipCMS': 1,
    'kind': 'weight',
    'type': 'shape',
    'samples': {
        'top': [
            '('+topGenPt+' * '+antitopGenPt+' <= 0.) * 1.0816 + ('+topGenPt+' * '+antitopGenPt+' > 0.)',
            '('+topGenPt+' * '+antitopGenPt+' <= 0.) * 0.9184 + ('+topGenPt+' * '+antitopGenPt+' > 0.)'
        ],
    }
}

#nuisances['Wjets_merge'] = {
#    'name': 'Wjets_merge_2017',
#    'skipCMS': 1,
#    'kind': 'weight',
#    'type': 'shape',
#    'samples': {
#        'Wjets': [
#            '(LHE_Vpt > 120)*1.05 + (LHE_Vpt < 120)*1.',
#            '(LHE_Vpt > 120)*0.95 + (LHE_Vpt < 120)*1.'
#        ],
#    }
#}

## Top pT reweighting uncertainty
#FIXME: correct?
nuisances['TopPtRew'] = {
    'name': 'CMS_topPtRew',   # Theory uncertainty
    'kind': 'weight',
    'type': 'shape',
    'samples': {
        'top': ["1./Top_pTrw", "Top_pTrw"],
    },
    #'symmetrize': True
}

nuisances['VgStar'] = {
    'name': 'CMS_hww_VgStarScale',
    'type': 'lnN',
    'samples': {
        'VgS_L': '1.25'
    }
}

nuisances['VZ'] = {
    'name': 'CMS_hww_VZScale',
    'type': 'lnN',
    'samples': {
        'VgS_H': '1.16'
    }
}

###### pdf uncertainties
nuisances['pdf_weight'] = { 
    'name'  : 'pdf_weight_2017_18',
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'samples' :  { s: [' Alt$(LHEPdfWeight['+str(i)+'], 1.)' for i in range(0,103)] for s in mc if s not in ['top','Wjets', 'WjetsNLO']},
    'AsLnN':  '1' 
}

#nuisances['pdf']  = {
#    'name'  : 'pdf',
#    'type'  : 'lnN',
#    'samples'  : {
#        'ggWW'    : '1.05',
#        'WW'      : '1.04',
#        'Vg'      : '1.04',
#        'VZ'      : '1.04',
#        'VgS'     : '1.04',
#        'Higgs'   : '1.04',
#        'DY'      : '1.002', 
#    },
#}

##### JER influence on MET
nuisances['jer_to_met']  = {
    'name'  : 'jer_to_met_2017',
    'type'  : 'lnN',
    'samples'  : {
        #'Wjets'  : '1.111188965734',
        'WZewk'  : '1.0416406210649',
        'WW'     : '1.0867795440189',
        'DY'     : '1.0857734597136',
        'Higgs'  : '1.0888269630155',
        'ZZ'     : '1.115497219168',
        'VVV'    : '1.0552524181545',
        'ggWW'   : '1.0995578177253',
        'Vg'     : '1.107239154749',
        'VBF-V'  : '1.0632429856598',
        #'top'    : '1.0441186588377',
        'WWewk'  : '1.0734941788458',
        'WZqcd'  : '1.0679605960038',
        'VgS_H'  : '1.0796135798483',
        'VgS_L'  : '1.107341335709',
    },
}


## pdf for separate Higgs
##FIXME: correct?
#valuesggh  = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','ggH', '125.09','pdf','sm')
#valuesggzh = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','ggZH','125.09','pdf','sm')
#valuesbbh  = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','bbH', '125.09','pdf','sm')
#
#nuisances['pdf_Higgs_gg'] = {
#    'name': 'pdf_Higgs_gg',
#    'samples': {
#        'ggH_hww' : valuesggh,
#        'ggH_htt' : valuesggh,
#        'ggZH_hww': valuesggzh,
#        'bbH_hww' : valuesbbh
#    },
#    'type': 'lnN',
#}
#
#valuesqqh = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','vbfH','125.09','pdf','sm')
#valueswh = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','WH','125.09','pdf','sm')
#valueszh = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','ZH','125.09','pdf','sm')
#
#nuisances['pdf_Higgs_qqbar'] = {
#    'name': 'pdf_Higgs_qqbar',
#    'type': 'lnN',
#    'samples': {
#        'qqH_hww': valuesqqh,
#        'qqH_htt': valuesqqh,
#        'WH_hww': valueswh,
#        'WH_htt': valueswh,
#        'ZH_hww': valueszh,
#        'ZH_htt': valueszh
#    },
#}
#
#values = HiggsXS.GetHiggsProdXSNP('YR4','13TeV','ttH','125.09','pdf','sm')
#
#nuisances['pdf_Higgs_ttH'] = {
#    'name': 'pdf_Higgs_ttH',
#    'samples': {
#        'ttH_hww': values
#    },
#    'type': 'lnN',
#}
#
#nuisances['pdf_Higgs_gg_ACCEPT'] = {
#    'name': 'pdf_Higgs_gg_ACCEPT',
#    'samples': {
#        'ggH_hww': '1.006',
#        'ggH_htt': '1.006',
#        'ggZH_hww': '1.006',
#        'bbH_hww': '1.006'
#    },
#    'type': 'lnN',
#}
#
#nuisances['pdf_gg_ACCEPT'] = {
#    'name': 'pdf_gg_ACCEPT',
#    'samples': {
#        'ggWW': '1.006',
#    },
#    'type': 'lnN',
#}
#
#nuisances['pdf_Higgs_qqbar_ACCEPT'] = {
#    'name': 'pdf_Higgs_qqbar_ACCEPT',
#    'type': 'lnN',
#    'samples': {
#        'qqH_hww': '1.002',
#        'qqH_htt': '1.002',
#        'WH_hww': '1.003',
#        'WH_htt': '1.003',
#        'ZH_hww': '1.002',
#        'ZH_htt': '1.002',
#    },
#}
#
#nuisances['pdf_qqbar_ACCEPT'] = {
#    'name': 'pdf_qqbar_ACCEPT',
#    'type': 'lnN',
#    'samples': {
#        'VZ': '1.001',
#    },
#}

###### Renormalization & factorization scales
## LHE scale variation weights (w_var / w_nominal)
#   [0] is ren_sf=0.5 fac_sf=0.5    DD
#   [1] is ren_sf=0.5 fac_sf=1      DN
#   [2] is ren_sf=0.5 fac_sf=2      DU  * -> Extreme case, we don't use
#   [3] is ren_sf=1   fac_sf=0.5    ND
#   [4] is ren_sf=1   fac_sf=1      NN -> Not always there
#   [5] is ren_sf=1   fac_sf=2      NU
#   [6] is ren_sf=2   fac_sf=0.5    UD  * -> Extreme case, we don't use
#   [7] is ren_sf=2   fac_sf=1      UN
#   [8] is ren_sf=2   fac_sf=2      UU

variations_ren  = ['LHEScaleWeight[Length$(LHEScaleWeight)-2]', 'LHEScaleWeight[1]']
variations_fac  = ['LHEScaleWeight[Length$(LHEScaleWeight)-4]', 'LHEScaleWeight[3]']
variations_both = ['LHEScaleWeight[Length$(LHEScaleWeight)-1]', 'LHEScaleWeight[0]']

variations = [
    'LHEScaleWeight[0]', 
    'LHEScaleWeight[1]', 
    'LHEScaleWeight[3]', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-4]', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-2]', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-1]'
]

nuisances['QCDscale_V'] = {
    'name': 'QCDscale_V',
    'kind': 'weight_envelope',
    'type': 'shape',
    'samples': {'DY': variations},
}

nuisances['QCDscale_VV'] = {
    'name': 'QCDscale_VV',
    'kind': 'weight_envelope',
    'type': 'shape',
    'samples': {
        'Vg': variations,
        # 'VZ': variations, #FIXME not all VZ have LHEScaleWeight
        'VgS': variations
    }
}

# Variations normalized to the top CR
# We have a rateParam to fix the normalization of the top
variations_top = [
    'LHEScaleWeight[0]*0.886682399441', 
    'LHEScaleWeight[1]*0.901663638362', 
    'LHEScaleWeight[3]*0.978408052981', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-4]*1.01695785314', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-2]*1.11080446146', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-1]*1.13460524274'
]
nuisances['QCDscale_top']  = {
    'name'  : 'QCDscale_top', 
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'samples'  : {
        #'top' : variations,
        'top' : variations_top,
    }
}

# Variations normalized to the W+jets CR
# We have a rateParam to fix the normalization of the W+jets
variations_wjets = [
    'LHEScaleWeight[0]*0.997971794172', 
    'LHEScaleWeight[1]*0.95528281743', 
    'LHEScaleWeight[3]*1.04654834119', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-4]*0.971535563141', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-2]*1.03909780122', 
    'LHEScaleWeight[Length$(LHEScaleWeight)-1]*1.00826470248'
]
variations_ren  = ['LHEScaleWeight[Length$(LHEScaleWeight)-2]*1.03909780122' , 'LHEScaleWeight[1]*0.95528281743']
variations_fac  = ['LHEScaleWeight[Length$(LHEScaleWeight)-4]*0.971535563141', 'LHEScaleWeight[3]*1.04654834119']
variations_both = ['LHEScaleWeight[Length$(LHEScaleWeight)-1]*1.00826470248' , 'LHEScaleWeight[0]*0.997971794172']
nuisances['QCDscale_Wjets']  = {
    'name'  : 'QCDscale_Wjets', 
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'samples'  : {
        #'Wjets' : variations,
        'Wjets' : variations_wjets,
        'WjetsNLO' : variations_wjets,
    }
}
#nuisances['QCDscale_ren_Wjets']  = {
#    'name'  : 'QCDscale_ren_Wjets', 
#    'kind'  : 'weight',
#    'type'  : 'shape',
#    'samples'  : {
#        'Wjets' : variations_ren,
#    }
#}
#nuisances['QCDscale_fac_Wjets']  = {
#    'name'  : 'QCDscale_fac_Wjets', 
#    'kind'  : 'weight',
#    'type'  : 'shape',
#    'samples'  : {
#        'Wjets' : variations_fac,
#    }
#}
#nuisances['QCDscale_both_Wjets']  = {
#    'name'  : 'QCDscale_both_Wjets', 
#    'kind'  : 'weight',
#    'type'  : 'shape',
#    'samples'  : {
#        'Wjets' : variations_both,
#    }
#}

nuisances['QCDscale_WWewk']  = {
    'name'  : 'QCDscale_WWewk',
    'type'  : 'lnN',
    'samples'  : {
        'WWewk' : '1.11',
    },
}

nuisances['QCDscale_ggVV'] = {
    'name': 'QCDscale_ggVV',
    'type': 'lnN',
    'samples': {
        'ggWW': '1.15',
    },
}

#### W+jets SF correction uncertainties
nuisances['Wjets_EWK_NLO_sf'] = {
    'name': 'Wjets_EWK_NLO_sf',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': {
        'Wjets': ['1./ewknloW', 'ewknloW'],
        'WjetsMjj': ['1./ewknloW', 'ewknloW'],
        'WjetsNLO': ['1./EWKnloW[0]', 'EWKnloW[0]'],
    }
}
nuisances['Wjets_QCD_NLO_sf_stat'] = {
    'name': 'Wjets_QCD_NLO_sf_stat_2017',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': {
        'Wjets': ['kfact[1]/kfact[0]', 'kfact[2]/kfact[0]'],
        'WjetsMjj': ['kfact[1]/kfact[0]', 'kfact[2]/kfact[0]'],
        'WjetsSP': ['kfactSP[1]/kfactSP[0]', 'kfactSP[2]/kfactSP[0]'],
    }
}
nuisances['Wjets_QCD_NLO_sf_jetsel'] = {
    'name': 'Wjets_QCD_NLO_sf_jetsel',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': {
        ##'Wjets': ['2*kfact[3]/kfact[0]', '1.'],
        #'Wjets': ['kfact[3]/kfact[0]', '(2.-kfact[3]/kfact[0])'],
        #'WjetsSP': ['kfactSP[3]/kfactSP[0]', '(2.-kfactSP[3]/kfactSP[0])'],

        # Variations normalized to the W+jets CR
        # We have a rateParam to fix the normalization of the W+jets
        'Wjets': ['(kfact[3]/kfact[0])*0.988003163255', '(2.-kfact[3]/kfact[0])*1.01229176125'],
        'WjetsMjj': ['(kfact[3]/kfact[0])*0.988003163255', '(2.-kfact[3]/kfact[0])*1.01229176125'],
    },
    #'symmetrize': True
}
nuisances['Wjets_QCD_NLO_sf_mjj'] = {
    #'name': 'Wjets_QCD_NLO_sf_mjj',
    'name': 'Wjets_QCD_NLO_sf_mjj_2017_18',
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples': {

        # Variations normalized to the W+jets CR
        # We have a rateParam to fix the normalization of the W+jets
        'Wjets': ['kfactMjj[0]', '1./kfactMjj[0]'],
        'WjetsMjj': ['kfactMjj[1]/kfactMjj[0]', 'kfactMjj[2]/kfactMjj[0]'],
    }
}


#nuisances['Wjets_corrSF'] = {
#    'name': 'Wjets_corrSF',
#    'kind'  : 'weight',
#    'type'  : 'shape',
#    'samples': {
#        #'Wjets': ['1./ewknloW', 'ewknloW'],
#        'Wjets_PuppiRW': ['1./Wjets_puppirw[0]', 'Wjets_puppirw[0]'],
#        'Wjets_HTsf': ['1./VptSF[0]', 'VptSF[0]'],
#    }
#}

# NLL resummation variations
nuisances['WWresum']  = {
  'name'  : 'CMS_hww_WWresum',
  'skipCMS' : 1,
  'kind'  : 'weight',
  'type'  : 'shape',
  'samples'  : {
     'WW'   : ['nllW_Rup/nllW', 'nllW_Rdown/nllW'],
   },
}

nuisances['WWqscale']  = {
   'name'  : 'CMS_hww_WWqscale',
   'skipCMS' : 1,
   'kind'  : 'weight',
   'type'  : 'shape',
   'samples'  : {
      'WW'   : ['nllW_Qup/nllW', 'nllW_Qdown/nllW'],
    },
}




## Use the following if you want to apply the automatic combine MC stat nuisances.
nuisances['stat'] = {
    'type': 'auto',
    'maxPoiss': '10',
    'includeSignal': '0',
    #  nuisance ['maxPoiss'] =  Number of threshold events for Poisson modelling
    #  nuisance ['includeSignal'] =  Include MC stat nuisances on signal processes (1=True, 0=False)
    'samples': {}
}

###rate parameters
nuisances['Topnorm']  = {
    'name'  : 'Topnorm_2017',
    'samples'  : {
        #'ttop' : '1.00',
        #'stop' : '1.00',
        'top' : '1.00',
    },
    'type'  : 'rateParam',
    'cuts'  : [
        'InCh_SR',
        'InCh_TCR',
        'InCh_SB',
    ]
}
    
nuisances['Wjetsnorm']  = {
    'name'  : 'Wjetsnorm_2017',
    'samples'  : {
        'Wjets' : '1.00',
        #'WjetsMjj' : '1.00',
    },
    'type'  : 'rateParam',
    'cuts'  : [
        'InCh_SR',
        'InCh_TCR',
        'InCh_SB',
    ]
}

for n in nuisances.values():
    n['skipCMS'] = 1

#print(' '.join(nuis['name'] for nname, nuis in nuisances.iteritems() if nname not in ('lumi', 'stat')))
