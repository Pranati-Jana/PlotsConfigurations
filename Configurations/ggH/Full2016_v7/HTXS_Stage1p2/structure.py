# structure configuration for datacard

#structure = {}

# keys here must match keys in samples.py    
#                    
structure['DY']  = {  
                  'isSignal' : 0,
                  'isData'   : 0
              }

structure['Dyemb']  = {
                  'isSignal' : 0,
                  'isData'   : 0
              }

structure['Dyveto']  = {
                  'isSignal' : 0,
                  'isData'   : 0,
                  'removeFromCuts' : [ k for k in cuts ],
              }

structure['Wjets']  = {  
                  'isSignal' : 0,
                  'isData'   : 0 
              }

structure['Fake']  = {  
                  'isSignal' : 0,
                  'isData'   : 0 
              }

structure['Fake_em']  = {  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'removeFromCuts' : [ k for k in cuts if 'me' in k],
              }

structure['Fake_me']  = {  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'removeFromCuts' : [ k for k in cuts if 'em' in k],
              }

structure['ttbar'] = {   
                  'isSignal' : 0,
                  'isData'   : 0 
                  }


structure['singletop'] = {   
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['top'] = {   
                  'isSignal' : 0,
                  'isData'   : 0 
                  }


structure['WW']  = {
                  'isSignal' : 0,
                  'isData'   : 0    
                  }

structure['WWewk']  = {
                  'isSignal' : 0,
                  'isData'   : 0
                  }

structure['ggWW']  = {
                  'isSignal' : 0,
                  'isData'   : 0    
                  }

structure['ggWW_Int']  = {
                  'isSignal' : 0,
                  'isData'   : 0    
                  }

structure['Wg']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['Vg']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['VgS'] = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['VgS_L'] = {
                  'isSignal' : 0,
                  'isData'   : 0
                  }

structure['VgS_H'] = {
                  'isSignal' : 0,
                  'isData'   : 0
                  }

structure['Zg']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['VZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['WZ']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }


structure['VVV']  = { 
                  'isSignal' : 0,
                  'isData'   : 0 
                  }

structure['ZZ']  = {
                  'isSignal' : 0,
                  'isData'   : 0    
                  }


structure['ggH'] = {
                  'isSignal' : 1,
                  'isData'   : 0    
                  }

for signal in signals:
    structure[signal] = {
        'isSignal' : 1,
        'isData'   : 0
    }

# data

print(signals)

import pickle
with open('vbfDipoleScaleSTXS.pkl', 'rb') as handle:
    vbfDipoleScale = pickle.load(handle)

with open('STXS_fraction_equalization.pkl') as handle:
    STXS_frac_eq = pickle.load(handle)

for signal in signals:
    if 'qqH_hww' in signal:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
            'scaleSampleForDatacard' : {cut : vbfDipoleScale[signal][cut] * STXS_frac_eq[signal] * 1.03621 for cut in cuts if cut in vbfDipoleScale[signal].keys()},
        }
    elif 'ggH_hww' in signal:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
            'scaleSampleForDatacard' : {cut : STXS_frac_eq[signal] * 1.03364 for cut in cuts},
        }
    elif 'WH_hww' in signal:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
            'scaleSampleForDatacard' : {cut : STXS_frac_eq[signal] * 1.01724 for cut in cuts},
        }
    elif 'ZH_hww' in signal:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
            'scaleSampleForDatacard' : {cut : STXS_frac_eq[signal] * 1.01994 for cut in cuts},
        }
    elif 'ggZH_hww' in signal:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
            'scaleSampleForDatacard' : {cut : STXS_frac_eq[signal] * 1.02494 for cut in cuts},
        }
    else:
        structure[signal] = {
            'isSignal' : 1,
            'isData'   : 0,
        }

structure['DATA']  = { 
                  'isSignal' : 0,
                  'isData'   : 1
              }


for nuis in nuisances.itervalues():
  if 'cutspost' in nuis:
    nuis['cuts'] = nuis['cutspost'](nuis, cuts)

