# plot configuration

# groupPlot = {}
# 
# Groups of samples to improve the plots.
# If not defined, normal plots is used
#
'''
groupPlot['top']  = {
                  'nameHR' : 'tW and t#bar{t}',
                  'isSignal' : 0,
                  'color': 400,   # kYellow
                  'samples'  : ['top']
              }

groupPlot['WW']  = {  
                  'nameHR' : 'WW',
                  'isSignal' : 0,
                  'color': 851, # kAzure -9 
                  'samples'  : ['WW']
              }


groupPlot['Fake']  = {
                  'nameHR' : 'nonprompt',
                  'isSignal' : 0,
                  'color': 921,    # kGray + 1
                  'samples'  : ['Fake_em' ,'Fake_me']
}
'''
'''
groupPlot['VV']  = {  
                  'nameHR' : 'Multiboson',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['Vg', 'VZ', 'VgS_L', 'VgS_H', 'VVV']
              }

groupPlot['Vg']  = {
                  'nameHR' : 'Vg',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['Vg']
              }

groupPlot['VZ']  = {
                  'nameHR' : 'VZ',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['VZ']
              }

groupPlot['VgS_L']  = {
                  'nameHR' : 'VgS_L',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['VgS_L']
              }

groupPlot['VgS_H']  = {
                  'nameHR' : 'VgS_H',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['VgS_H']
              }


groupPlot['VVV']  = {
                  'nameHR' : 'VVV',
                  'isSignal' : 0,
                  'color': 617,   # kViolet +1
                  'samples'  : ['VVV']
              }

groupPlot['DY']  = {  
                  'nameHR' : "DY",
                  'isSignal' : 0,
                  'color': 418,    # kGreen+2
                  'samples'  : ['DY'],
                  'scale'    : 1,               
              }


groupPlot['VZ']  = {  
                  'nameHR' : "VZ",
                  'isSignal' : 0,
                  'color'    : 617,   # kViolet + 1  
                  'samples'  : ['VZ']
              }

groupPlot['Vg']  = {  
                  'nameHR' : "V#gamma",
                  'isSignal' : 0,
                  'color'    : 810,   # kOrange + 10
                  'samples'  : ['Vg']
              }

groupPlot['VgS']  = {
                  'nameHR' : "V#gamma*",
                  'isSignal' : 0,
                  'color'    : 409,   # kGreen - 9
                  'samples'  : ['VgS_L','VgS_H']
              }

groupPlot['VVV']  = {  
                  'nameHR' : 'VVV',
                  'isSignal' : 0,
                  'color': 857, # kAzure -3  
                  'samples'  : ['VVV']
              }


'''

groupPlot['Higgs']  = {
                  'nameHR' : 'Gluon Fusion',
                  'isSignal' : 1,
                  'color': 409, # kRed 
                  #'samples'  : ['H_htt', 'H_hww', 'ZH_hww', 'ggZH_hww', 'WH_hww', 'ggH_hww','bbH_hww','ttH_hww','ZH_htt', 'ggZH_htt', 'WH_htt', 'ggH_htt','bbH_htt','ttH_htt' ]
                  'samples' : ['ZH_hww', 'ggZH_hww', 'WH_hww', 'ggH_hww','ttH_hww','ZH_htt', 'WH_htt', 'ggH_htt','qqH_htt' ]
                  #'samples'  : ['H_htt', 'H_hww', 'ZH_hww', 'ggZH_hww', 'WH_hww', 'ggH_hww','bbH_hww','ttH_hww', 'ggH_htt' ]
              }
groupPlot['VBF']  = {
                  'nameHR' : 'VBF',
                  'isSignal' : 2,
                  'color': 632,
                  'samples'  : ['qqH_hww']
              }
'''

# 3 puts signal on top of stack
# 2 hist only, included in ratio?
# 1 on top of stack and hist also, included in ratio?

groupPlot['HSM']  = {
                  'nameHR' : 'SM h',
                  'isSignal' : 0,
                  'color': 632,
                  'samples'  : ['H0PM','ZH_H0PM','WH_H0PM','VBF_H0PM']
              }

groupPlot['HBSM']  = {
                  'nameHR' : '0^{-}',
                  'isSignal' : 2,
                  'color': 1,
                  'samples'  : ['H0M','ZH_H0M','WH_H0M','VBF_H0M']
              }

groupPlot['HBSM2']  = {
                  'nameHR' : '0^{+}',
                  'isSignal' : 2,
                  'color': 2,
                  'samples'  : ['H0PH','ZH_H0PH','WH_H0PH','VBF_H0PH']
              }
'''
groupPlot['HBSM3']  = {
                  'nameHR' : 'L1',
                  'isSignal' : 2,
                  'color': 3,
                  'samples'  : ['H0L1','ZH_H0L1','WH_H0L1','VBF_H0L1']
              }
'''
#plot = {}

# keys here must match keys in samples.py    
#
                    
plot['DY']  = {  
                  'color': 418,    # kGreen+2
                  'isSignal' : 0,
                  'isData'   : 0, 
                  'scale'    : 0.77,
                  #'cuts'  : {
                       #'hww2l2v_13TeV_of0j'      : 0.95 ,
                       #'hww2l2v_13TeV_top_of0j'  : 0.95 , 
                       #'hww2l2v_13TeV_dytt_of0j' : 0.95 ,
                       #'hww2l2v_13TeV_em_0j'     : 0.95 , 
                       #'hww2l2v_13TeV_me_0j'     : 0.95 , 
                       ##
                       #'hww2l2v_13TeV_of1j'      : 1.08 ,
                       #'hww2l2v_13TeV_top_of1j'  : 1.08 , 
                       #'hww2l2v_13TeV_dytt_of1j' : 1.08 ,
                       #'hww2l2v_13TeV_em_1j'     : 1.08 , 
                       #'hww2l2v_13TeV_me_1j'     : 1.08 , 
                        #},

              }
'''
'''
if useEmbeddedDY:
  plot['Dyemb']  = {  
                  'color': 418,    # kGreen+2
                  'isSignal' : 0,
                  'isData'   : 0, 
                  'scale'    : 0.77,
              }
'''
'''
plot['Fake_me']  = {
                  'color': 921,    # kGray + 1
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
              }


plot['Fake_em']  = {
                  'color': 921,    # kGray + 1
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
              }

'''            
plot['top'] = {   
                  'nameHR' : 'tW and t#bar{t}',
                  'color': 400,   # kYellow
                  'isSignal' : 0,
                  'isData'   : 0, 
                  'scale'    : 1.0,
                  #'cuts'  : {
                       #'hww2l2v_13TeV_of0j'      : 0.94 ,
                       #'hww2l2v_13TeV_top_of0j'  : 0.94 , 
                       #'hww2l2v_13TeV_dytt_of0j' : 0.94 ,
                       #'hww2l2v_13TeV_em_0j'     : 0.94 , 
                       #'hww2l2v_13TeV_me_0j'     : 0.94 , 
                       ##
                       #'hww2l2v_13TeV_of1j'      : 0.86 ,
                       #'hww2l2v_13TeV_top_of1j'  : 0.86 , 
                       #'hww2l2v_13TeV_dytt_of1j' : 0.86 ,
                       #'hww2l2v_13TeV_em_1j'     : 0.86 , 
                       #'hww2l2v_13TeV_me_1j'     : 0.86 , 
                        #},
                  }
'''

plot['WW']  = {
                  'color': 851, # kAzure -9 
                  'isSignal' : 0,
                  'isData'   : 0,    
                  'scale'    : 1.0   # ele/mu trigger efficiency   datadriven
                  }
'''
plot['ggWW']  = {
                  'color': 850, # kAzure -10
                  'isSignal' : 0,
                  'isData'   : 0,    
                  'scale'    : 1.0
                  }

plot['Vg']  = {
                  'color': 859, # kAzure -10
                  'isSignal' : 0,
                  'isData'   : 0,    
                  'scale'    : 1.0
                  }
plot['VgS_L'] = {
                  'color'    : 617,   # kViolet + 1  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
                  }
plot['VgS_H'] = {
                  'color'    : 618,   # kViolet + 1  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
                  }


plot['VZ']  = { 
                  'color': 858, # kAzure -2  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
                  }

plot['VVV']  = { 
                  'color': 857, # kAzure -3  
                  'isSignal' : 0,
                  'isData'   : 0,
                  'scale'    : 1.0
                  }

# Htautau


plot['ZH_htt'] = {
                  'nameHR' : 'ZHtt',
                  'color': 632+3, # kRed+3 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }


#plot['bbH_htt'] = {
#                  'nameHR' : 'bbHtt',
#                  'color': 632-1, # kRed-1 
#                  'isSignal' : 1,
#                  'isData'   : 0,
#                  'scale'    : 1    #
#                  }
#
#plot['ttH_htt'] = {
#                  'nameHR' : 'bbHtt',
#                  'color': 632-2, # kRed-1 
#                  'isSignal' : 1,
#                  'isData'   : 0,
#                  'scale'    : 1    #
#                  }
#
#
#plot['ggZH_htt'] = {
#                  'nameHR' : 'ggZHtt',
#                  'color': 632+4, # kRed+4
#                  'isSignal' : 1,
#                  'isData'   : 0,    
#                  'scale'    : 1    #
#                  }
#


plot['WH_htt'] = {
                  'nameHR' : 'WHtt',
                  'color': 632+2, # kRed+2 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }


#plot['qqH_htt'] = {
#                  'nameHR' : 'qqHtt',
#                  'color': 632+1, # kRed+1 
#                  'isSignal' : 1,
#                  'isData'   : 0,    
#                  'scale'    : 1    #
#                  }


plot['ggH_htt'] = {
                  'nameHR' : 'ggHtt',
                  'color': 632, # kRed 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }

# HWW 

plot['H0PM']  =   {
                      'nameHR' : 'h',
                      'color' : 620+1,
                      'isSignal' : 0, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['ZH_H0PM']  =   {
                      'nameHR' : 'ZH h',
                      'color' : 620+2,
                      'isSignal' : 0, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['WH_H0PM']  =   {
                      'nameHR' : 'WH h',
                      'color' : 620+3,
                      'isSignal' : 0, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['VBF_H0PM']  =   {
                      'nameHR' : 'VBF h',
                      'color' : 620+4,
                      'isSignal' : 0, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }


plot['H0M']  =   {
                      'nameHR' : '0^{-}',
                      'color' : 632+1,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['ZH_H0M']  =   {
                      'nameHR' : 'ZH 0^{-}',
                      'color' : 632+2,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['WH_H0M']  =   {
                      'nameHR' : 'WH 0^{-}',
                      'color' : 632+3,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['VBF_H0M']  =   {
                      'nameHR' : 'VBF 0^{-}',
                      'color' : 632+4,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }


plot['H0PH']  =   {
                      'nameHR' : '0^{+}',
                      'color' : 632+1,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['ZH_H0PH']  =   {
                      'nameHR' : 'ZH 0^{+}',
                      'color' : 632+2,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['WH_H0PH']  =   {
                      'nameHR' : 'WH 0^{+}',
                      'color' : 632+3,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['VBF_H0PH']  =   {
                      'nameHR' : 'VBF 0^{+}',
                      'color' : 632+4,
                      'isSignal' : 2, 
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['H0L1']  =   {
                      'nameHR' : 'L1',
                      'color' : 632+1,
                      'isSignal' : 2,
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['ZH_H0L1']  =   {
                      'nameHR' : 'ZH L1',
                      'color' : 632+2,
                      'isSignal' : 2,
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['WH_H0L1']  =   {
                      'nameHR' : 'WH L1',
                      'color' : 632+3,
                      'isSignal' : 2,
                      'isData'   : 0,
                      'scale'    : 1,
                     }

plot['VBF_H0L1']  =   {
                      'nameHR' : 'VBF L1',
                      'color' : 632+4,
                      'isSignal' : 2,
                      'isData'   : 0,
                      'scale'    : 1,
                     }


#plot['H_hww'] = {
#                  'nameHR' : 'Hww',
#                  'color': 632, # kRed 
#                  'isSignal' : 1,
#                  'isData'   : 0,    
#                  'scale'    : 1    #
#                  }


plot['ZH_hww'] = {
                  'nameHR' : 'ZH',
                  'color': 632+3, # kRed+3 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }

plot['ggZH_hww'] = {
                  'nameHR' : 'ggZH',
                  'color': 632+4, # kRed+4
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }

plot['WH_hww'] = {
                  'nameHR' : 'WH',
                  'color': 632+2, # kRed+2 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }


plot['qqH_hww'] = {
                  'nameHR' : 'qqH',
                  'color': 632+1, # kRed+1 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }


plot['ggH_hww'] = {
                  'nameHR' : 'ggH',
                  'color': 632, # kRed 
                  'isSignal' : 1,
                  'isData'   : 0,    
                  'scale'    : 1    #
                  }

#plot['bbH_hww'] = {
#                  'nameHR' : 'bbH',
#                  'color': 632+5, # kRed+5 
#                  'isSignal' : 1,
#                  'isData'   : 0,
#                  'scale'    : 1    #
#                  }

plot['ttH_hww'] = {
                  'nameHR' : 'ttH',
                  'color': 632+6, # kRed+6
                  'isSignal' : 1,
                  'isData'   : 0,
                  'scale'    : 1    #
                  }



# data

plot['DATA']  = { 
                  'nameHR' : 'Data',
                  'color': 1 ,  
                  'isSignal' : 0,
                  'isData'   : 1 ,
                  'isBlind'  : 0
              }


# additional options

legend['lumi'] = 'L =  41.86/fb'

legend['sqrt'] = '#sqrt{s} = 13 TeV'




