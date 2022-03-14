# variables

#variables = {}
    
#'fold' : # 0 = not fold (default), 1 = fold underflowbin, 2 = fold overflow bin, 3 = fold underflow and overflow

variables['events']  = {   'name': '1',
                        'range' : (1,0,2),
                        'xaxis' : 'events',
                         'fold' : 3
                        }
'''
variables['BDTG_WH'] = { 'name': 'hww_WH3l_OSSF_mvaBDTG(Entry$,0)',
                        'range' : (10,-1,1),
                        'xaxis' : 'MVA_WH',
                        'fold' : 3,
                        'linesToAdd' : ['.L %s/src/PlotsConfigurations/Configurations/EFT/WH_AC/Full2017/hww_WH3l_OSSF_mvaBDTG.C+' % os.getenv('CMSSW_BASE')] 
                   }
'''
'''
variables['BDTG_ZH'] = { 'name': 'hww_WH3l_OSSF_mvaBDTG(Entry$,0)',
                        'range' : (10,-1,1),
                        'xaxis' : 'MVA_ZH',
                        'fold' : 3,
                       'linesToAdd' :  ['.L %s/src/PlotsConfigurations/Configurations/EFT/VBF/v7/Full2017/hww_WH3l_OSSF_mvaBDTG_old.C+' % os.getenv('CMSSW_BASE')]                       
                       }
'''
