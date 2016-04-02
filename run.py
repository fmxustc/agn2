import os
import galaxy as gl
import pandas as pd
import warnings

warnings.filterwarnings('ignore')
info = pd.read_csv('list.csv')
t = open('type2.txt', 'w')
t.write('NAME,G,M,A,C\n')
for i in range(437, len(info)):
    # fg = info.ix[i]['NAME2']+'_g.fits'
    fr = info.ix[i]['NAME2']+'_r.fits'
    if os.path.exists('/Users/franky/Desktop/galaxy_fits/type2cut/'+fr):
        print(i, fr)
        a = gl.Galaxy('/Users/franky/Desktop/galaxy_fits/type2cut/'+fr)
        a.truncate()
        a.find_pollutions()
        a.eliminate_pollution()
        a.calculate_parameters()
        t.write(fr+',')
        b = a.galaxy_information['parameter']
        t.write(str(b['G'])+','+str(b['M'])+','+str(b['A'])+','+str(b['C'])+'\n')
        t.flush()
        os.system('mv truncate.fits ./truncate/'+fr)
        os.system('mv background.fits ./check/'+fr)
        os.system('mv treated.fits ./treated/'+fr)
        os.system('mv galaxy.fits ./galaxy/'+fr)
        os.system('mv result.txt ./result/'+fr+'.txt')
t.close()
