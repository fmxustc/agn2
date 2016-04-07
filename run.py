import os
import galaxy as gl
import pandas as pd
import warnings
"""
When you want to run the script in your PC, please change the path of files
"""


warnings.filterwarnings('ignore')
info = pd.read_csv('list.csv')
t = open('type1r.csv', 'w')
t.write('NAME,G,M,A,C\n')
for i in range(0, len(info), 5):
    # fg = info.ix[i]['NAME2']+'_g.fits'
    fr = info.ix[i]['NAME1']+'_r.fits'
    if os.path.exists('/home/fmxustc/Desktop/type1cut/'+fr):
        print(i, fr)
        a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/'+fr)
        a.truncate()
        a.find_pollutions()
        a.eliminate_pollution()
        a.calculate_parameters()
        b = a.galaxy_information['parameter']
        for k in range(5):
            t.write(fr + ',')
            t.write(str(b['G'])+','+str(b['M'])+','+str(b['A'])+','+str(b['C'])+'\n')
            t.flush()
        os.system('mv truncate.fits /home/fmxustc/Desktop/galaxy/type1r/truncate/'+fr)
        os.system('mv pollutions.fits /home/fmxustc/Desktop/galaxy/type1r/pollution/'+fr)
        os.system('mv treated.fits /home/fmxustc/Desktop/galaxy/type1r/treated/'+fr)
        os.system('mv galaxy.fits /home/fmxustc/Desktop/galaxy/type1r/galaxy/'+fr)
        os.system('mv result.txt /home/fmxustc/Desktop/galaxy/type1r/result/'+fr+'.txt')
t.close()
