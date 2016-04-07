import os
import galaxy as gl
import pandas as pd
import warnings
"""
When you want to run the script in your PC, please change the path of files
"""


warnings.filterwarnings('ignore')
info = pd.read_csv('list.csv')
t = open('type2.txt', 'w')
# t.write('NAME,G,M,A,C\n')
for i in range(0, len(info)):
    # fg = info.ix[i]['NAME2']+'_g.fits'
    fr = info.ix[i]['NAME2']+'_r.fits'
    if os.path.exists('/home/fmxustc/Desktop/type2cut/'+fr):
        print(i, fr)
        a = gl.Galaxy('/home/fmxustc/Desktop/type2cut/'+fr)
        a.truncate()
        a.find_pollutions()
        a.eliminate_pollution()
        # a.calculate_parameters()
        # t.write(fr+',')
        b = a.galaxy_information['parameter']
        # t.write(str(b['G'])+','+str(b['M'])+','+str(b['A'])+','+str(b['C'])+'\n')
        # t.flush()
        os.system('mv truncate.fits /home/fmxustc/Desktop/galaxy/type2/truncate/'+fr)
        os.system('mv pollutions.fits /home/fmxustc/Desktop/galaxy/type2/pollution/'+fr)
        # os.system('mv treated.fits /home/fmxustc/Desktop/galaxy/type2/treated/'+fr)
        # os.system('mv galaxy.fits /home/fmxustc/Desktop/galaxy/type2/galaxy/'+fr)
        os.system('mv result.txt /home/fmxustc/Desktop/galaxy/type2/result/'+fr+'.txt')
t.close()
