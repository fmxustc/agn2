import galaxy
import warnings
"""
When you want to run the script in your PC, please change the path of files
"""


warnings.filterwarnings('ignore')


a = galaxy.Galaxy('/home/fmxustc/Desktop/type1cut/J080538.66+261005.5_r.fits')
a.truncate()
a.find_pollutions()
a.eliminate_pollution()
a.calculate_parameters()
a.show_edge_of_galaxy()

# r = open('type12_mc.txt', 'r')
# w = open('list.csv', 'w')
# w.write('NAME1,RA1,DEC1,Z1,NAME2,RA2,DEC2,Z2\n')
# for line in r.readlines()[1:]:
#     item = line.split()
#     for i in range(len(item)-1):
#         w.write(item[i]+',')
#     w.write(item[len(item)-1]+'\n')





