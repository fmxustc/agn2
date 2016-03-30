import galaxy as gl
import warnings

warnings.filterwarnings('ignore')


a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/J083732.70+284218.7_r.fits')
a.truncate()
a.find_pollutions()
a.eliminate_pollution()
a.calculate_parameters()



