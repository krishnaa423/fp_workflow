from builtins import input
import numpy as np
import os

# # Enter the number of irr. q-points
# prefix = input('Enter the prefix used for PH calculations (e.g. diam)\n')

# # Enter the number of irr. q-points
# nqpt = input('Enter the number of irreducible q-points\n')
# try:
#   nqpt = int(nqpt)
# except ValueError:
#   raise Exception('The value you enter is not an integer!')

def get_nqpt(prefix):
    fname = './tmp/_ph0/' + prefix + '.phsave/control_ph.xml'

    fid = open(fname, 'r')
    lines = fid.readlines()
    # these files are relatively small so reading the whole thing shouldn't
    # be an issue
    fid.close()

    line_number_of_nqpt = 0
    while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]:
        # increment to line of interest
        line_number_of_nqpt += 1
    line_number_of_nqpt += 1  # its on the next line after that text

    nqpt = int(lines[line_number_of_nqpt])

    return nqpt

prefix='struct'
nqpt=get_nqpt(prefix)

os.system('mkdir save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp ./tmp/_ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r ./tmp/_ph0/'+prefix+'.phsave save/')
  else:
    os.system('cp ./tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('rm ./tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )