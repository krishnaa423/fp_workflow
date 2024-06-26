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

prefix='struct'
nqpt=1

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