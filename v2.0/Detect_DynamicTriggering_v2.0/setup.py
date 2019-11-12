import os
from distutils.sysconfig import get_python_lib

# Set path
print ('>>> Set path...')

pwd=os.getcwd()
lib_path=get_python_lib()
pth_file=os.path.join(lib_path,'dyntrigger.pth')
with open(pth_file,'w') as f:
    f.write(pwd)

# Install packages
print ('>>> Install packages...')

os.system('pip install numpy')
os.system('pip install pandas')
os.system('pip install scipy')
os.system('pip install obspy')
os.system('pip install seaborn')
os.system('pip install statsmodels')

#
seaborn_path=os.path.join(lib_path,'seaborn')
os.system('cp dyntrigger/utils/distributions.py %s'%seaborn_path)

print('Finish! Welcome to dyntrigger!')
