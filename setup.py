import os
from distutils.sysconfig import get_python_lib

# Set path
print('>>> Set path...')

pwd = os.getcwd()
lib_path = get_python_lib()
pth_file = os.path.join(lib_path, 'dyntripy.pth')
with open(pth_file, 'w') as f:
    f.write(pwd)

# Install packages
print('>>> Install packages...')
try:
    os.system('pip install -r requirements.txt')
except BaseException:
    print('Install failed!')
