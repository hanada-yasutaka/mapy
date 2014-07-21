from distutils.core import setup, Extension
import glob
#from setuptools import find_packages
#import os
#PKG = ['mapy','mapy.QtMapy','mapy.QtMapy.Design']
#pkgs = find_packages()
#pkgs.append('mapy.shared')
#print(pkgs)

#os.system("cd mapy/; bash make.sh")

#cfiles = ["wrapper_HusimiRep.c", "CoherentState.c", "HusimiRep.c", "c_complex.c"]
hsmsources = glob.glob("mapy/source_files/libhsm/*.c")
headers = glob.glob("mapy/source_files/libhsm/*.h")
#hsmsources=['mapy/source_files/libhsm/' + cf for cf in cfiles ]
#hsmsources=['mapy/source_files/libhsm/' + cf for cf in cfiles ]
module1 = Extension('libhsm', 
#    include_dirs = ['mapy/source_files/libhsm/'],
#    library_dirs = ['mapy/source_files/libhsm/'],
#    headers=headers,
    sources = hsmsources)


setup(
    name = 'mapy',
    version='0.3',
    description='Calculator of symplectic system',
    author='Yasutaka Hanada',
    author_email='hanada.yasutaka@gmail.com',
    url='http://www.abc.com',
    packages=['mapy','mapy.Classical','mapy.Quantum','mapy.ctypes_wrapper',\
              'mapy.QtMapy','mapy.QtMapy.Design'],
    ext_package='mapy.shared',
    ext_modules = [module1]
)