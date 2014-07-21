from distutils.core import setup, Extension
import glob
#from setuptools import find_packages

hsmsources = glob.glob("mapy/source_files/libhsm/*.c")
hsm_module = Extension('libhsm', 
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
    ext_modules = [hsm_module]
)