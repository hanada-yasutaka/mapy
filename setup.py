from distutils.core import setup
from setuptools import find_packages
#PKG = ['mapy','mapy.QtMapy','mapy.QtMapy.Design']
pkgs = find_packages()

setup(
    name = 'mapy',
    version='0.1',
    description='Calculator of symplectic system',
    author='Yasutaka Hanada',
    author_email='hanada.yasutaka@gmail.com',
    url='http://www.abc.com',
    packages=pkgs,
)