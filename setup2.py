from distutils.core import setup, Extension
from os import getcwd
cwd = getcwd()
module1 = Extension('KM_TreeAlign', sources =['alignment.c','suffix_tree.c'],extra_compile_args = ['-w'],extra_link_args = ['-lm'])
setup(name = 'KM_TreeAlign', version = '1.0',description = 'Suffix Tree based overlap locator',ext_modules = [module1])
module2 = Extension('kvqsplitForClust', sources =['kamivq.c','fvqe.c'],extra_compile_args = ['-w'],extra_link_args = ['-lm'])
setup(name = 'kvqsplitForClust', version = '1.0',description = 'Linde, Buzo, and Gray Clustering Algorithm by K. Sayood',ext_modules = [module2])
