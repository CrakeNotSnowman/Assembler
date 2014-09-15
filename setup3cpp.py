from distutils.core import setup, Extension

# the c++ extension module
module3 = Extension('JBGenSuffTree', sources =['EdgeBag.cpp', 'Edge.cpp', 'GeneralizedSuffixTree.cpp', 'Node.cpp', 'NodeLabel.cpp', 'VkObject.cpp', 'Assembler.cpp', 'fns.cpp', 'main.cpp'],extra_compile_args = ['-O2', '-std=c++0x', '-static', '-ggdb'],extra_link_args = ['-lm'])
setup(name = 'JBGenSuffTree', version = '1.0',description = 'A Generalized Suffix Tree based Sequence Assembler',ext_modules = [module3])

