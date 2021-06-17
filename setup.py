from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(['./drive/pre_shared_segments_analysis_scripts/shared_segment_detection/combine_output/*.pyx'], language_level = "3"))