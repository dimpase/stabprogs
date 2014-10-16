from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules=[
   Extension("wlcy",
             ["wlcy.pyx"],
	     library_dirs = ['.'],
             libraries=["wl"]) # Unix-like specific
]
setup(
   name = "Demos",
   cmdclass = {"build_ext": build_ext},
   ext_modules = ext_modules
)
