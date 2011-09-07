from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name = "faster hgp",
    cmdclass = {"build_ext": build_ext},
    ext_modules = [Extension("hgp_cython",["hgp_cython.pyx"])],
    )
