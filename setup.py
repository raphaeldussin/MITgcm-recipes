from numpy.distutils.core import setup, Extension

akima1d = Extension(name = 'MITgcm_recipes.akima1d',
                sources = ['MITgcm_recipes/f90/mod_akima_1d.f90'])

drown    = Extension(name = 'MITgcm_recipes.mod_drown_sosie',
                sources = ['MITgcm_recipes/f90/mod_drown_sosie.f90'])

setup(
    name = "MITgcm_recipes",
    version = "0.1",
    author = "Raphael Dussin",
    author_email = "rdussin@ldeo.columbia.edu",
    description = ("MITgcm tools for Real World Applications" ),
    license = "GPLv3",
    keywords = "ocean modeling, BGC",
    url = "",
    ext_modules = [akima1d, drown],
    packages=['MITgcm_recipes']
)
