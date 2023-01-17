from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from setuptools import find_packages
PACKAGE_DIR = "."
setup(
    package_dir = {"": PACKAGE_DIR},
    install_requires = ["numpy"],
    packages = find_packages(where=PACKAGE_DIR),
    include_package_data=True,
    package_data = {
        'annopro.data': ['*'],
        'annopro.data_procession._libprofeatconfig': ['*.sdf', '*.dat'],
        'annopro.model_param': ['*.h5']
    },
    ext_modules = [Extension(
        name="annopro.data_procession._libprofeat",
        sources=[f"{PACKAGE_DIR}/annopro/data_procession/libprofeat.f"],
        extra_f77_compile_args=[
            "-fallow-argument-mismatch",
            "-w"]
    )]
)