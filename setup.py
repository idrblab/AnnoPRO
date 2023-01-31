from pathlib import Path
from numpy.distutils.core import Extension, setup, numpy_cmdclass
from setuptools import find_packages
import versioneer

PACKAGE_DIR = "."
install_requires = Path("requirements.txt").read_text().split("\n")

setup(
    package_dir = {"": PACKAGE_DIR},
    install_requires = install_requires,
    packages = find_packages(where=PACKAGE_DIR),
    include_package_data=True,
    package_data = {
        'annopro.data_procession._libprofeatconfig': ['*.sdf', '*.dat']
    },
    ext_modules = [Extension(
        name="annopro.data_procession._libprofeat",
        sources=[f"{PACKAGE_DIR}/annopro/data_procession/_libprofeat.f"],
        extra_f77_compile_args=[
            ## Uncomment this arg if you using latest gfortran
            "-fallow-argument-mismatch",
            "-w"]
    )],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(numpy_cmdclass.copy())
)