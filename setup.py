from pathlib import Path
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
from setuptools import find_packages
import sys

PACKAGE_DIR = "."
install_requires = Path("requirements.txt").read_text().split("\n")

if sys.platform == "linux":
    install_requires = ["fasta" if r.startswith("fasta") else r for r in install_requires]

setup(
    package_dir = {"": PACKAGE_DIR},
    install_requires = install_requires,
    packages = find_packages(where=PACKAGE_DIR),
    include_package_data=True,
    package_data = {
        'annopro.data': ['*'],
        'annopro.data_procession._libprofeatconfig': ['*.sdf', '*.dat'],
        'annopro.model_param': ['*.h5']
    },
    ext_modules = [Extension(
        name="annopro.data_procession._libprofeat",
        sources=[f"{PACKAGE_DIR}/annopro/data_procession/_libprofeat.f"],
        extra_f77_compile_args=[
            "-fallow-argument-mismatch",
            "-w"]
    )],
    python_requires=">=3.8",
    entry_points=dict(console_scripts=["annopro = annopro:console_main"])
)