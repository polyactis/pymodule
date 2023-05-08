#!/usr/bin/env python3g

import os
import sys
from setuptools import setup, find_packages

src_dir = os.path.dirname(__file__)

install_requires = ['h5py', 'matplotlib', "numpy", 
    'scipy', "sqlalchemy", 'statsmodels', 'tables',
    'future', 'future-fstrings']

if sys.version_info.major>2:
    #biopython, Pegaflow is Python3 only
    install_requires.append("biopython")
    install_requires.append("pegaflow")

with open(os.path.join(src_dir, 'README.md')) as readme_file:
    README = readme_file.read()
#
# Create Manifest file to exclude tests, and service files
#
def create_manifest_file():
    f = None
    try:
        f = open('MANIFEST.in', 'w')
        f.write(u"global-exclude *.py[cod]\n")
    finally:
        if f:
            f.close()

#
# Install conditional dependencies
#
def setup_installer_dependencies():
    global install_requires

def find_package_data(dirname):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items

    items = find_paths(dirname)
    return [path.replace(dirname, "") for path in items]




if __name__ == '__main__':
    create_manifest_file()
    setup_installer_dependencies()
    setup(
        name="Palos",
        version="0.1.35",
        author="Yu S. Huang",
        author_email="polyactis@gmail.com",
        description="Misc Python modules developed and used by the yfish group",
        long_description_content_type="text/markdown",
        long_description=README,
        license="Apache2",
        url="https://github.com/polyactis/pymodule",
        python_requires=">=2.7",
        classifiers=[
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Operating System :: Unix",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering",
            "Topic :: Utilities",
        ],
        packages=find_packages(include=['palos', 'palos.*'], exclude=['.test*']),
        package_data={
            # If any package contains *.sh files, include them:
            "": ["*.sh", "*.md"],
        },
        include_package_data=True,
        zip_safe=False,
        install_requires = install_requires,
    )
