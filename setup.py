import setuptools
import TFEA

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tfea",
    version=TFEA.__version__, #Version read from __init__.py
    python_requires=">=3",
    description="Transcription Factor Enrichment Analysis",
    url="https://github.com/Dowell-Lab/TFEA.git",
    author="Jonathan Rubin",
    author_email="jonathan.rubin@colorado.edu",
    license="GPL-3.0",
    packages=setuptools.find_packages(),
    package_data={"": ["test/test_files/*", "*sbatch"]},
    long_description=long_description,
    long_description_content_type="text/markdown",
    scripts=["bin/TFEA", "bin/TFEA-annotate", "bin/TFEA-simulate"],
    install_requires=[
        "matplotlib>=3.1.1",
        "scipy",
        "numpy",
        "pybedtools",
        "htseq",
        "psutil",
        "ujson",
        "pysam==0.15.2"
    ],
    zip_safe=False,
)
