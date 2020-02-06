import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tfea",
    version="1.0.1",
    python_requires=">=3",
    description="Transcription Factor Enrichment Analysis",
    url="https://github.com/jdrubin91/TFEA.git",
    author="Jonathan Rubin",
    author_email="jonathan.rubin@colorado.edu",
    license="CU Boulder Dowell Lab",
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
    ],
    zip_safe=False,
)
