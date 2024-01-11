from setuptools import setup, find_packages
from genebe import version

setup(
    name="genebe",
    version=version.__version__,
    packages=find_packages(),
    install_requires=["pymmh3", "tinynetrc", "pandas", "requests", "tqdm", "jsonlines"],
    extras_require={
        "cpp": ["mmh3", "cyvcf2"],
    },
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Piotr Stawinski",
    description="GeneBe Client: A user-friendly system for annotating genetic variants",
    url="https://github.com/pstawinski/pygenebe",
    entry_points={
        "console_scripts": [
            "genebe=genebe.entrypoint:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
