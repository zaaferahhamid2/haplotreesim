"""
Setup script for HaploTreeSim package.
"""

from setuptools import setup, find_packages

setup(
    name="haplotreesim",
    version="0.1.0",
    description="A controlled, end-to-end benchmark for low-pass scDNA-seq with haplotype-resolved CNA ground truth",
    author="name",
    author_email="email",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
        ],
        "analysis": [
            "pandas>=1.3.0",
            "matplotlib>=3.4.0",
            "seaborn>=0.11.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
