from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()


setup(
    name="fragmentembedding",
    version="0.0.1",
    author="Miquel Duran-Frigola",
    author_email="miquel@ersilia.io",
    url="https://github.com/github.com/fragment-embedding",
    description="Fully Functionalized Fragment Embedding for the Ligand Discovery project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.10",
    install_requires=install_requires,
    packages=find_packages(exclude=("utilities")),
    entry_points={"console_scripts": []},
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    keywords="artificial-intelligence chemistry-embedding fragment-based-drug-discovery",
    project_urls={
        "Source Code": "https://github.com/ligand-discovery/fragment-embedding"
    },
    include_package_data=True,
)
