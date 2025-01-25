from setuptools import setup, find_packages

setup(
    name="phenofun",
    version="0.1.0",
    author="Guilherme Azevedo",
    description="A CLI tool to calculate the probability of fixation of differences in a hypothetical nuclear locus that controls phenotype under neutral divergence.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ghfazevedo/phenofun",
    packages=["phenofun"],
    package_dir={"phenofun": "src/"},  
    include_package_data=True,
    install_requires=[
        'dendropy'
    ],
    entry_points={
        'console_scripts': [
            'phenofun = phenofun.phenofun:main',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)

