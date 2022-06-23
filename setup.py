from setuptools import setup, find_packages


setup(
    name = "ForConX",
    version = "0.1",
    packages = find_packages(),

    # requires dependencies

    author = "Christian Schroeder",
    author_email = "christian.schroeder@univie.ac.at",

    # license = "GPLv3",
    
    entry_points={
        'console_scripts': [
            'forconx = ForConX.forconx:main',
        ]
    }
)
