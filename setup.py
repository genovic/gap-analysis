from setuptools import find_packages, setup

setup(
    name='gap_analysis',
    version='0.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'find_gaps = gap_analysis.find_gaps:main',
            'find_coverage = gap_analysis.find_coverage_at_position:main',
            'filter_genes = gap_analysis.filter_genes:main',
        ]
    }
)
