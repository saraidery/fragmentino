from setuptools import setup

extra_requirements = {
	"tests": ["pytest", "coverage", "pytest-cov"]
}

setup(
	name="framol",
	author="Sarai D. Folkestad",
	description="Description of package",
	install_requires=["numpy", "pandas"],
	extras_require=extra_requirements
)
