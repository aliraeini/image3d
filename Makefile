build:
	#pip install nox
	nox -s build

build_pyi:
	pip install .
	pybind11-stubgen voxlib --output-dir src

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: build tests
