build_update:
	nox -s build
	pip install .
	pybind11-stubgen voxlib --output-dir src

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: build tests
