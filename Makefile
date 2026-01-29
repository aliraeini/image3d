build_update:
	# nox -s build
	pip install .
	pybind11-stubgen image3d --output-dir src

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: build tests
