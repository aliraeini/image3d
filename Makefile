build:
	#pip install nox
	nox -s build

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: build tests
