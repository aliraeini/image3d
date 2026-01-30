
help:
	@echo "Install prerequisite using commands:"
	@echo "  python -m venv ./.venv"
	@echo "  source .venv/bin/activate"
	@echo "  pip install pybind11-stubgen numpy  pytest"
	@echo
	@echo "Available target commands:"
	@echo "  make build_update # Install and update pybind11 stubs"
	@echo "  make install      # Install the package"
	@echo "  make test         # Run tests"

build_update:
	pip install .
	pybind11-stubgen image3d --output-dir src

install:
	pip install .

test:
	# .venv/bin/pip install pytest numpy
	python -m pytest

.PHONY: build tests
