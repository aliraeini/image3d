# Archived: see [image3kit](https://github.com/image3kit/image3kit) instead

## Image3d 

**EXPERIMENTAL: not ready for public use**

The `image3d` is Python packaging of a collection of C++ libraries for processing
and analysing 3D images such as those obtained using X-ray micro-tomography.

## Installation

### Linux / pip

Prerequisites: Cmake and C++ development tools.

```bash
pip install git+https://github.com/aliraeini/image3d.git
```

Note: This repository will most likely moved and the above address will stop working.

### macOS
Prerequisites: `libomp` and `ninja` build system:

```bash
brew install libomp ninja
export OpenMP_ROOT=$(brew --prefix libomp)

# activate a venv: 
# python -m venv .venv
# .venv/bin/activate
pip install .[test]

# To test your installation, run:
pytest
```

### Windows
You need to install `zlib` e.g. using `vcpkg`:

```powershell
git clone https://github.com/microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg install zlib:x64-windows

$env:CMAKE_TOOLCHAIN_FILE = "$PWD/vcpkg/scripts/buildsystems/vcpkg.cmake"
$env:SKBUILD_CMAKE_ARGS = "-DCMAKE_TOOLCHAIN_FILE=$env:CMAKE_TOOLCHAIN_FILE"

pip install .[test]
pytest
```

## Licenses

License is subject to change.

This package uses several third-party codes, see [pkgs/](pkgs/) for details.

