## Installation

### macOS
You need to install `libomp` required by the project:
```bash
brew install libomp ninja

# activate a venv: 
# python -m venv .venv
# .venv/bin/activate
pip install --use-feature=in-tree-build .

# to test your installation, run:
pip install pytest
python -m pytest
```

### Windows
You need to install `zlib` and `libtiff`. We recommend using `vcpkg`:
```powershell
git clone https://github.com/microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg install zlib:x64-windows
```
And set `CMAKE_TOOLCHAIN_FILE` environment variable to `[vcpkg root]/scripts/buildsystems/vcpkg.cmake`.
