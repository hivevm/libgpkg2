# Boost C++ Libraries

The Boost project provides free, peer-reviewed, portable C++ source libraries.

We emphasize libraries that work seamlessly with the C++ Standard Library. Boost libraries are designed to be widely useful and applicable across a broad spectrum of applications. The Boost license encourages both commercial and non-commercial use without requiring attribution for binary distributions.

Visit our project website at www.boost.org for more information and to [download](https://www.boost.org/users/download/) the current release.

-----

## Reducing Library Size

The **bcp** utility extracts subsets of Boost, making it valuable for:
- Boost authors who want to distribute their library independently
- Developers who need to include only specific Boost components with their application

The tool can also analyze your code's Boost dependencies and report the licenses associated with those components.

### Installation

**Linux/Unix:**
```bash
sudo apt install libboost-tools-dev
```

**Windows:**

On Windows, bcp must be compiled from source. Assuming your Boost distribution is extracted to `C:\development\Tools\libs\boost_1_XX_0\`:

1. Open a command prompt and navigate to the `boost_1_XX_0` root directory
2. Run `bootstrap.bat` (this creates `b2.exe` in the `tools\build` directory)
3. Run `tools\build\b2.exe tools\bcp` (this creates `bcp.exe` in the `dist\bin` directory)

For detailed Windows build instructions, see [this guide](https://medium.com/@biswa8998/building-c-boost-and-using-bcp-exe-f89881b2cc60).

### Usage

Basic command syntax:
```bash
bcp --boost=BOOST_HOME MODULES_OR_HEADERS OUTPUT_PATH
```
or
```bash
bcp --boost=BOOST_HOME --scan FILES_USING_BOOST OUTPUT_PATH
```

**Parameters:**
- `BOOST_HOME` — Directory containing the "boost" folder (optional, defaults to `./`)
- `MODULES_OR_HEADERS` — List of Boost modules and/or header files to extract
- `FILES_USING_BOOST` — List of your source files that include Boost headers
- `OUTPUT_PATH` — Destination directory for extracted files (**must exist before running**)

**Important:** All paths (BOOST_HOME, FILES_USING_BOOST, and OUTPUT_PATH) must be on the same physical drive, or bcp will terminate with an error.


## Example: Linux

The implementation expects the boost headers in the _gpkg/boost_ directory. Currently _libgpkg_ uses version 1.90.

The headers are already included in the sources. If an update is required following steps are necessary:

```sh
# Remove current headers
rm -rf gpkg/boost
rm -rf gpkg/libs

# Download and extract Boost
wget -P /tmp https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.bz2
tar xjpf /tmp/boost_1_90_0.tar.bz2 -C /tmp

### Use the development header to deploy Boost
bcp --boost=/tmp/boost_1_90_0/ --scan boost.hpp ./boost/
```


### Additional Resources

- [Official bcp documentation](http://www.boost.org/doc/libs/1_90_0/tools/bcp/doc/html/index.html)
- [Windows build guide](https://medium.com/@biswa8998/building-c-boost-and-using-bcp-exe-f89881b2cc60)