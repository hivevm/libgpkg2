# Boost C++ Libraries

The Boost project provides free peer-reviewed portable C++ source libraries.

We emphasize libraries that work well with the C++ Standard Library. Boost
libraries are intended to be widely useful, and usable across a broad spectrum
of applications. The Boost license encourages both commercial and non-commercial use
and does not require attribution for binary use.

The project website is www.boost.org, where you can obtain more information and
[download](https://www.boost.org/users/download/) the current release.

-----

## Reduce Library size

The *bcp* utility is a tool for extracting subsets of Boost, it's useful for Boost 
authors who want to distribute their library separately from Boost, and for 
Boost users who want to distribute a subset of Boost with their application.

*bcp* can also report on which parts of Boost your code is dependent on, and what l
icences are used by those dependencies.

The command to install it under *unix* is:
```
sudo apt install libboost-tools-dev
```
On Windows the tool bcp must be compiled from the boost distribution. 
Assuming that the dir where the boost distribution has been extracted is "c:\development\Tools\libs\boost_1_78_0\", The steps are:

- open a cmd, go in the root dir "boost_1_78_0"
- issue the command 'bootstrap.bat' After completion in the directory 'tools\build' now there is the executable 'b2.exe'
- issue the command 'tools\build\b2.exe tools\bcp'
After completion in the directory 'dist\bin' now there is the executable 'bcp.exe'

The most common usages are quite simple:
```
bcp --boost=BOOST_HOME MODULES_OR_HEADERS OUTPUT_PATH
or
bcp --boost=BOOST_HOME --scan FILES_USING_BOOST OUTPUT_PATH
```
where
* BOOST_HOME is the folder containing one more "boost" directory (in this case the location of this readme file). This parameter is **optional** and defaults to "./"
* MODULES_OR_HEADERS is a list of boost modules and/or single boost header files we want to keep
* FILES_USING_BOOST is a list of [non-boost] files which include boost headers
* OUTPUT_PATH is the directory where the necessary files are copied (**this path must exist**)

Note: BOOST_HOME, FILES_USING_BOOST and OUTPUT_PATH must be on the same physical disk, otherwise the BCP tools terminates with an error.

### Linux example
~~~shell
// Scan boost.hpp and copy the dependencies in the current folder.

rm -rf gpkg/boost
rm -rf gpkg/libs
wget https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2
tar xvjpf boost_1_88_0.tar.bz2
bcp --boost=boost_1_88_0/ --scan boost.hpp ./gpkg/
~~~


For more information read the [documentation](http://www.boost.org/doc/libs/1_88_0/tools/bcp/doc/html/index.html).
For more information about building bcp on Windows see [this post](https://medium.com/@biswa8998/building-c-boost-and-using-bcp-exe-f89881b2cc60)
