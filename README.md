GPStudio
========
Repository for paper *Knowledge-Based Trajectory Completion from Sparse GPS Samples*, by Yang Li, Yangyan Li, Dimitrios Gunopulos, and Leonidas Guibas, published in Proceedings of ACM SIGSPATIAL 2016.


Requirements and Dependencies
------------
Tested on Ubuntu 14.04 and 16.04
Minimum requirement: 8G ram

GPStudio depends on the following libraries:
* gflags
* opencv
* boost
* protobuf
* QT5
* openscenegraph
* geographiclib

To Compile
----------

1. Configure

Go to build directory and call cmake

    cd build
    cmake ..

Specify include and cmake paths if nessecessary.
On Yang's workstation:

    cmake -DCMAKE_PREFIX_PATH=/opt/Qt5.6.0/5.6/gcc_64 ..

On Yangyan's workstation:

```cmake
cmake -DBoost_INCLUDE_DIR=/home/yang/boost/include -DBoost_LIBRARY_DIR=/home/yang/boost/lib -DCV_INCLUDE_DIR=/home/yang/opencv/include -DCV_LIB_DIR=/home/yang/opencv/build/lib ..
```


2. Make protobuf and update source file

```
make proto
./setup_proto.sh
```

3. Make executable
```
make GPStudio
```

To Run
------

## GUI mode
```
bin/GPStudio
```
## Command line mode
```
bin/GPStudio --interactive=false --trace_fname=<gps_data> --bbox_fname=<bound_box>
```
Alternatively, use `guiless.sh` to invoke GPStudio and process output.

See `doc/DEMO.md` for the tutorial.

Known issues
------------
- Invalid pointer error when exiting the GUI.
- Not optimized for performance and memory
