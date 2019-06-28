# correlation-ws

This repository contains code for loading correlation and expression data, and retrieving it from a web service.

*NOTE:* to build this code, you will need the Java HDF5 library. You can get this from https://portal.hdfgroup.org/display/support
To install it, run this command: 

```bash
mvn install:install-file -Dfile=./jarhdf5-1.10.5.jar -DgroupId=hdf5 -DartifactId=jarhdf5 -Dversion=1.10.5 -Dpackaging=jar
```