# Stability Prediction

## Overview

This repository holds a BV-BRC module for doing stability prediction.  It currently only supports the use of a custom version of [ThermoMPNN-D](https://github.com/Tinyman392/ThermoMPNN-D/tree/main) that adds CPU support.  Currently, the module will only run ThermoMPNN in CPU mode (as opposed to GPU mode).  There runtime for inference between CPU and GPU is not huge.  

## Running

Currently, this service can only be run if you are using an alpha version of the dev container since ThermoMPNN-D and it's python environments are currently not part of the production BV-BRC runtime.

``` bash
singularity shell --bind $(pwd) --bind /vol/patric3/cli/ubuntu-runtime/ /vol/patric3/production/containers/ubuntu-dev-093-12.sif
```

You'll also want to bind the location of the dev_container and any tests you'll be running in the singularity.  Once ThermoMPNN-D is added to the BV-BRC runtime (non-alpha), then the use of this container won't be necessary any longer.  

Like normal you'll want to source the user environment in the dev container.  

```bash
source /PATH/TO/dev_container/user-env.sh
```

After that you'll want to adjust your PATH variable too.

```bash
PATH=$PATH:/opt/patric-common/runtime/bin
```

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

