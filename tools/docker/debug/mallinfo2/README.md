
Scripts to build/deploy EH on a glibc 2.33 OS in order to use new mallinfo2 API.
Note this isn't completely push-button - a few edits to external libraries may be
required to complete the build. 

To generate docker image and EH binary, use:

    setup-docker-image-and-build.bash

To further convert this image to singularity:

    sudo singularity build mallinfo2.sif docker-daemon://ehdebug:mallinfo2

