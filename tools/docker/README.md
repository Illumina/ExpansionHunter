# ExpansionHunter docker utilities

### debug/

Docker based debug facilities, where the guest OS provides special
debug build and/or runtime diagnostics capabilities.

### deployment/

These are docker image instructions and scripts associated with
issuing a conventional binary release. These images typically
use an older OS (for wide compatibility) with a newer compiler (for
binary performance, new c++ features, etc.). These images are not
meant to assist users who want to run in docker.

