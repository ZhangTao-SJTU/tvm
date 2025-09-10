#!/bin/bash

/usr/local/bin/singularity shell $1 <<EOT
cd $(pwd)
$2 ${@:3}
EOT

exit 0
