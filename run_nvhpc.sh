#!/bin/bash

docker compose up -d fp-nvhpc-cc86-linux-amd64

docker container exec -it fp-nvhpc /bin/bash