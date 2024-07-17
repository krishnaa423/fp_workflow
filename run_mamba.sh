#!/bin/bash

docker compose up -d mamba-mpich-linux-amd64

docker container exec -it mamba /bin/bash