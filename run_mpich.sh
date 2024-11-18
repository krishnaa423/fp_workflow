#!/bin/bash

docker compose up -d fp-mpich-gpu-cc86-linux-amd64

docker container exec -it fp-mpich-gpu-cc86 /bin/bash