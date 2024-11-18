#!/bin/bash

docker compose up -d fp-testing-amd64

docker container exec -it fp-testing-amd64 /bin/bash