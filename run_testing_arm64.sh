#!/bin/bash

docker compose up -d fp-testing-arm64

docker container exec -it fp-testing-arm64 /bin/bash