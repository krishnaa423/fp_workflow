#!/bin/bash

docker compose up -d fp-testing

docker container exec -it fp-testing /bin/bash