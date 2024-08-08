#!/bin/bash
salloc --account=m3571 --qos=interactive --job-name=struct_job --constraint=gpu --nodes=4 --time=01:45:00
