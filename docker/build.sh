#!/bin/bash

docker build -t uncalled:latest --build-arg git_checkout=$(git rev-parse --short HEAD) --build-arg git_url=$(git config --get remote.origin.url) .
