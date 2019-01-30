#!/usr/bin/env bash

kubectl -n kube-system describe secret $(kubectl -n kube-system get secret | grep admin | awk '{print $1}')