#! /bin/bash

cd /home/ubuntu/r-projekt/ || exit
git fetch origin main
git pull
# git merge -X theirs
# git add .
# git commit -m "Sync @ $(date)"
# git push
