#!/bin/bash

execute_jobs() {
    for file in $(find . -type f -name "job.sh"); do
        echo "submitting: $file"
        qsub "$file"
    done
}

execute_jobs