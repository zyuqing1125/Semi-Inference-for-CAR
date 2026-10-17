#!/bin/bash


for i in {1..1000}

do
    sbatch main_job_content.sh $i
done
