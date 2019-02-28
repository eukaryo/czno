#!/bin/sh

python3 random_dataset_generator.py isolen > dataset-iso-len.txt
python3 iso_len_time_experiment.py dataset-iso-len.txt > result-isolen-time.txt
python3 random_dataset_generator.py isodist > dataset-iso-dist.txt
python3 iso_dist_time_experiment.py dataset-iso-dist.txt > result-isodist-time.txt
python3 s151Rfam-localminima-generator.py > s151-localminima-dataset.txt
tar cvJf s151-localminima-dataset.txt.tar.xz s151-localminima-dataset.txt
tar xvJf s151-localminima-dataset.txt.tar.xz
python3 directpath-comparison-experiment.py > directpath-result.txt