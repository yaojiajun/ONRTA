# ONRTA: Non-rejection aware Online Task Assignment in Spatial Crowdsourcing
This repository stores the source code of the solutions to the problem called ONRTA.

## Usage of the algorithms
Environment
gcc/g++ version: 7.4.0

OS: Ubuntu

Compile the algorithms
cd algorithm && make all

## Run the algorithms
./Random ./synthetic/workers/100_3000_1_6_0.5_6_10/data_00.txt

./Greedy ./synthetic/workers/100_3000_1_6_0.5_6_10/data_00.txt

./ONRTA-RT ./synthetic/workers/100_3000_1_6_0.5_6_10/data_00.txt

./ONRTA-Base ./synthetic/workers/100_3000_1_6_0.5_6_10/data_00.txt

./ONRTA-OP ./real/EverySender_cap1/800/data_00.txt

./ONRTA-Greedy ./real/EverySender_cap1/800/data_00.txt

./OPT ./real/EverySender_cap1/800/data_00.txt

## Description of the datasets

Environment
Python: 2.7

Synthetic dataset
dataset/synthetic: a sample of our synthetic dataset (#2)

dataset/genDataSynthetic.py: a script to generate the synthetic datasets

Please refer to genDataSynthetic.py for the format of the dataset.

Real dataset
dataset/real/Didi*: includes the datasets of Didi

dataset/real/EverySender*: incldues the datasets of EverySender

Please refer to the source code for the format of the dataset.

## Contact
If you have any questions or concerns, please raise an issue or email: 201810102795@mail.scut.edu.cn
