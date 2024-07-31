# Time-Optimal Route Planning for Non-Linear Recharging Electric Vehicles on Road Networks
Thanks for all of you for your interest in our work.

This repository is the experiment code for paper: "Time-Optimal Route Planning for Non-Linear Recharging Electric Vehicles on Road Networks".

## Brief Introduction
The cpp file "Query-exact.cpp" is the exact algorithm for route planning of any two vertices on road networks, while the "Query-aprox.cpp" is the approximate algorithm.

## Road Networks
You can download the road networks from [9th DIMACS Implementation Challenge - Shortest Paths](https://www.diag.uniroma1.it/challenge9/download.shtml). Besides, since in our paper, each edge has two properties: travel time and electricity consumption, you can use "graph_process.cpp" to generate a graph that fits our papers.

## Complie & Run
After downloading the dataset, you can use `g++ -O3 Query-exact.cpp -std=c++11 -o exact.out` to compile the code and generate a "exact.out" file. You can use `./exact.out` to run the output file.

## Contact
- Qinzhou Xiao: qinzhou.xiao@stu.ecnu.edu.cn
- Yu Shao: yushao@stu.ecnu.edu.cn
- Peng Cheng: pcheng@sei.ecnu.edu.cn
