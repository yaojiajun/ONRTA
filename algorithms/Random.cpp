/**
	\author: Trasier
	\date: 2017.04.02
*/
#include <bits/stdc++.h>
using namespace std;
//#pragma comment(linker,"/STACK:102400000,1024000")
#include "monitor.h"
#include "input.h"
#include <algorithm>

enum rule_t {
	worker, task
};

union W_un {
	double cost;
	double pay;
	double rate;
};
	
struct node_t {
	rule_t type;		// 0: task, 1: worker
	pair<double, double> loc;	// location
    float vel;			// velocity
	int flow;			// flow
	double rad;			// radius	
	W_un cost;			// cost
    int id;
	int begTime, endTime;	// time interval

	void print() {
		if (type == worker)
		 	printf("w: loc = (%.2lf, %.2lf), rad = %.2lf, cap = %d, time = (%d, %d), ratio = %.2lf\n",
		 			loc.first, loc.second, rad, vel, begTime, endTime, cost.rate);
		else
			printf("t: loc = (%.2lf, %.2lf), time = (%d, %d), pay = %.2lf\n",
		 			loc.first, loc.second, begTime, endTime, cost.pay);
	}
};

typedef long long LL;
int n, m, sumC;
double p=0.04;
int matched_worker=0;
double umax;
double utility;
int usedMemory;
vector<string> task_done;

void init(int taskN, int workerN, double Umax, int SumC) {
	n = workerN;
	m = taskN;
	umax = Umax;
	sumC = SumC;
	utility = 0;
	usedMemory = 0;
}

void nextSeq(ifstream& fin, node_t& nd) {
	int timeId;
	string stype;

	fin >> nd.id>>nd.begTime >> stype;
	if (stype[0] == 'w') {
		nd.type = worker;
		fin >> nd.loc.first >> nd.loc.second >> nd.rad >> nd.vel >> nd.endTime >> nd.cost.rate;
		nd.endTime += nd.begTime;
        //nd.vel = 2;
	} else {
		nd.type = task;
		fin >> nd.loc.first >> nd.loc.second >> nd.endTime >> nd.cost.pay;
		nd.endTime += nd.begTime;
	}

	nd.flow = 0;
}

inline double Length(pair<double,double> pa, pair<double,double> pb) {
	return sqrt( (pa.first-pb.first)*(pa.first-pb.first) + (pa.second-pb.second)*(pa.second-pb.second) );
}

inline double Length2(pair<double,double> pa, pair<double,double> pb) {
	return (pa.first-pb.first)*(pa.first-pb.first) + (pa.second-pb.second)*(pa.second-pb.second);
}


inline bool satisfyLoc(const node_t& worker, const node_t& task) {
	// 4. condition of location
	if (Length2(worker.loc, task.loc) > worker.rad * worker.rad)
		return false;
	
	return true;
}

inline bool satisfyCap(const node_t& worker, const node_t& task) {
	// 2&3. capacity of worker & task
	if (1<=worker.flow || 1<=task.flow)
		return false;
	return true;
}

inline bool satisfyTime(const node_t& worker, const node_t& task) {
	// 1. condition of deadline
	//if (!(worker.begTime<=task.endTime && task.begTime<=worker.endTime))
    double aa=Length(worker.loc, task.loc);
    int a=aa/worker.vel+1;

    if (!(worker.begTime+a<task.endTime && worker.begTime+a<worker.endTime && task.begTime<worker.endTime))
        return false;

	return true;
}

inline bool satisfyCost(const node_t& worker, const node_t& task) {
    // 1. condition of deadline
    //if (!(worker.begTime<=task.endTime && task.begTime<=worker.endTime))
    if (task.cost.rate * worker.cost.pay-(Length(worker.loc, task.loc)/worker.vel)*p<=0)
        return false;
    return true;
}

bool satisfy(const node_t& worker, const node_t& task) {
	return satisfyCap(worker, task) && satisfyTime(worker, task) && satisfyLoc(worker, task)&& satisfyCost(worker, task);
	//return satisfyCap(worker, task) && satisfyLoc(worker, task);
}

inline double calcCost(const node_t& task, const node_t& worker) {
    double cost=(Length(worker.loc, task.loc)/worker.vel)*p;
    return task.cost.pay * worker.cost.rate-cost;
}

vector<string> split(const string &s, const string &seperator){
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while(i != s.size()){
        //找到字符串中首个不等于分隔符的字母；
        int flag = 0;
        while(i != s.size() && flag == 0){
            flag = 1;
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[i] == seperator[x]){
                    ++i;
                    flag = 0;
                    break;
                }
        }
        //找到又一个分隔符，将两个分隔符之间的字符串取出；
        flag = 0;
        string_size j = i;
        while(j != s.size() && flag == 0){
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[j] == seperator[x]){
                    flag = 1;
                    break;
                }
            if(flag == 0)
                ++j;
        }
        if(i != j){
            result.push_back(s.substr(i, j-i));
            i = j;
        }
    }
    return result;
}


int chosenNextTask(const vector<node_t>& tasks, const node_t& worker) {
	int taskN = tasks.size();
	double tmpCost;
	double mxCost = 0.0;
	int ret = -1;
    vector<int> temp;
    for (int i = 0; i < taskN; ++i)
    {
        temp.push_back(i);
    }
    random_shuffle(temp.begin(), temp.end());
    while(taskN--){
        tmpCost = calcCost(tasks[temp[taskN]], worker);
        if (satisfy(worker, tasks[temp[taskN]]) && tmpCost>0) {
            return ret = temp[taskN];
        }
    }

	return ret;
}

int chosenNextWorker(const vector<node_t>& workers, const node_t& task) {
	int workerN = workers.size();
	double tmpCost;
	int ret = -1;
    vector<int> temp;
    for (int i = 0; i < workerN; ++i)
    {
        temp.push_back(i);
    }
    random_shuffle(temp.begin(), temp.end());
    while(workerN--){
        tmpCost = calcCost(task, workers[temp[workerN]]);
        if (satisfy(workers[temp[workerN]], task) && tmpCost>0) {
            return ret = temp[workerN];
        }
    }
	return ret;
}

void addOneMatch(node_t& task, node_t& worker) {
	// add cost to utility
	utility += calcCost(task, worker);
	// update the capacity of task & worker
	++task.flow;
	++worker.flow;
}

void Pure_Greedy(ifstream& fin, int seqN) {
	node_t node;
	vector<node_t> tasks, workers,workingworker, tworkingworker;
    vector<node_t> temptasks,tempworkers;
	int taskId, workerId;
    node_t nd;
    vector<node_t> Sequence;
    vector<node_t> batch;

    while (seqN--) {
        nextSeq(fin, node);
        Sequence.push_back(node);
    }
    for (int id = 0; id < Sequence.size(); ++id) {
        batch.push_back(Sequence[id]);
        node_t nd = Sequence[id];
        if (id == Sequence.size() - 1 || nd.begTime != Sequence[id + 1].begTime) {
            for(int i=0; i< workingworker.size(); ++i){
                if(workingworker[i].begTime<=nd.begTime&&workingworker[i].endTime>=nd.begTime){
                    workingworker[i].flow=0;
                    batch.push_back(workingworker[i]);
                }
            }
            tworkingworker.clear();
            for (int t = 0; t < workingworker.size(); t++) {
                if (workingworker[t].flow >0 &&
                    workingworker[t].endTime > node.begTime)//&&tasks[t].endTime>=node.begTime
                    tworkingworker.push_back(workingworker[t]);
            }
            workingworker.clear();
            for (int t = 0; t < tworkingworker.size(); t++) {
                workingworker.push_back(tworkingworker[t]);
            }
            for (int ii = 0; ii < batch.size(); ++ii) {
                node = batch[ii];
                workerId = taskId = -1;
                if (node.type == task) { // node is task
                    taskId = tasks.size();
                    tasks.push_back(node);
                    workerId = chosenNextWorker(workers, node);
                } else {
                    workerId = workers.size();
                    workers.push_back(node);
                    taskId = chosenNextTask(tasks, node);
                }

                if (workerId >= 0 && taskId >= 0) {
                    addOneMatch(tasks[taskId], workers[workerId]);
                    matched_worker++;
                    cout << "<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility << " " << matched_worker << endl;
                    workers[workerId].begTime=workers[workerId].begTime+sqrt(Length2(workers[workerId].loc, tasks[taskId].loc))/workers[workerId].vel+1;
                    workers[workerId].loc.first=tasks[taskId].loc.first;
                    workers[workerId].loc.second=tasks[taskId].loc.second;
                    if (workers[workerId].begTime < workers[workerId].endTime)
                        workingworker.push_back(workers[workerId]);
                    // task_done.push_back(taskId);
                }
                temptasks.clear();
                tempworkers.clear();
                for (int t = 0; t < tasks.size(); t++) {
                    if (tasks[t].flow <1 &&
                        tasks[t].endTime > node.begTime)//&&tasks[t].endTime>=node.begTime
                        temptasks.push_back(tasks[t]);
                }
                for (int t = 0; t < workers.size(); t++) {
                    if (workers[t].flow < 1 &&
                        workers[t].endTime > node.begTime)//&&workers[t].endTime>=node.begTime
                        tempworkers.push_back(workers[t]);
                }
                tasks.clear();
                workers.clear();
                for (int t = 0; t < temptasks.size(); t++) {
                    tasks.push_back(temptasks[t]);
                }
                for (int t = 0; t < tempworkers.size(); t++) {
                    workers.push_back(tempworkers[t]);
                }
            }
            batch.clear();
		}
	}
	
#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
}

void solve(string fileName) {
	int taskN, workerN, seqN, sumC;
	double Umax;
	ifstream fin(fileName.c_str(), ios::in);

	if (!fin.is_open()) {
		printf("Error openning FILE %s.\n", fileName.c_str());
		exit(1);
	}

	fin >> workerN >> taskN >> Umax >> sumC;
	seqN = taskN + workerN;
	init(taskN, workerN, Umax, sumC);
	Pure_Greedy(fin, seqN);
}

int main(int argc, char* argv[]) {
	cin.tie(0);
	ios::sync_with_stdio(false);

	string edgeFileName;
	program_t begProg, endProg;
    double sumUtility=0,sumUsedTime=0;
    int sumUsedMemory=0, sumWorker=0;

    for(int i=0;i<100;i++) {
        edgeFileName.clear();
        edgeFileName += "./data/synthetic2/workers/1000_3000_1_10_0.5_6_10/data_";
        //edgeFileName += "./data/synthetic1/task/500_2500_1_10_0.5_6_10/data_";
        //edgeFileName += "./data/synthetic1/cw/500_2500_3_10_0.5_6_10/data_";
        //edgeFileName += "./data/scaleData/syn2/2000_10000_2_10_0.5_6_10/data_";
        //edgeFileName += "./data/real/EverySender_cap1/800/data_80";

        edgeFileName += to_string(static_cast<long long>(i));;
        edgeFileName += ".txt";
        /*if (argc > 1) {
            edgeFileName = string(argv[1]);
        }else{
            edgeFileName="./data/synthetic1/worker/100_2500_1_10_0.5_6_10/data_00.txt";
        }*/

        save_time(begProg);
        solve(edgeFileName);
        watchSolutionOnce(getpid(), usedMemory);
        save_time(endProg);

        double usedTime = calc_time(begProg, endProg);
#ifdef WATCH_MEM
        printf("Greedy %.6lf %.6lf %d\n", utility, usedTime, usedMemory/1024);
#else
        printf("Random %d %.6lf %d  %.6lf %d\n", i+1, utility, matched_worker, usedTime,usedMemory);
        sumUtility=sumUtility+utility;
        sumUsedTime=sumUsedTime+usedTime;
        sumUsedMemory=sumUsedMemory+usedMemory;
        sumWorker=sumWorker+matched_worker;
        matched_worker=0;
        utility=0;
        usedTime=0;
        usedMemory=0;
        task_done.clear();
#endif
        fflush(stdout);
    }
    printf("Random  %.6lf %d %.6lf %d\n",  sumUtility/100, sumWorker/100 , sumUsedTime/100, sumUsedMemory/100);
    return 0;
}
