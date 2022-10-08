/**
	\author: Trasier
	\date: 2017.04.02
*/
#include <bits/stdc++.h>
using namespace std;
//#pragma comment(linker,"/STACK:102400000,1024000")
#include "monitor.h"
#include "input.h"


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
    int id;             //
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
double umax;
double utility;
int usedMemory = 0;
int matched_worker=0;
float p=0.05; //moving cost of worker   p=0.1
vector<string> task_done;

void init(int taskN, int workerN, double Umax, int sumc) {
	n = workerN;
	m = taskN;
	umax = Umax;
	utility = 0;
	usedMemory = 0;
	sumC = sumc;
}

void nextSeq(ifstream& fin, node_t& nd) {
	int timeId;
	string stype;

	fin >> nd.id>>nd.begTime >> stype;
	if (stype[0] == 'w') {
		nd.type = worker;
		fin >> nd.loc.first >> nd.loc.second >> nd.rad >> nd.vel >> nd.endTime >> nd.cost.pay;
		nd.endTime += nd.begTime;
        //nd.vel=2;
	} else {
		nd.type = task;
		fin >> nd.loc.first >> nd.loc.second >> nd.endTime >> nd.cost.rate;
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
    double aa=Length(worker.loc, task.loc);
    int a=(aa/worker.vel)+1;

    if (!(worker.begTime+a<task.endTime && worker.begTime+a<worker.endTime && task.begTime<worker.endTime))
        return false;
    return true;
}
inline bool satisfyCost(const node_t& worker, const node_t& task) {
    // 1. condition of deadline
    //if (!(worker.begTime<=task.endTime && task.begTime<=worker.endTime))
    double a=task.cost.pay * worker.cost.rate-(Length(worker.loc, task.loc)/worker.vel)*p;
    if (task.cost.pay * worker.cost.rate-(Length(worker.loc, task.loc)/worker.vel)*p<=0)
        return false;
    return true;
}

inline double calcCost(const node_t& task, const node_t& worker) {
    double cost=(Length(worker.loc, task.loc)/worker.vel)*p;
    return task.cost.pay * worker.cost.rate-cost;
}

bool satisfy(const node_t& worker, const node_t& task) {
	return satisfyCap(worker, task) && satisfyTime(worker, task) && satisfyLoc(worker, task)&& satisfyCost(worker, task);
}

int chosenNextTask(const vector<node_t>& tasks, const node_t& worker, double costBound) {
	int taskN = tasks.size();
	double tmpCost;
	double mxCost = -1e8;
	int ret = -1;

	for (int i=0; i<taskN; ++i) {
        tmpCost = calcCost(tasks[i], worker);
        if (tmpCost > costBound && satisfy(worker, tasks[i])) {
            ret = i;
            break;
        }
    }
	if (ret<0)
    {
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
    }
	return ret;
}

int chosenNextWorker(const vector<node_t>& workers, const node_t& task, double costBound) {
	int workerN = workers.size();
	double tmpCost;
	double mxCost = -1e8;
	int ret = -1;
    for (int i=0; i<workerN; ++i) {
        tmpCost = calcCost(task, workers[i]);
        if (tmpCost>costBound && satisfy(workers[i], task)) {
            ret = i;
            break;
        }
    }
    if (ret<0)
    {
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

void Extend_Greedy_RT(ifstream& fin, int seqN, int k) {
	double costBound = (k==0) ? 0.0 : pow(exp(1.0), k);
	node_t node;
	vector<node_t> tasks,temptasks, workers,tempworkers;
    vector<node_t> workingworker, tworkingworker;

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
            for (int i = 0; i < workingworker.size(); ++i) {
                if (workingworker[i].begTime <= nd.begTime&&workingworker[i].endTime>=nd.begTime) {
                    workingworker[i].flow = 0;
                    batch.push_back(workingworker[i]);
                }
            }
            tworkingworker.clear();
            for (int t = 0; t < workingworker.size(); t++) {
                if (workingworker[t].flow > 0 &&
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
                    workerId = chosenNextWorker(workers, node, costBound);
                    if (workerId >= 0 && taskId >= 0) {
                        addOneMatch(tasks[taskId], workers[workerId]);
                        matched_worker++;
                        cout << "<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility << " "
                             << matched_worker << endl;
                        int das=sqrt(Length2(workers[workerId].loc, tasks[taskId].loc))/workers[workerId].vel;
                        workers[workerId].begTime=workers[workerId].begTime+sqrt(Length2(workers[workerId].loc, tasks[taskId].loc))/workers[workerId].vel+1;
                        workers[workerId].loc.first=tasks[taskId].loc.first;
                        workers[workerId].loc.second=tasks[taskId].loc.second;
                        if (workers[workerId].begTime<workers[workerId].endTime)
                            workingworker.push_back(workers[workerId]);
                    }
                } else {
                    workerId = workers.size();
                    workers.push_back(node);
                    for (int i = 0; i < 1; i++) {
                        taskId = chosenNextTask(tasks, node, costBound);
                        if (workerId >= 0 && taskId >= 0) {
                            addOneMatch(tasks[taskId], workers[workerId]);
                            matched_worker++;
                            cout << "<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility << " "
                                 << matched_worker << endl;
                            workers[workerId].begTime=workers[workerId].begTime+sqrt(Length2(workers[workerId].loc, tasks[taskId].loc))/workers[workerId].vel+1;
                            workers[workerId].loc.first=tasks[taskId].loc.first;
                            workers[workerId].loc.second=tasks[taskId].loc.second;
                            if (workers[workerId].begTime<workers[workerId].endTime)
                                workingworker.push_back(workers[workerId]);
                        }
                    }
                }
                temptasks.clear();
                tempworkers.clear();
                for (int t = 0; t < tasks.size(); t++) {
                    if (tasks[t].flow < 1 &&
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
#ifdef WATCH_MEM
		watchSolutionOnce(getpid(), usedMemory);
#endif
	}
}

double calcUtility(const vector<double>& utilities) {
	double ret;
	const int sz = utilities.size();

	ret = 0.0;
	for (int i=0; i<sz; ++i)
		ret += utilities[i];
	ret /= sz;

	return ret;
}

void solve(string fileName) {
	int taskN, workerN, seqN, sumC;
	double Umax;
	{// get Umax to calculate theta
		ifstream fin(fileName.c_str(), ios::in);

		if (!fin.is_open()) {
			printf("Error openning FILE %s.\n", fileName.c_str());
			exit(1);
		}

		fin >> workerN >> taskN >> Umax >> sumC;
		fin.close();
	}
	
	int theta = ceil(log(Umax + 1.0));
	vector<double> utilities;
	
	for (int i=0; i<theta; ++i) {
		int k = (theta==0) ? 0 : i;
		ifstream fin(fileName.c_str(), ios::in);

		if (!fin.is_open()) {
			printf("Error openning FILE %s.\n", fileName.c_str());
			exit(1);
		}

		fin >> workerN >> taskN >> Umax >> sumC;
		seqN = taskN + workerN;
		init(taskN, workerN, Umax, sumC);
		Extend_Greedy_RT(fin, seqN, k);
		
		utilities.push_back(utility);
		
		printf("k = %d, utility = %.6lf\n", k, utility);
		
		fin.close();
	}
	
	utility = calcUtility(utilities);
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
        edgeFileName += "./data/synthetic2/tasks/500_5000_1_10_0.5_6_10/data_";
        //edgeFileName += "./data/real/EverySender_cap1/cw/1/data_cw";
        //edgeFileName+="./data/scaleData1/syn2/10000_50000_1_10_0.5_6_10_0.5_3/data_" ;
        //edgeFileName += "./data/synthetic6/bra/500_2500_1_10_0.5_6_10_0.5_5/data_";
        //edgeFileName += "./data/scaleData/syn2/2000_10000_2_10_0.5_6_10/data_";
        // edgeFileName += "./data/real/EverySender_cap1/100/data_10";
        edgeFileName += to_string(static_cast<long long>(i));
        edgeFileName += ".txt";
        /*if (argc > 1) {
            edgeFileName = string(argv[1]);
        } else{
            edgeFileName="./data/greedy/data_00.txt";
        }*/

        save_time(begProg);
        solve(edgeFileName);
        save_time(endProg);

        watchSolutionOnce(getpid(), usedMemory);
        double usedTime = calc_time(begProg, endProg);
        int theta = ceil(log(umax + 1.0));
#ifdef WATCH_MEM
       // printf("Ext-GRT %.6lf %.6lf %d\n", utility, usedTime/theta, usedMemory/1024);
#else
        printf("Ext-GRT %d %.6lf %d %.6lf %d\n",i+1, utility, matched_worker/4, usedTime / theta, usedMemory);
#endif
        fflush(stdout);
        sumUtility=sumUtility+utility;
        sumUsedTime=sumUsedTime+usedTime/theta;
        sumUsedMemory=sumUsedMemory+usedMemory;
        sumWorker=sumWorker+matched_worker/3;
        utility=0;
        usedTime=0;
        usedMemory=0;
        matched_worker=0;
        task_done.clear();
    }
    printf("Ext-GRT  %.6lf %d %.6lf %d\n",  sumUtility/100, sumWorker/100, sumUsedTime/100, sumUsedMemory/100);
	return 0;
}
