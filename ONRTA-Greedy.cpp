/*
	\author: Trasier
	\date: 2017.04.02
*/
#include <bits/stdc++.h>
using namespace std;
//#pragma comment(linker,"/STACK:102400000,1024000")
#include "input.h"
#include "monitor.h"

vector<string> task_done;

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
    int cap;			// velocity
    int flow;			// flow
    double rad;			// radius
    W_un cost;			// cost
    int id;
    int begTime, endTime;	// time interval

    void print() {
        if (type == worker)
            printf("w: loc = (%.2lf, %.2lf), rad = %.2lf, cap = %d, time = (%d, %d), ratio = %.2lf\n",
                   loc.first, loc.second, rad, cap, begTime, endTime, cost.rate);
        else
            printf("t: loc = (%.2lf, %.2lf), time = (%d, %d), pay = %.2lf\n",
                   loc.first, loc.second, begTime, endTime, cost.pay);
    }
};
int matched_worker=0;

bool satisfyLoc(const node_t& worker, const node_t& task);
bool satisfyTime(const node_t& worker, const node_t& task);



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
    if (worker.cap<=worker.flow || task.cap<=task.flow)
        return false;
    return true;
}

inline bool satisfyTime(const node_t& worker, const node_t& task) {
    // 1. condition of deadline
    //if (!(worker.begTime<=task.endTime && task.begTime<=worker.endTime))
    if (!(worker.begTime<task.endTime && task.begTime<worker.endTime))
        return false;
    return true;
}

inline double calcCost(const node_t& task, const node_t& worker) {
    return task.cost.pay * worker.cost.rate;
}

bool satisfy(const node_t& worker, const node_t& task) {
    return satisfyCap(worker, task) && satisfyTime(worker, task) && satisfyLoc(worker, task);
}

const double eps = 1e-6;
inline double dcmp(double a) {
	if (fabs(a) < eps)
		return 0;
	return a>0 ? 1:-1;
}

const double INF = 1e18;
struct vertex_t {
	int u, v;
	double w;
	
	vertex_t(int u=0, int v=0, double w=0):
		u(u), v(v), w(w) {}
		
	bool operator<(const vertex_t& o) const {
		if (w != o.w)
			return w > o.w;
		if (u != o.u)
			return u < o.u;
		return v < o.v;
	}
};
set<vertex_t> est;
struct Greedy_t {
	vector<int> yx, xy;

	Greedy_t() {
		init(0, 0);
	}
	
	void init(int Wsz, int Tsz) {
		clear();
		yx.resize(Tsz, -1);
		xy.resize(Wsz, -1);
	}
	
	void clear() {
		yx.clear();
		xy.clear();
	}

    void build(const vector<int>& T_delta, const vector<int>& W_delta, const vector<node_t>& tasks,
               const vector<node_t>& workers, const node_t& node) {
        const int Tsz = T_delta.size();
        const int Wsz = W_delta.size();
        double cost;
        est.clear();
        init(Wsz, Tsz);
        if (node.type == worker) {
            //const int i = Wsz - 1;
            for(int i=0;i<Wsz;i++) {
                const int workerId = W_delta[i];
                for (int j = 0; j < Tsz; ++j) {
                    const int taskId = T_delta[j];
                    if (satisfyLoc(workers[workerId], tasks[taskId]) && satisfyTime(workers[workerId], tasks[taskId])) {
                        cost = calcCost(tasks[taskId], workers[workerId]);
                        est.insert(vertex_t(i, j, cost));
                    }
                }
            }
        } else {
            for (int j = 0; j < Tsz; j++) {
                //if(Tsz!=0) {
                const int taskId = T_delta[j];
                for (int i = 0; i < Wsz; ++i) {
                    const int workerId = W_delta[i];
                    if (satisfyLoc(workers[workerId], tasks[taskId]) && satisfyTime(workers[workerId], tasks[taskId])) {
                        cost = calcCost(tasks[taskId], workers[workerId]);
                        est.insert(vertex_t(i, j, cost));
                    }
                }
            }
        }
    }
	
	void GreedyMatch(const vector<int>& T_delta, const vector<int>& W_delta, const vector<node_t>& tasks, 
				const vector<node_t>& workers, const node_t& node) {
		const int Tsz = T_delta.size();
		const int Wsz = W_delta.size();
		vertex_t ver;
		int taskId, workerId;
		
		if (node.type == worker) {
			for (set<vertex_t>::iterator iter=est.begin(); iter!=est.end(); ++iter) {
				ver = *iter;

				if (xy[ver.u]>=0 || yx[ver.v]>=0) continue;

				workerId = W_delta[ver.u];
				taskId = T_delta[ver.v];
				const node_t& workerNode = workers[workerId];
				const node_t& taskNode = tasks[taskId]; 

				xy[ver.u] = ver.v;
				yx[ver.v] = ver.u;
			}
		} else {
			for (set<vertex_t>::iterator iter=est.begin(); iter!=est.end(); ++iter) {
				ver = *iter;

				if (xy[ver.u]>=0 || yx[ver.v]>=0) continue;

				workerId = W_delta[ver.u];
				taskId = T_delta[ver.v];
				const node_t& workerNode = workers[workerId];
				const node_t& taskNode = tasks[taskId];

				xy[ver.u] = ver.v;
				yx[ver.v] = ver.u;
			}
		}
	}
	
	void match(const vector<int>& T_delta, const vector<int>& W_delta, const vector<node_t>& tasks, 
				const vector<node_t>& workers, const node_t& node) {
		GreedyMatch(T_delta, W_delta, tasks, workers, node);
	}
};

typedef long long LL;
int n, m, sumC;
double umax;
double utility;
int usedMemory;
Greedy_t greedy;

void init(int taskN, int workerN, double Umax, int SumC) {
	n = workerN;
	m = taskN;
	umax = Umax;
	sumC = SumC;
	utility = 0;
	usedMemory = 0;
	est.clear();
}

void nextSeq(ifstream& fin, node_t& nd) {
	int timeId;
	string stype;

	fin >> nd.id>>nd.begTime >> stype;
	if (stype[0] == 'w') {
		nd.type = worker;
		fin >> nd.loc.first >> nd.loc.second >> nd.rad >> nd.cap >> nd.endTime >> nd.cost.rate;
		nd.endTime += nd.begTime;
	} else {
		nd.type = task;
		fin >> nd.loc.first >> nd.loc.second >> nd.endTime >> nd.cost.pay;
		nd.endTime += nd.begTime;
		nd.cap = 1;
	}

	nd.flow = 0;
}

int chosenNextTask(const vector<node_t>& tasks, node_t& worker) {
	int taskN = tasks.size();
	double tmpCost;
	double mxCost = 0;
	int ret = -1;
	
	for (int i=0; i<taskN; ++i) {
		tmpCost = calcCost(tasks[i], worker);
		if (satisfy(worker, tasks[i]) && tmpCost>mxCost) {
			mxCost = tmpCost;
			ret = i;
		}
	}
	
	return ret;
}

int chosenNextWorker(const vector<node_t>& workers, node_t& task) {
	int workerN = workers.size();
	double tmpCost;
	double mxCost = 0;
	int ret = -1;
	
	for (int i=0; i<workerN; ++i) {
		tmpCost = calcCost(task, workers[i]);
		if (satisfy(workers[i], task) && tmpCost>mxCost) {
			mxCost = tmpCost;
			ret = i;
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

void TGOA_Greedy(ifstream& fin, int seqN) {
	int k = sumC / 2,count=0;
	vector<int> W_delta, T_delta;
	node_t node;
	vector<node_t> tasks, workers,temptasks,tempworkers;
	int taskId, workerId;
	bool isSecondHalf = false;
	vector <node_t> Sequence;
    vector <node_t> batch;

	while (seqN--) {
		nextSeq(fin, node);
		Sequence.push_back(node);
	}
	for (int id = 0; id < Sequence.size(); ++ id) {
        batch.push_back(Sequence[id]);
        node_t nd = Sequence[id];
        if (id == Sequence.size() - 1 || nd.begTime != Sequence[id + 1].begTime) {
            int Tdelsz = T_delta.size();
            int Wdelsz = W_delta.size();
            if (!isSecondHalf) {
                for (int ii = 0; ii < batch.size(); ++ ii) {
                    node = batch[ii];
                    int cap = node.cap;
                    node.cap = 1;
                    while (cap--) {
                        count++;
                        if (node.type == task) { // node is task
                            taskId = tasks.size();
                            tasks.push_back(node);
                            workerId = chosenNextWorker(workers, node);
                            if (workerId >= 0 && taskId >= 0) {
                                addOneMatch(tasks[taskId], workers[workerId]);
                                matched_worker++;
                                cout << "<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility << " "
                                     << matched_worker << endl;
                            }
                        } else {
                            workerId = workers.size();
                            workers.push_back(node);
                            taskId = chosenNextTask(tasks, node);
                            if (workerId >= 0 && taskId >= 0) {
                                addOneMatch(tasks[taskId], workers[workerId]);
                                matched_worker++;
                                cout << "<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility << " "
                                     << matched_worker << endl;
                            }
                        }
                    }
                }
            } else {
                for (int ii = 0; ii < batch.size(); ++ii) {
                    node = batch[ii];
                    int cap = node.cap;
                    node.cap = 1;
                    while (cap--) {
                        count++;
                        workerId = taskId = -1;
                        if (node.type == task) { // node is task
                            taskId = tasks.size();
                            tasks.push_back(node);
                            T_delta.push_back(taskId);
                        } else {
                            workerId = workers.size();
                            workers.push_back(node);
                            for (int i = 0; i < 1; ++i) {
                                W_delta.push_back(workerId);
                            }
                        }
                    }
                }
                        greedy.build(T_delta, W_delta, tasks, workers, node);
                        greedy.match(T_delta, W_delta, tasks, workers, node);
                        const int Tsz = T_delta.size();
                        const int Wsz = W_delta.size();
                        for (int i = Wdelsz; i < Wsz; ++i) {
                            workerId = W_delta[i];
                            if (greedy.xy[i] >= 0) {
                                taskId = T_delta[greedy.xy[i]];
                                if (satisfy(workers[workerId], tasks[taskId])) {
                                    /* valid, do nothing*/
                                } else {
                                    taskId = -1;
                                }
                            }
                            if (workerId >= 0 && taskId >= 0) {
                                if (satisfy(workers[workerId], tasks[taskId])) {
                                    addOneMatch(tasks[taskId], workers[workerId]);
                                    matched_worker++;
                                    cout << "2<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility
                                         << " "
                                         << matched_worker << endl;
                                }
                            }
                        }
                        for (int i = Tdelsz; i < Tsz; ++i) {
                            taskId = T_delta[i];
                            if (greedy.yx[i] >= 0) {
                                workerId = W_delta[greedy.yx[i]];
                                if (satisfy(workers[workerId], tasks[taskId])) { ;
                                } else {
                                    workerId = -1;
                                }
                            }
                            if (workerId >= 0 && taskId >= 0) {
                                if (satisfy(workers[workerId], tasks[taskId])) {
                                    addOneMatch(tasks[taskId], workers[workerId]);
                                    matched_worker++;
                                    cout << "2<" << tasks[taskId].id << "," << workers[workerId].id << ">" << utility
                                         << " "
                                         << matched_worker << endl;
                                }
                            }
                        }
//                    }
//                }
            }
            if (!isSecondHalf && count + batch.size() >= k) {
                isSecondHalf = true;
            }
            batch.clear();
            temptasks.clear();
            tempworkers.clear();
            for(int t=0;t<tasks.size();t++){
                if(tasks[t].flow<1 &&
                        tasks[t].endTime>node.begTime)//tasks[t].flow<1 && tasks[t].endTime>node.begTime
                    temptasks.push_back(tasks[t]);
            }
            for(int t=0;t<workers.size();t++){
                if(workers[t].flow<1 &&
                        workers[t].endTime>node.begTime)//workers[t].flow<1 && workers[t].endTime>node.begTime
                    tempworkers.push_back(workers[t]);
            }
            tasks.clear();
            workers.clear();
            T_delta.clear();
            W_delta.clear();
            for(int t=0;t<temptasks.size();t++){
                int taskId=tasks.size();
                tasks.push_back(temptasks[t]);
                for (int i=0; i<1; ++i) {
                    T_delta.push_back(taskId);
                }
            }
            for(int t=0;t<tempworkers.size();t++){
                int workerId=workers.size();
                workers.push_back(tempworkers[t]);
                for (int i=0; i<1; ++i) {
                    W_delta.push_back(workerId);
                }
            }
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
	TGOA_Greedy(fin, seqN);
}
#define VMRSS_LINE 17
#define VMSIZE_LINE 13
#define PROCESS_ITEM 14
unsigned int get_proc_mem(unsigned int pid){

    char file_name[64]={0};
    FILE *fd;
    char line_buff[512]={0};
    sprintf(file_name,"/proc/%d/status",pid);

    fd =fopen(file_name,"r");
    if(nullptr == fd){
        return 0;
    }

    char name[64];
    int vmrss;
    for (int i=0; i<VMRSS_LINE-1;i++){
        fgets(line_buff,sizeof(line_buff),fd);
    }

    fgets(line_buff,sizeof(line_buff),fd);
    sscanf(line_buff,"%s %d",name,&vmrss);
    fclose(fd);

    return vmrss;
}
int get_pid(const char* process_name, const char* user = nullptr)
{
    if(user == nullptr){
        user = getlogin();
    }

    char cmd[512];
    if (user){
        sprintf(cmd, "pgrep %s -u %s", process_name, user);
    }

    FILE *pstr = popen(cmd,"r");

    if(pstr == nullptr){
        return 0;
    }

    char buff[512];
    ::memset(buff, 0, sizeof(buff));
    if(NULL == fgets(buff, 512, pstr)){
        return 0;
    }

    return atoi(buff);
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
        edgeFileName += "./data/synthetic2/tasks/500_5000_1_6_0.5_10_10/data_";
        edgeFileName += to_string(static_cast<long long>(i));;
        edgeFileName += ".txt";

        save_time(begProg);
        solve(edgeFileName);
        save_time(endProg);
        watchSolutionOnce(getpid(), usedMemory);
        // int SysMem=curMem-freeMen;
        double usedTime = calc_time(begProg, endProg);

#ifdef WATCH_MEM
        printf("TGOA-Greedy %.6lf %.6lf %d\n", utility, usedTime, usedMemory/1024);
#else
        printf("ONRTA-Greedy %d %.6lf %d %.6lf %d\n", i+1, utility, matched_worker, usedTime, usedMemory);

#endif
        sumUtility=sumUtility+utility;
        sumUsedTime=sumUsedTime+usedTime;
        sumUsedMemory=sumUsedMemory+usedMemory;
        sumWorker=sumWorker+matched_worker;
        matched_worker=0;
        utility=0;
        usedTime=0;
        usedMemory=0;

        fflush(stdout);
    }
    printf("ONRTA-Greedy  %.6lf  %.6lf %d\n",  sumUtility/100, sumUsedTime/100, sumUsedMemory/100);
	return 0;
}
