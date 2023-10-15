#!/usr/bin/env python
# encoding=utf-8
import numpy as np
from random import randint, sample, shuffle
import sys
import os
import math
import random
import xlwt

#dependy = {

#dependy = []
#dependy_des = []
class constForDataSet:
	caseN = 1
	locRng = [0, 100]     #syn set [0-500]
	#TList = [10000, 20000, 30000, 40000 ,50000]
	#WList = [2000, 4000, 6000, 8000, 10000]
	TList = [1000, 2000, 3000, 4000, 5000]
	WList = [100, 200, 300, 500, 500]
	cwList = [1, 2, 3, 4, 5]     #速度
	rwList = [2, 5, 10, 15, 20]  #可达半径
	dwList = [0.1, 0.3, 0.5, 0.7, 0.9] #完成质量
	ddlList = [2, 4, 6, 8, 10]   #工人截止日期
	ptList = [2, 5, 10, 20, 50]
	lambdaList = [2, 5, 10, 20, 50]
	muList = [2, 5, 10, 20, 50]
	defaultValues = {
		"T": TList[len(TList)/2],
		"W": WList[len(WList)/2],
		"cw": cwList[len(cwList)/2],
		"rw": rwList[len(rwList)/2],
		"dw": dwList[len(dwList)/2],
		"ddl": ddlList[len(ddlList)/2],
		"pt": ptList[len(ptList)/2],
		"lambda": lambdaList[len(lambdaList)/2],
		"muList": muList[len(muList)/2],
	}
	umax = 10
	umin = 0
	tmax = 1000  #syn4 and syn6 set is 600  scla set 3000
	per = 50

class CFDS(constForDataSet):
	pass



class baseGenerator:

	def gen(self, n):
		pass

class randomGenerator(baseGenerator):

	def __init__(self, mx):
		self.mx = mx

	def setMx(self, mx):
		self.mx = mx

	def gen(self, n):
		ret = [0] * n
		for i in xrange(n):
			x = randint(1, self.mx)
			ret[i] = x
		return ret

class normalGenerator(baseGenerator):

	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma

	def gen(self, n, lb = None, rb = None):
		# print lb, rb
		ret = np.random.normal(self.mu, self.sigma, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		# print ret
		return ret

	def setMu(self, mu):
		self.mu = mu

	def setSigma(self, sigma):
		self.sigma = sigma


class uniformGenerator(baseGenerator):

	def __init__(self, low, high):
		self.low = low
		self.high = high

	def gen(self, n, lb = None, rb = None):
		ret = np.random.uniform(self.low, self.high, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setLow(self, low):
		self.low = low

	def setHigh(self, high):
		self.high = high

class expGenerator(baseGenerator):

	def __init__(self, mu):
		self.mu = mu

	def gen(self, n, lb = None, rb = None):
		ret = np.random.exponential(self.mu, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setMu(self, mu):
		self.mu = mu


def genLoc(n, low=0, high=100):
	ret = []
	st = set()
	for i in xrange(n):
		while True:
			x = 1.0*randint(low*10**6, high*10**6) / 10**6
			y = 1.0*randint(low*10**6, high*10**6) / 10**6
			t = (x, y)
			if t not in st:
				break
		ret.append(t)
		st.add(t)
	return ret


def genTmpFileName(W, T, cw, rw, dw, ddl, pt,):
	return "%s_%s_%s_%s_%s_%s_%s" % (W, T, cw, rw, dw, ddl, pt)


def genData(desFileName, _oids_w, _oids_t, wids, wlocs, tids, tlocs, cw, rw, dws, ddl, pts):
	# print type(_oids_w), _oids_w
	# print type(_oids_t), _oids_t
	with open(desFileName, "w") as fout:
		W, T, umax, sumc = len(wids), len(tids), CFDS.umax, len(tids)+len(wids)
		fout.write("%d %d %d %d\n" %(W, T, umax, sumc))
		i, j = 0, 0
		tids.sort();
		time=wids[W-1]
		generator = randomGenerator(time)
		ntids = sorted(generator.gen(T))
		while i<W and j<T:
			if wids[i] <= tids[j]:
				k = _oids_w[i]
				fout.write("%d %d w %.6lf %.6lf %.6lf %d %d %.6lf\n" % (i, wids[i], wlocs[k][0], wlocs[k][1], rw, cw, ddl, dws[k]))
				i += 1
			else:
				k = _oids_t[j]
				fout.write("%d %d t %.6lf %.6lf %d %.6lf\n" % (j, tids[j], tlocs[k][0], tlocs[k][1], ddl, pts[k]))
				j += 1
		while i < W:
			k = _oids_w[i]
			fout.write("%d %d w %.6lf %.6lf %.6lf %d %d %.6lf\n" % (i, wids[i], wlocs[k][0], wlocs[k][1], rw, cw, ddl, dws[k]))
			i += 1
		while j < T:
			k = _oids_t[j]
			fout.write("%d %d t %.6lf %.6lf %d %.6lf\n" % (j, tids[j], tlocs[k][0], tlocs[k][1], ddl, pts[k]))
			j += 1


def genOids(W, T):
	oids_w = range(0, W)
	oids_t = range(0, T)
	return oids_w, oids_t


def nextRandOrder(oids_w, oids_t):
	_oids_w, _oids_t = oids_w, oids_t
	shuffle(_oids_w)
	shuffle(_oids_t)
	return _oids_w, _oids_t


def nextRestrictOrder(oids_w, oids_t):
	_oids_w, _oids_t = oids_w, []
	m, n = len(oids_t), len(oids_w)
	curPer = int( math.ceil(len(oids_t)*1.0 / len(oids_w)) )
	shuffle(_oids_w)
	for i in xrange(n):
		wid = _oids_w[i]
		fr, to = wid*curPer, (wid+1)*curPer
		if fr >= len(oids_t):
			fr = len(oids_t)
		if to >= len(oids_t):
			to = len(oids_t)
		tmpList = oids_t[fr:to]
		shuffle(tmpList)
		_oids_t += tmpList
	assert len(_oids_t) == len(oids_t)
	return _oids_w, _oids_t


def nextOrder(oids_w, oids_t):
	# _oids_w, _oids_t = nextRandOrder(oids_w, oids_t)
	_oids_w, _oids_t = nextRestrictOrder(oids_w, oids_t)
	return _oids_w, _oids_t


def genWorkers(n):

	generator = randomGenerator(n/50)
	generator = randomGenerator(CFDS.tmax)
	#tids = []
	#for i in xrange(n):
		#tids.append(i+1)
	# #for i in xrange(n/2):
	# 	#tids.append(2*i+1)
	# #for i in xrange(n/2):
	# 	#tids.append(n+i+1)
	tids = sorted(generator.gen(n))
	locs = genLoc(n, CFDS.locRng[0], CFDS.locRng[1])
	return tids, locs


def genDws(W, dw, sigma = 0.05):
	generator = normalGenerator(dw, sigma)
	#generator = randomGenerator(dw)
	dws = generator.gen(W, 0.01, 1.0)
	return dws

	
def genPts(T, m, n, minPts, maxPts):
	ig, jg = 0, 0
	pts = []
	for i in xrange(n):
		for j in xrange(CFDS.per):
			if j>=int(CFDS.per*0.5):
				pt = maxPts[ig]
				ig = (ig+1) % len(maxPts)
			else:
				pt = minPts[jg]
				jg = (jg+1) % len(minPts)
			pts.append(pt)
	return pts[:T]


def genPts_Normal(T, m, n, pt, sigma=3.75):
	minGenerator = normalGenerator(pt*0.1, sigma/40.0)
	maxGenerator = normalGenerator(pt, sigma/10.0)
	minPts = minGenerator.gen(m*CFDS.per, 10**-4, CFDS.umax)
	maxPts = maxGenerator.gen(m*CFDS.per, pt*0.1, CFDS.umax)
	return genPts(T, m, n, minPts, maxPts)

def genPts_Uniform(T, m, n, pt):
	minGenerator = uniformGenerator(0, pt*0.5)
	maxGenerator = uniformGenerator(0, pt*1.0)
	minPts = minGenerator.gen(m*CFDS.per, 10**-4, CFDS.umax)
	maxPts = maxGenerator.gen(m*CFDS.per, pt*0.1, CFDS.umax)
	return genPts(T, m, n, minPts, maxPts)
	
def genPts_Exponential(T, m, n, pt):
	minGenerator = expGenerator(pt*0.1)
	maxGenerator = expGenerator(pt)
	minPts = minGenerator.gen(m*CFDS.per, 10**-4, CFDS.umax)
	maxPts = maxGenerator.gen(m*CFDS.per, pt*0.1, CFDS.umax)
	return genPts(T, m, n, minPts, maxPts)
	

'''
def sampleTasks(m, n, tids, tlocs, pts):
	per = CFDS.per
	curPer = int( math.ceil(m*1.0/n) )
	_tids = []
	_tlocs = []
	_pts = []
	for i in xrange(n):
		for j in xrange(curPer):
			if j>=int(curPer*0.5):
				k = i*per + int(per*0.7) + j - int(curPer*0.5)
			else:
				k = i*per + j
			# print k, len(tids), len(tlocs), m, n, curPer
			tid, tloc, pt = tids[k], tlocs[k], pts[k]
			_tids.append(tid)
			_tlocs.append(tloc)
			_pts.append(pt)
	return _tids[:m], _tlocs[:m], _pts[:m]
'''

'''
#工人和任务的位置随机产生
def genTasks(m, n, wids, wlocs, rw):
	per = CFDS.per
	tids = []
	tlocs = []
	locs = genLoc(m, CFDS.locRng[0], CFDS.locRng[1])
	k=0
	for i in xrange(n):
		for j in xrange(per):
			tid = wids[i];
			tloc_x = locs[i][0]
			tloc_y = locs[i][1]
			tloc = (tloc_x, tloc_y)
			tids.append(tid)
			tlocs.append(tloc)
	return tids, tlocs

def sampleTasks(m, n, tids, tlocs, pts):
	per = CFDS.per
	_tids = []
	_tlocs = []
	_pts = []
	temp = 0
	temp1=0.0  #m:task n:worker
	curPer =(m*1.0)/n
	flagcurPer=curPer
	for i in xrange(n):
		if flagcurPer >= 1:
			temp=int(flagcurPer//1)
			flagcurPer=flagcurPer-temp
			for j in xrange(temp):
				ranv=randint(1, 2)
				#l=randint(1, len(tlocs)-1)
				if ranv==1:
					k = i*per + randint(1, 35)
					# print k, len(tids), len(tlocs), m, n, curPer
					tid, tloc, pt = tids[k], tlocs[k], pts[k]+random.random()
					_tids.append(tid)
					_tlocs.append(tloc)
					_pts.append(pt)
				else:
					k = i*per + randint(36, 49)
					# print k, len(tids), len(tlocs), m, n, curPer
					tid, tloc, pt = tids[k], tlocs[k], pts[k]-random.random()
					_tids.append(tid)
					_tlocs.append(tloc)
					_pts.append(pt)
			flagcurPer = flagcurPer + curPer
		else:
			flagcurPer = flagcurPer + curPer
	if flagcurPer>0:
		k = i * per + 1
		tid, tloc, pt = tids[k], tlocs[k], pts[k]
		_tids.append(tid)
		_tlocs.append(tloc)
		_pts.append(pt)
	return _tids[:m], _tlocs[:m], _pts[:m]

'''
def genTasks(m, n, wids, wlocs, rw):
	per = CFDS.per
	tids = []
	tlocs = []
	for i in xrange(n):
		for j in xrange(per):
			tid = wids[i];
			if j>=int(per*0.5):
				#l=len(wids)iii
				#k = randint(wids[i], wids[l-1])
				tloc_x = wlocs[i][0] + 1.0*randint(-0.2*rw*10**6, 0.2*rw*10**6) / 10**6
				tloc_y = wlocs[i][1] + 1.0*randint(-0.2*rw*10**6, 0.2*rw*10**6) / 10**6
			else:
				tloc_x = wlocs[i][0] + 1.0*randint(-0.2*rw*10**6, 0.2*rw*10**6) / 10**6
				tloc_y = wlocs[i][1] + 1.0*randint(-0.2*rw*10**6, 0.2*rw*10**6) / 10**6
			tloc = (tloc_x, tloc_y)
			tids.append(tid)
			tlocs.append(tloc)
	return tids, tlocs
	
def sampleTasks(m, n, tids, tlocs, pts):
	per = CFDS.per
	_tids = []
	_tlocs = []
	_pts = []
	temp = 0
	temp1=0.0  #m:task n:worker
	curPer =(m*1.0)/n
	flagcurPer=curPer
	for i in xrange(n):
		if flagcurPer >= 1:
			temp=int(flagcurPer//1)
			flagcurPer=flagcurPer-temp
			for j in xrange(temp):
				ranv=randint(1, 2)
				if ranv==1:
					k = i*per + randint(1, 35)
					# print k, len(tids), len(tlocs), m, n, curPer
					tid, tloc, pt = tids[k], tlocs[k], pts[k]+random.random()
					_tids.append(tid)
					_tlocs.append(tloc)
					_pts.append(pt)
				else:
					k = i*per + randint(36, 49)
					# print k, len(tids), len(tlocs), m, n, curPer
					tid, tloc, pt = tids[k], tlocs[k], pts[k]-random.random()
					_tids.append(tid)
					_tlocs.append(tloc)
					_pts.append(pt)
			flagcurPer = flagcurPer + curPer
		else:
			flagcurPer = flagcurPer + curPer
	if flagcurPer>0:
		k = i * per + 1
		tid, tloc, pt = tids[k], tlocs[k], pts[k]
		_tids.append(tid)
		_tlocs.append(tloc)
		_pts.append(pt)
	return _tids[:m], _tlocs[:m], _pts[:m]



def _genPts(T, W, pt, sigma=3.75):
	per = CFDS.per
	curPer = int( math.ceil(T*1.0/W) )
	minGenerator = normalGenerator(pt/10.0, sigma/40.0)
	maxGenerator = normalGenerator(pt, sigma/10.0)
	minPts = minGenerator.gen(T*per, 10**-4, CFDS.umax)
	maxPts = maxGenerator.gen(T*per, pt/10.0, CFDS.umax)
	ig, jg = 0, 0
	pts = []
	for i in xrange(W):
		for j in xrange(curPer):
			if j>=int(curPer*0.7):
				pt = maxPts[ig]
				ig = (ig+1) % len(maxPts)
			else:
				pt = minPts[jg]
				jg = (jg+1) % len(minPts)
			pts.append(pt)
	return pts[:T]


def genDataSet(desFilePath, orderN = 10):
	# vary number of workers
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]
	for W in CFDS.WList:
		wids, wlocs = genWorkers(W+T)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)


	# vary number of tasks
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]
	for T in CFDS.TList:
		wids, wlocs = genWorkers(W+T/10)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		#_d, _des = genTaskDependy(pro, bra)

		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of cw
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]

	for cw in CFDS.cwList:
		wids, wlocs = genWorkers(W)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t

		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of rw
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]

	for rw in CFDS.rwList:
		wids, wlocs = genWorkers(W)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of quel
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]

	for dw in CFDS.dwList:
		wids, wlocs = genWorkers(W+T)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		_dws = genDws(W, dw)
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of ddl
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]
	for ddl in CFDS.ddlList:
		wids, wlocs = genWorkers(W)
		tids, tlocs = genTasks(T, W, wids, wlocs, max(CFDS.rwList))
		pts = genPts_Uniform(W * CFDS.per, T, W, CFDS.defaultValues["pt"])
		dws = genDws(W, CFDS.defaultValues["dw"])
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of pt under normal distribution
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]

	for pt in CFDS.ptList:
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		_pts = _genPts(T, W, pt)
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of pt under exponetial distribution
	#pts = genPts_Exponential(wmax*CFDS.per, tmax, wmax, CFDS.defaultValues["pt"])
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]
	for pt in CFDS.ptList:
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt+100)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		_pts = _genPts(T, W, pt)
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)

	# vary value of pt under uniform distribution
	#pts = genPts_Uniform(wmax*CFDS.per, tmax, wmax, CFDS.defaultValues["pt"])
	W = CFDS.defaultValues["W"]
	T = CFDS.defaultValues["T"]
	cw = CFDS.defaultValues["cw"]
	rw = CFDS.defaultValues["rw"]
	dw = CFDS.defaultValues["dw"]
	ddl = CFDS.defaultValues["ddl"]
	pt = CFDS.defaultValues["pt"]
	for pt in CFDS.ptList:
		tmpFilePath = genTmpFileName(W, T, cw, rw, dw, ddl, pt+1000)
		tmpFilePath = os.path.join(desFilePath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		_wids, _wlocs, _dws = wids[:W], wlocs[:W], dws[:W]
		_tids, _tlocs, _pts = sampleTasks(T, W, tids, tlocs, pts)
		__oids_w, __oids_t = genOids(W, T)
		_oids_w, _oids_t = __oids_w, __oids_t
		_pts = _genPts(T, W, pt)
		#_d, _des = genTaskDependy(pro, bra)
		for oid in xrange(orderN):
			desFileName = "data_%02d.txt" % (oid)
			desFileName = os.path.join(tmpFilePath, desFileName)
			if os.path.exists(desFileName):
				continue
			genData(desFileName, _oids_w, _oids_t, _wids, _wlocs, _tids, _tlocs, cw, rw, _dws, ddl, _pts)
			if oid < orderN-1:
				_oids_w, _oids_t = nextOrder(__oids_w, __oids_t)
					
				
def exp1():
	desFilePath = "./synthetic4"
	if not os.path.exists(desFilePath):
		os.mkdir(desFilePath)
	#genTaskDependy()
	genDataSet(desFilePath, CFDS.caseN)
	print "DONE."
	
if __name__ == "__main__":
	exp1()
