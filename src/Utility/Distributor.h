#pragma once
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
#include "../data.h"
#include <algorithm>
#include <random>
#include "../Transitions.h"
#include "../HarmonicBayes/HarmonicFunctions.h"
#include "../HarmonicTree/HarmonicNetwork.h"
class ErroredBinomial;


class Distributor
{

	private:
		
		std::vector<int> Register;
		std::vector<std::thread> Threads;
		void Worker(int id);

		std::vector<int> ebStart;
		std::vector<int> ebEnd;
		void EB_Task(int id);
		
		void Assignment_Task(int id);

		void Harmonic_Task(int id);
		void SpoolDown(int id);
		int CheckRegister(int id);
		void SetRegister(int value, int id);
		Data & DataCopy;
		std::vector<std::vector<int>> ChromAssignments;
		double alpha;
		int L;
		int accelerator;
		double nu;
		double gamma;
	public:
		std::vector<int> MainChromAssigment;
		ErroredBinomial * EB;
		std::vector<HarmonicNetwork> & Ns;
		Transitions * Output;
		int WorkerCount;
		void Signal(int value);
		void Gather();
		void UpdateEB(ErroredBinomial * eb);
		void UpdateAssigner(Transitions * Output, double alpha, int L, int accelerator, double nu);
		
		void UpdateParameters(double nu, double gamma);
		
		Distributor(int nWorkers, Data & d,std::vector<HarmonicNetwork> & Ns);

		~Distributor();
};