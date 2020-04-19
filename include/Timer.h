#ifndef INCLUDE_TIMER_H_
#define INCLUDE_TIMER_H_

#include <chrono>
#include <iostream>

class Timer {
public:
	Timer() {
		m_startPoint = std::chrono::high_resolution_clock::now();
	}

	~Timer() {
		Stop();
	}

	void Stop() {
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startPoint).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();

		auto duration = end - start;

		std::cout << std::fixed << duration * 0.001 << " (ms)\n";
	}
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> m_startPoint;
};

#endif /* INCLUDE_TIMER_H_ */
