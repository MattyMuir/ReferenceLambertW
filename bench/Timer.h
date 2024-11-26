#pragma once
#include <chrono>

// Printing logic
#ifdef TIMER_NPRINT
#include <iostream>
#include <format>
#define PRINT(fmt, ...) std::cout << std::format(fmt, __VA_ARGS__);
#else
#include <print>
#define PRINT(fmt, ...) std::print(fmt, __VA_ARGS__);
#endif

#define TIMER(name) Timer name
#define STOP_LOG(name) { name.Stop(); PRINT("{} took: ", #name); name.Log(); }

class Timer
{
protected:
	using Clock = std::chrono::steady_clock;
	using Time = Clock::time_point;
	using Duration = decltype(Time{} - Time{});

	static consteval int64_t CountPerSec()
	{
		using namespace std::chrono;
		return duration_cast<Duration>(1s).count();
	}

public:
	inline Timer()
		: start(Clock::now()) {}

	inline void Stop(bool log = false)
	{
		Time end = Clock::now();
		elapsed = end - start;
		if (log) Log();
	}

	inline void Log() const
	{
		double seconds = (double)elapsed.count() / CountPerSec();
		PRINT("{}ms\n", seconds * 1000.0);
	}

	inline double GetSeconds() const
	{
		return (double)elapsed.count() / CountPerSec();
	}

protected:
	Time start;
	Duration elapsed;
};