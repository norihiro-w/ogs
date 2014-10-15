/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MPIPrettyUnitTestResultPrinter.h"

#include <cstdarg>
#include <string>

namespace testing {

enum GTestColor
{
	COLOR_DEFAULT, COLOR_RED, COLOR_GREEN, COLOR_YELLOW
};

#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE

// Returns the character attribute for the given color.
WORD GetColorAttribute(GTestColor color)
{
	switch (color)
	{
		case COLOR_RED: return FOREGROUND_RED;
		case COLOR_GREEN: return FOREGROUND_GREEN;
		case COLOR_YELLOW: return FOREGROUND_RED | FOREGROUND_GREEN;
		default: return 0;
	}
}

#else

// Returns the ANSI color code for the given color.  COLOR_DEFAULT is
// an invalid input.
const char* GetAnsiColorCode(GTestColor color)
{
	switch (color)
	{
	case COLOR_RED:
		return "1";
	case COLOR_GREEN:
		return "2";
	case COLOR_YELLOW:
		return "3";
	default:
		return NULL;
	};
}

#endif  // GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE

// Helpers for printing colored strings to stdout. Note that on Windows, we
// cannot simply emit special characters and have the terminal change colors.
// This routine must actually emit the characters rather than return a string
// that would be colored when printed, as can be done on Linux.
void ColoredPrintf(GTestColor color, const char* fmt, ...)
{
	va_list args;
	va_start(args, fmt);

#if GTEST_OS_WINDOWS_MOBILE || GTEST_OS_SYMBIAN || GTEST_OS_ZOS || GTEST_OS_IOS
	const bool use_color = false;
#else
	static const bool in_color_mode = true;
//      = ShouldUseColor(posix::IsATTY(posix::FileNo(stdout)) != 0);
	const bool use_color = in_color_mode && (color != COLOR_DEFAULT);
#endif  // GTEST_OS_WINDOWS_MOBILE || GTEST_OS_SYMBIAN || GTEST_OS_ZOS
	// The '!= 0' comparison is necessary to satisfy MSVC 7.1.

	if (!use_color)
	{
		vprintf(fmt, args);
		va_end(args);
		return;
	}

#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
	const HANDLE stdout_handle = GetStdHandle(STD_OUTPUT_HANDLE);

	// Gets the current text color.
	CONSOLE_SCREEN_BUFFER_INFO buffer_info;
	GetConsoleScreenBufferInfo(stdout_handle, &buffer_info);
	const WORD old_color_attrs = buffer_info.wAttributes;

	// We need to flush the stream buffers into the console before each
	// SetConsoleTextAttribute call lest it affect the text that is already
	// printed but has not yet reached the console.
	fflush(stdout);
	SetConsoleTextAttribute(stdout_handle,
			GetColorAttribute(color) | FOREGROUND_INTENSITY);
	vprintf(fmt, args);

	fflush(stdout);
	// Restores the text color.
	SetConsoleTextAttribute(stdout_handle, old_color_attrs);
#else
	printf("\033[0;3%sm", GetAnsiColorCode(color));
	vprintf(fmt, args);
	printf("\033[m");  // Resets the terminal to default.
#endif  // GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
	va_end(args);
}

// Formats a countable noun.  Depending on its quantity, either the
// singular form or the plural form is used. e.g.
//
// FormatCountableNoun(1, "formula", "formuli") returns "1 formula".
// FormatCountableNoun(5, "book", "books") returns "5 books".
static std::string FormatCountableNoun(int count, const char * singular_form, const char * plural_form)
{
	return internal::StreamableToString(count) + " " + (count == 1 ? singular_form : plural_form);
}

// Formats the count of tests.
static std::string FormatTestCount(int test_count)
{
	return FormatCountableNoun(test_count, "test", "tests");
}

// Formats the count of test cases.
static std::string FormatTestCaseCount(int test_case_count)
{
	return FormatCountableNoun(test_case_count, "test case", "test cases");
}

// Converts a TestPartResult::Type enum to human-friendly string
// representation.  Both kNonFatalFailure and kFatalFailure are translated
// to "Failure", as the user usually doesn't care about the difference
// between the two when viewing the test result.
static const char * TestPartResultTypeToString(TestPartResult::Type type)
{
	switch (type)
	{
	case TestPartResult::kSuccess:
		return "Success";

	case TestPartResult::kNonFatalFailure:
	case TestPartResult::kFatalFailure:
#ifdef _MSC_VER
		return "error: ";
#else
		return "Failure\n";
#endif
	default:
		return "Unknown result type";
	}
}

static const char kTypeParamLabel[] = "TypeParam";
static const char kValueParamLabel[] = "GetParam()";

void PrintFullTestCommentIfPresent(const TestInfo& test_info)
{
	const char* const type_param = test_info.type_param();
	const char* const value_param = test_info.value_param();

	if (type_param != NULL || value_param != NULL)
	{
		printf(", where ");
		if (type_param != NULL)
		{
			printf("%s = %s", kTypeParamLabel, type_param);
			if (value_param != NULL)
				printf(" and ");
		}
		if (value_param != NULL)
		{
			printf("%s = %s", kValueParamLabel, value_param);
		}
	}
}

// Prints a TestPartResult to an std::string.
static std::string PrintTestPartResultToString(const TestPartResult& test_part_result)
{
	return (Message() << internal::FormatFileLocation(test_part_result.file_name(), test_part_result.line_number()) << " " << TestPartResultTypeToString(test_part_result.type()) << test_part_result.message()).GetString();
}

// Prints a TestPartResult.
static void PrintTestPartResult(const TestPartResult& test_part_result)
{
	const std::string& result = PrintTestPartResultToString(test_part_result);
	printf("%s\n", result.c_str());
	fflush(stdout);
	// If the test program runs in Visual Studio or a debugger, the
	// following statements add the test part result message to the Output
	// window such that the user can double-click on it to jump to the
	// corresponding source code location; otherwise they do nothing.
#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
	// We don't call OutputDebugString*() on Windows Mobile, as printing
	// to stdout is done by OutputDebugString() there already - we don't
	// want the same message printed twice.
	::OutputDebugStringA(result.c_str());
	::OutputDebugStringA("\n");
#endif
}

// Fired before each iteration of tests starts.
void MPIPrettyUnitTestResultPrinter::OnTestIterationStart(const UnitTest& unit_test, int iteration)
{
	if (_mpi_rank != 0)
		return;
	if (GTEST_FLAG(repeat) != 1)
		printf("\nRepeating all tests (iteration %d) . . .\n\n", iteration + 1);

	const char* const filter = GTEST_FLAG(filter).c_str();
	const std::string str_filter(filter);

	// Prints the filter if it's not *.  This reminds the user that some
	// tests may be skipped.
	if (!str_filter.empty())
	{
//  if (!String::CStringEquals(filter, kUniversalFilter)) {
		ColoredPrintf(COLOR_YELLOW, "Note: %s filter = %s\n", GTEST_NAME_, filter);
	}

#if 0
	if (internal::ShouldShard(kTestTotalShards, kTestShardIndex, false))
	{
		const Int32 shard_index = Int32FromEnvOrDie(kTestShardIndex, -1);
		ColoredPrintf(COLOR_YELLOW,
				"Note: This is test shard %d of %s.\n",
				static_cast<int>(shard_index) + 1,
				internal::posix::GetEnv(kTestTotalShards));
	}
#endif

	if (GTEST_FLAG(shuffle))
	{
		ColoredPrintf(COLOR_YELLOW, "Note: Randomizing tests' orders with a seed of %d .\n", unit_test.random_seed());
	}

	ColoredPrintf(COLOR_GREEN, "[==========] ");
	printf("Running %s from %s.\n", FormatTestCount(unit_test.test_to_run_count()).c_str(), FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnEnvironmentsSetUpStart(const UnitTest& /*unit_test*/)
{
	if (_mpi_rank != 0)
		return;
	ColoredPrintf(COLOR_GREEN, "[----------] ");
	printf("Global test environment set-up.\n");
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnTestCaseStart(const TestCase& test_case)
{
	if (_mpi_rank != 0)
		return;
	const std::string counts = FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
	ColoredPrintf(COLOR_GREEN, "[----------] ");
	printf("%s from %s", counts.c_str(), test_case.name());
	if (test_case.type_param() == NULL)
	{
		printf("\n");
	}
	else
	{
		printf(", where %s = %s\n", kTypeParamLabel, test_case.type_param());
	}
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnTestStart(const TestInfo& test_info)
{
	if (_mpi_rank != 0)
		return;
	ColoredPrintf(COLOR_GREEN, "[ RUN      ] ");
	PrintTestName(test_info.test_case_name(), test_info.name());
	printf("\n");
	fflush(stdout);
}

// Called after an assertion failure.
void MPIPrettyUnitTestResultPrinter::OnTestPartResult(const TestPartResult& result)
{
	if (_mpi_rank != 0)
		return;
	// If the test part succeeded, we don't need to do anything.
	if (result.type() == TestPartResult::kSuccess)
		return;

	// Print failure message from the assertion (e.g. expected this and got that).
	PrintTestPartResult(result);
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnTestEnd(const TestInfo& test_info)
{
	if (_mpi_rank != 0)
		return;
	if (test_info.result()->Passed())
	{
		ColoredPrintf(COLOR_GREEN, "[       OK ] ");
	}
	else
	{
		ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
	}
	PrintTestName(test_info.test_case_name(), test_info.name());
	if (test_info.result()->Failed())
		PrintFullTestCommentIfPresent(test_info);

	if (GTEST_FLAG(print_time))
	{
		printf(" (%s ms)\n", internal::StreamableToString(test_info.result()->elapsed_time()).c_str());
	}
	else
	{
		printf("\n");
	}
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnTestCaseEnd(const TestCase& test_case)
{
	if (!GTEST_FLAG(print_time))
		return;
	if (_mpi_rank != 0)
		return;

	const std::string counts = FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
	ColoredPrintf(COLOR_GREEN, "[----------] ");
	printf("%s from %s (%s ms total)\n\n", counts.c_str(), test_case.name(), internal::StreamableToString(test_case.elapsed_time()).c_str());
	fflush(stdout);
}

void MPIPrettyUnitTestResultPrinter::OnEnvironmentsTearDownStart(const UnitTest& /*unit_test*/)
{
	if (_mpi_rank != 0)
		return;
	ColoredPrintf(COLOR_GREEN, "[----------] ");
	printf("Global test environment tear-down\n");
	fflush(stdout);
}

// Internal helper for printing the list of failed tests.
void MPIPrettyUnitTestResultPrinter::PrintFailedTests(const UnitTest& unit_test)
{
	const int failed_test_count = unit_test.failed_test_count();
	if (failed_test_count == 0)
	{
		return;
	}

	for (int i = 0; i < unit_test.total_test_case_count(); ++i)
	{
		const TestCase& test_case = *unit_test.GetTestCase(i);
		if (!test_case.should_run() || (test_case.failed_test_count() == 0))
		{
			continue;
		}
		for (int j = 0; j < test_case.total_test_count(); ++j)
		{
			const TestInfo& test_info = *test_case.GetTestInfo(j);
			if (!test_info.should_run() || test_info.result()->Passed())
			{
				continue;
			}
			ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
			printf("%s.%s", test_case.name(), test_info.name());
			PrintFullTestCommentIfPresent(test_info);
			printf("\n");
		}
	}
}

void MPIPrettyUnitTestResultPrinter::OnTestIterationEnd(const UnitTest& unit_test, int /*iteration*/)
{
	if (_mpi_rank != 0)
		return;
	ColoredPrintf(COLOR_GREEN, "[==========] ");
	printf("%s from %s ran.", FormatTestCount(unit_test.test_to_run_count()).c_str(), FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
	if (GTEST_FLAG(print_time))
	{
		printf(" (%s ms total)", internal::StreamableToString(unit_test.elapsed_time()).c_str());
	}
	printf("\n");
	ColoredPrintf(COLOR_GREEN, "[  PASSED  ] ");
	printf("%s.\n", FormatTestCount(unit_test.successful_test_count()).c_str());

	int num_failures = unit_test.failed_test_count();
	if (!unit_test.Passed())
	{
		const int failed_test_count = unit_test.failed_test_count();
		ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
		printf("%s, listed below:\n", FormatTestCount(failed_test_count).c_str());
		PrintFailedTests(unit_test);
		printf("\n%2d FAILED %s\n", num_failures, num_failures == 1 ? "TEST" : "TESTS");
	}

	int num_disabled = unit_test.reportable_disabled_test_count();
	if (num_disabled && !GTEST_FLAG(also_run_disabled_tests))
	{
		if (!num_failures)
		{
			printf("\n");  // Add a spacer if no FAILURE banner is displayed.
		}
		ColoredPrintf(COLOR_YELLOW, "  YOU HAVE %d DISABLED %s\n\n", num_disabled, num_disabled == 1 ? "TEST" : "TESTS");
	}
	// Ensure that Google Test output is printed before, e.g., heapchecker output.
	fflush(stdout);
}

}

