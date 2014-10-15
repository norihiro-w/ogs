/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MPIPRETTYUNITTESTRESULTPRINTER_H_
#define MPIPRETTYUNITTESTRESULTPRINTER_H_

#include "gtest/gtest.h"

namespace testing
{

/// gtest printer for MPI computing
class MPIPrettyUnitTestResultPrinter: public TestEventListener
{
public:
	explicit MPIPrettyUnitTestResultPrinter(int mrank)
	: _mpi_rank(mrank)
	{}

	static void PrintTestName(const char * test_case, const char * test)
	{
		printf("%s.%s", test_case, test);
	}

	// The following methods override what's in the TestEventListener class.
	virtual void OnTestProgramStart(const UnitTest& /*unit_test*/) {}
	virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
	virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
	virtual void OnEnvironmentsSetUpEnd(const UnitTest& /*unit_test*/) {}
	virtual void OnTestCaseStart(const TestCase& test_case);
	virtual void OnTestStart(const TestInfo& test_info);
	virtual void OnTestPartResult(const TestPartResult& result);
	virtual void OnTestEnd(const TestInfo& test_info);
	virtual void OnTestCaseEnd(const TestCase& test_case);
	virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
	virtual void OnEnvironmentsTearDownEnd(const UnitTest& /*unit_test*/) {}
	virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
	virtual void OnTestProgramEnd(const UnitTest& /*unit_test*/) {}

private:
	static void PrintFailedTests(const UnitTest& unit_test);

	const int _mpi_rank;
};

}

#endif /* MPIPRETTYUNITTESTRESULTPRINTER_H_ */
