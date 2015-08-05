/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-16
 * \brief  Implementation of string helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StringTools.h"

#include <algorithm>
#include <cctype>
#include <iomanip>

#include "logog/include/logog.hpp"

#include <boost/property_tree/json_parser.hpp>


namespace BaseLib
{
std::list<std::string> splitString(const std::string &str, char delim)
{
	std::list<std::string> strList;
	std::stringstream ss(str);
	std::string item;
	while(getline(ss, item, delim))
		strList.push_back(item);
	return strList;
}

std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace)
{
	std::string::size_type pos = stringToReplace.find(searchString, 0);
	std::string::size_type const intLengthSearch = searchString.length();

	while (std::string::npos != pos) {
		stringToReplace.replace(pos, intLengthSearch, replaceString);
		pos = stringToReplace.find(searchString, 0);
	}
	return stringToReplace;
}

void trim(std::string &str, char ch)
{
	std::string::size_type pos = str.find_last_not_of(ch);
	if(pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(ch);
		if(pos != std::string::npos)
			str.erase(0, pos);
	}
	else
		str.erase(str.begin(), str.end());
}

void simplify(std::string &str)
{
	trim (str);
	str.erase(
	    std::unique(str.begin(), str.end(), [](char a, char b) { return a == ' ' && b == ' '; }),
		str.end()
	);
}

std::string padLeft(std::string const& str, int maxlen, char ch)
{
	std::stringstream ss(str);
	ss << std::right << std::setw(maxlen) << std::setfill(ch) << str;
	return ss.str();
}

std::string propertyTreeToString(boost::property_tree::ptree const& tree)
{
	std::ostringstream ss;
	boost::property_tree::write_json(ss, tree);
	return ss.str();
}


std::string leftPadding(const std::string &src, size_t total_len, const char paddingChar)
{
    std::string str(src);
    if(total_len > str.size())
        str.insert(0, total_len - str.size(), paddingChar);

    return str;
}

std::string rightPadding(const std::string &src, size_t total_len, const char paddingChar)
{
    std::string str(src);
    if(total_len > str.size())
        str.insert(str.size(), total_len - str.size(), paddingChar);

    return str;
}

std::string bothPadding(const std::string &src, size_t total_len, const char paddingChar)
{
    std::string str(src);
    if(total_len == str.size())
        return str;

    size_t pad_len = total_len - str.size();
    str = leftPadding(str, pad_len / 2 + str.size(), paddingChar);
    str = rightPadding(str, total_len, paddingChar);

    return str;
}


} // end namespace BaseLib

#ifdef MSVC
void correctScientificNotation(std::string filename, size_t precision)
{
	std::ifstream stream;
	std::ofstream outputStream;

	stream.open(filename.c_str());
	std::string tmpFilename = filename + ".tmp";
	outputStream.open(tmpFilename.c_str());

	if (!stream)
	{
		WARN("correctScientificNotation: fstream is not open.");
		return;
	}

	std::string line;

	// Iterate over lines in stream
	while (getline(stream, line))
	{
		std::string word;
		std::istringstream iss(line);
		// Iterate over all words in line
		while (iss >> word)
		{
			// Search for e+0
			std::size_t exponentPosition = word.find("e+0", precision);
			if (exponentPosition == std::string::npos)
				// If not found search for e-0
				exponentPosition = word.find("e-0", precision);
			if (exponentPosition != std::string::npos)
			{
				std::size_t wordSize = word.size();
				std::size_t exponentSize = wordSize - exponentPosition;

				if(exponentSize > 4)
				{
					// Erase the leading zero considering trailing characters
					int i = wordSize - 1;
					while (!isdigit(word[i]))
						--i;

					size_t erasePos = wordSize - 3 - (wordSize - 1 - i);
					std::string eraseString = word.substr(erasePos, 1);
					if (eraseString.find("0") != std::string::npos)
						word.erase(erasePos, 1);
				}
			}

			outputStream << word << " ";
		}
		outputStream << "\n";
	}

	stream.close();
	outputStream.close();

	remove(filename.c_str());
	rename(tmpFilename.c_str(), filename.c_str());
}
#endif

