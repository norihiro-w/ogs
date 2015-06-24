/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "Ogs5FileTools.h"

#include "makros.h"

std::string GetUncommentedLine(std::string& line)
{
    std::string zeile = "";
    int i = 0, j = 0;
    //----------------------------------------------------------------------
    i = (int) line.find_first_not_of(" ",0);
    j = (int) line.find(";",i);
    if((i != -1))
        zeile = line.substr(i,j - i);
    i = (int) zeile.find_last_not_of(" ");
    if(i >= 0)
    {
        line = zeile.substr(0,i + 1);
        zeile = line;
    }

    return zeile;
}

std::ios::pos_type GetNextSubKeyword(std::ifstream* file,std::string* line, bool* keyword)
{
    char buffer[MAX_ZEILE];
    std::ios::pos_type position = file->tellg();
    *keyword = false;
    std::string line_complete;
    int i,j;
    // Look for next subkeyword
    while(!(line_complete.find("$") != std::string::npos) && (!file->eof()))
    {
        file->getline(buffer,MAX_ZEILE);
        line_complete = buffer;
        if(line_complete.find("#") != std::string::npos)
        {
            *keyword = true;
            return position;
        }
        i = (int) line_complete.find_first_not_of(" ",0);
        j = (int) line_complete.find(";",i);
        if(j < 0)
            j = (int)line_complete.length();
        //if(j!=i) break;
        if(i != -1)
            *line = line_complete.substr(i,j - i);
    }
    return position;
}

/**
 * read a non blank line from given input stream
 * @param in the input stream
 * @return read line into a string
 */
std::string readNonBlankLineFromInputStream(std::istream & in)
{
    std::string line;

    bool not_finished (true);
    while (not_finished)
    {
        // read line
        getline(in, line);
        if (!in.fail())
        {
            // skip initial space characters
            std::string::size_type i (line.find_first_not_of(" ", 0));
            // search comment symbol ;
            std::string::size_type j (line.find(";", i));
            if (j == i) // first non space character is equal to the comment symbol
                not_finished = true;
            else
            {
                if ((i != std::string::npos))
                    // cut string from first non blank space character until the first comment character
                    line = line.substr(i, j - i);

                // remove last blank spaces
                i = line.find_last_not_of(" ");
                if (i != std::string::npos)
                    line = line.substr(0, i + 1);

                not_finished = false;
            }
        }
        else
        {
            line = "";
            not_finished = false;
        }
    }
    return line;
}

