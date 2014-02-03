// Copyright (C) 2014 Alexey Anikeenko
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CFGLINESTREAM_H
#define CFGLINESTREAM_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <exception>

// Simple configuration file parser.
// Reads non-comment lines one-by-one converting them into istringstreams
// Example of usage:
//   Cfglinestream cfg(filename);
//   cfg >> var1 >> var2; // reads two values from the first line
//   cfg >> output_stream; // reads file name from next line and opens stream
class Cfglinestream 
{
public:
    explicit Cfglinestream(const char *filename);
    template<class T> std::istream& operator >> (T & t);

protected:
    std::ifstream m_in;
    std::istringstream m_out;
    void next_line();
};

// Initializes input stream
Cfglinestream::Cfglinestream(const char *filename) : m_in(filename)
{
    if (!m_in) {
        std::cerr << "Can't read config file : " << filename << std::endl;
        throw std::exception();
    }
}

// Saves next non-comment line into internal istringstream object
// skipping leading whitespaces.
void Cfglinestream::next_line()
{
    const std::string delims(" \t");

    std::string str;
    std::string::size_type i;
    bool get_next_line = true;

    while (get_next_line) {
        str.clear();
        if (std::getline(m_in, str)) {
            i = str.find_first_not_of(delims);
            if (i != std::string::npos && str[i] != '#')
                get_next_line = false;
        } else {
            get_next_line = false;  // end of file reached
        }
    }

    m_out.clear();
    m_out.str(str);
}

// Gets next line, reads first value from it, returns the remaining stringstream
template<class T> std::istream& Cfglinestream::operator >> (T & t)
{
    next_line();
    return (m_out >> t);
}

// Specialization of >> operator for strings in order to
// handle quoted values with spaces inside.
template<> std::istream& Cfglinestream::operator >> (std::string & s)
{
    const std::string delims(" \t");

    std::string::size_type i, j;
    next_line();

    i = m_out.str().find_first_not_of(delims);
    if (i != std::string::npos && m_out.str()[i] == '"') {
        j = m_out.str().find('"', i + 1);
        m_out.str(m_out.str().substr(i + 1, j - i - 1));
        s = m_out.str();
    } else {
        m_out >> s;
    }

    return m_out;
}

// Specialization of >> operator for ofstream type.
// Reads filename and opens file for writing if name doesn't start with '!'
template<> std::istream& Cfglinestream::operator >> (std::ofstream & s)
{
    std::string fname;
    this->operator >> (fname);

    if (!fname.empty() && fname[0] != '!') {
        s.open(fname.c_str());

        if (!s)
            std::cerr << "Error opening file for writing : " << fname << std::endl;
    }

    return m_out;
}

// Specialization of >> operator for ifstream type.
// Reads filename and opens file for reading if name doesn't start with '!'
template<> std::istream& Cfglinestream::operator >> (std::ifstream & s)
{
    std::string fname;
    this->operator >> (fname);

    if (!fname.empty() && fname[0] != '!') {
        s.open(fname.c_str());

        if (!s)
            std::cerr << "Error opening file for reading : " << fname << std::endl;
    }

    return m_out;
}

#endif /* CFGLINESTREAM_H */
