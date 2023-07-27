// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

// Create and change directories
// None of these are thread safe.
// Author: Shayan Hoshyari

#include <string>

namespace foldertools
{

// Create a folder
bool makedir(const char * name);

// Go to a folder
bool setdir(const char * name);

// Create a folder and then go to it.
bool makeandsetdir(const char * name);

// Equivalent of the command pushd in bash.
// Remember the current directory in the directory stack
void pushd();

// Equivalent of the command pwd in bash.
// Goes back to thre previous directory in directoy stack
void popd();

// get current working directory
// the int is the maximum size of the array in the second version
std::string cwd();
void cwd(char[], const int); 

// Delete a folder if it exists.
bool exists(const std::string);

} // end of foldertools

