// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

// Create and change directories
// Author: Shayan Hoshyari

#include "foldertools.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#ifdef _WIN32
#include "Windows.h"
#else // assuming mac or linux
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#define MAX_STRING_LEN 4048
#define MAX_STACK_DEPTH 20

namespace foldertools {

#if _WIN32

namespace {
bool makesingledir(const char *path) {
  int err = CreateDirectoryA(path, NULL);
  return err != 0;
}
} // namespace

bool setdir(const char *name) {
  SetCurrentDirectoryA(name);
  return true;
}

void cwd(char ans[], const int max_size) {
  GetCurrentDirectoryA(max_size, ans);
}

bool exists(const std::string iname) {
  // Check file
  {
    std::ifstream ifstr(iname);
    if (ifstr.is_open())
      return true;
  }
  DWORD dwAttrib = GetFileAttributes(iname.c_str());
  return (dwAttrib != INVALID_FILE_ATTRIBUTES) &&
         (dwAttrib & FILE_ATTRIBUTE_DIRECTORY);
}

#else /* mac and linux */

namespace {
bool makesingledir(const char *name) {
  char tmp[MAX_STRING_LEN];
  char *p = NULL;
  size_t len;

  snprintf(tmp, sizeof(tmp), "%s", name);
  len = strlen(tmp);
  if (tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for (p = tmp + 1; *p; p++)
    if (*p == '/') {
      *p = 0;
      mkdir(tmp, S_IRWXU);
      *p = '/';
    }
  mkdir(tmp, S_IRWXU);

  return true;
}
} // namespace

bool setdir(const char *name) {
  int success = chdir(name);
  return (success == 0);
}

void cwd(char ans[], const int max_size) {
  char *success = getcwd(ans, max_size);
  (void)(success); // unused
}

bool exists(const std::string iname) {
  // Check file
  {
    std::ifstream ifstr(iname);
    if (ifstr.is_open())
      return true;
  }
  // Check folder
  DIR *dir = opendir("mydir");
  if (dir) {
    closedir(dir);
    return true;
  } else {
    return false;
  }
}

#endif // UNIX

namespace {
bool is_valid_folder(std::string in) {

  if (in == ".")
    return false;
  if (in == "")
    return false;
  if (in == "\\")
    return false;
  if (in == "/")
    return false;
  if (in == ".\\")
    return false;
  if (in == "./")
    return false;
  if (in == "")
    return false;

  return true;
}

static char stack[MAX_STACK_DEPTH][MAX_STRING_LEN];
static int stack_pos = 0;
} // namespace

// A fixed version of the buggy code in.
// https://stackoverflow.com/questions/1530760/how-do-i-recursively-create-a-folder-in-win32
bool makedir(const char *name) {
  std::string path(name);
  std::string::size_type pos = 0;
  std::string::size_type pos_prev = 0;
  while (true) {
    pos_prev = pos;

    pos = path.find_first_of("\\/", pos + 1);

    std::string current = path.substr(0, pos);
    // trim initial '/' and '\\'.
    while (true) {
      if (current.empty())
        break;
      if (current[0] != '/' && current[0] != '\\')
        break;
      ++pos;
    }

    // No slashes,
    if (pos == std::string::npos) {
      // valid dir, create and then get out.
      if (is_valid_folder(current)) {
        makesingledir(current.c_str());
        break;
      }
      // Invalid nonesense, get out
      else {
        break;
      }
    }
    // Found a slash
    else {
      // valid dir, create
      if (is_valid_folder(current)) {
        makesingledir(current.c_str());
      }
    }
  }

  return true;
}

bool makeandsetdir(const char *name) {
  makedir(name);
  setdir(name);
  return true;
}

void pushd() {

  assert(stack_pos < MAX_STACK_DEPTH - 1);
  cwd(stack[stack_pos], MAX_STRING_LEN);
  ++stack_pos;
}

void popd() {
  if (stack_pos > 0) {
    --stack_pos;
    setdir(stack[stack_pos]);
  }
}

std::string cwd() {
  char ans[MAX_STRING_LEN];
  cwd(ans, MAX_STRING_LEN);
  return std::string(ans);
}

} // namespace foldertools
